#!/usr/bin/env python

import os
import pybel
import shlex
import subprocess32
import tempfile
import numpy as np
from collections import defaultdict

from astex_xcms import runPkcombu
from apoc_inputs import ApocInput
from run_control_apoc import ApocResultParer, PkcombuAtomMatchParser
from single_complex_cms import contactVector
from scipy.stats import spearmanr


def readLigandCoords(path, mlist, format="sdf"):
    mol = pybel.readfile(format, path).next()
    mol.removeh()
    all_coords = [atom.coords for atom in mol if not atom.OBAtom.IsHydrogen()]
    coords = []
    for idx in mlist:
        coords.append(all_coords[idx - 1])

    return coords


class FastSearch:
    def __init__(self,
                 lig_path,
                 index="/ddnB/work/jaydy/dat/BioLip/ligand_nr.fs",
                 tani=0.5,
                 minimum_size=6,
                 max_tani=1.0,
                 maximum_search_results=100):
        self.lig_path = lig_path
        self.index = index
        self.tani = tani
        self.max_tani = max_tani
        self.minimum_size = minimum_size
        self.maximum_search_results = maximum_search_results

    def search(self, verbose=True):
        try:
            ofn = tempfile.mkstemp(suffix='.sdf')[1]
            cmd = "babel %s %s -s %s -at%f,%f" % (
                self.index, ofn, self.lig_path, self.tani, self.max_tani)
            if verbose:
                print cmd
            args = shlex.split(cmd)
            subprocess32.call(args)
            mols = [mol
                    for mol in pybel.readfile(format="sdf",
                                              filename=ofn)
                    if len(mol.atoms) >= self.minimum_size]
            # mols were implicitly ranked by their similarity to the query ligand
            if len(mols) > self.maximum_search_results:
                mols = mols[:self.maximum_search_results]
            return mols
        except Exception as detail:
            print "FAIL Fastsearch for %s" % (self.lig_path)
            print detail
        finally:
            os.remove(ofn)


class BioLipQuery(FastSearch):
    def correspondingPrtPdb(self, lig):
        self.prt_dir = "/work/jaydy/dat/BioLip/prt/"
        mid_two, code = lig.title.split('/')
        pdb = code.split('_')[0] + code.split('_')[2] + ".pdb"
        path = os.path.join(self.prt_dir, mid_two, pdb)
        return path


class BioLipReferencedSpearmanR:
    def __init__(self, lig_path, prt_path):
        self.lig_path = lig_path
        self.prt_path = prt_path
        apoc_input = ApocInput(self.lig_path, self.prt_path, threshold=7.0)
        self.apoc_input = apoc_input.input4Apoc()
        self.pkt = apoc_input.pocketSection()
        suffix = lig_path.split('.')[-1]
        # only read the first molecule
        self.lig = pybel.readfile(suffix, self.lig_path).next()

    def alignProteins(self, pkt1, res1, pkt2, res2):
        def prt_coords(pkt):
            res_coords = {}
            for line in pkt.splitlines():
                if len(line) > 56:
                    fullname = line[12:16]
                    # get rid of whitespace in atom names
                    split_list = fullname.split()
                    name = ""
                    if len(split_list) != 1:
                        # atom name has internal spaces, e.g. " N B ", so
                        # we do not strip spaces
                        name = fullname
                    else:
                        # atom name is like " CA ", so we can strip spaces
                        name = split_list[0]
                    resseq = int(line[22:26].split()[0])  # sequence identifier
                    if name == "CA" or name == "CB":  # only count the backbone
                        residue_id = (name, resseq)
                        try:
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            res_coords[residue_id] = np.array((x, y, z), "f")
                        except Exception:
                            raise Exception(
                                "Invalid or missing coordinates:\n %s" % line)
            return res_coords

        def align(pc1, pc2, res1, res2):
            reorded_pc1, reorded_pc2 = [], []
            for r1, r2 in zip(res1, res2):
                if ("CA", r1) in pc1 and ("CA", r2) in pc2:
                    reorded_pc1.append(pc1[("CA", r1)])
                    reorded_pc2.append(pc2[("CA", r2)])
                if ("CB", r1) in pc1 and ("CB", r2) in pc2:
                    reorded_pc1.append(pc1[("CB", r1)])
                    reorded_pc2.append(pc2[("CB", r2)])
            return reorded_pc1, reorded_pc2

        pc1 = prt_coords(pkt1)
        pc2 = prt_coords(pkt2)
        pc1, pc2 = align(pc1, pc2, res1, res2)
        if len(pc1) != len(pc2):
            raise Exception("Unequal pocket residue number %d vs %d" %
                            (len(pc1), len(pc2)))
        else:
            return pc1, pc2

    def alignLigands(self, lig1, lig2):
        try:
            sdf1 = tempfile.mkstemp()[-1]
            sdf2 = tempfile.mkstemp()[-1]
            oam_path = tempfile.mkstemp()[-1]

            lig1.removeh()
            lig2.removeh()

            lig1.write("sdf", sdf1, overwrite=True)
            lig2.write("sdf", sdf2, overwrite=True)

            runPkcombu(sdf1, sdf2, oam_path)
            kcombu = PkcombuAtomMatchParser(oam_path)
            list1, list2 = kcombu.getMatchingSerialNums()
            c1 = readLigandCoords(sdf1, list1)
            c2 = readLigandCoords(sdf2, list2)
            return kcombu, c1, c2
        finally:
            os.remove(oam_path)
            os.remove(sdf1)
            os.remove(sdf2)

    def calculate(self,
                  maximum_search_results=100,
                  max_tani=1.0,
                  minimum_size=6):
        self.pkt_path = tempfile.mkstemp()[1]
        with open(self.pkt_path, 'w') as ofs:
            ofs.write(self.apoc_input)

        ref_pkt_path = tempfile.mkstemp()[1]
        query = BioLipQuery(self.lig_path,
                            index="/ddnB/work/jaydy/dat/BioLip/ligand_nr.fs",
                            max_tani=max_tani,
                            minimum_size=minimum_size,
                            maximum_search_results=maximum_search_results)

        results = defaultdict(dict)
        for ref_lig in query.search():
            ref_prt_pdb = query.correspondingPrtPdb(ref_lig)
            if os.path.exists(ref_prt_pdb):
                try:
                    apoc_input = ApocInput(ref_lig, ref_prt_pdb, threshold=7.0)
                    ref_apoc_input = apoc_input.input4Apoc()
                    ref_pkt = apoc_input.pocketSection()
                    with open(ref_pkt_path, 'w') as ofs:
                        ofs.write(ref_apoc_input)

                    cmd = "apoc -plen 1 %s %s" % (self.pkt_path, ref_pkt_path)
                    args = shlex.split(cmd)
                    apoc_result = subprocess32.check_output(args)
                    apoc_parser = ApocResultParer(apoc_result)
                    if apoc_parser.numPocketSections() == 1:
                        pocket_alignment = apoc_parser.queryPocket()
                        global_alignment = apoc_parser.queryGlobal()
                        self_res = pocket_alignment.template_res
                        ref_res = pocket_alignment.query_res

                        pc1, pc2 = self.alignProteins(self.pkt, self_res,
                                                      ref_pkt, ref_res)
                        kcombu, lc1, lc2 = self.alignLigands(self.lig, ref_lig)
                        template_lig_atn, query_lig_atn = kcombu.getMatchingSerialNums(
                        )

                        vec1 = contactVector(lc1, pc1)
                        vec2 = contactVector(lc2, pc2)
                        tc = kcombu.getTc()
                        rho, pval = spearmanr(vec1, vec2)

                        my_result = {
                            "Tc": tc,
                            "spearmanr": rho,
                            "ps_score_pval": pocket_alignment.p_value,
                            "num_binding_res": len(ref_res),
                            "num_lig_atoms": len(lc1),
                            "pval": pval,
                            "ps_score": pocket_alignment.ps_score,
                            "tc_times_ps": tc * pocket_alignment.ps_score,
                            "TM-score": global_alignment.tm_score,
                            "rmsd": global_alignment.rmsd,
                            "seq_identity": global_alignment.seq_identity,
                            "query_res": pocket_alignment.query_res,
                            "query_lig_atom_num": query_lig_atn,
                            "template_lig_atom_num": template_lig_atn,
                            "template_res": pocket_alignment.template_res,
                        }

                        results[ref_lig.title] = my_result
                except Exception as detail:
                    print detail

        results = sorted(results.items(), key=lambda d: d[1].get('pval', 1))

        os.remove(self.pkt_path)
        os.remove(ref_pkt_path)

        return results


class FixedPocketBioLipReferencedSpearmanR(BioLipReferencedSpearmanR):
    def __init__(self, lig_path, prt_path, native_lig_path):
        self.lig_path = lig_path
        self.prt_path = prt_path
        apoc_input = ApocInput(native_lig_path, self.prt_path, threshold=7.0)
        self.apoc_input = apoc_input.input4Apoc()
        self.pkt = apoc_input.pocketSection()
        suffix = lig_path.split('.')[-1]
        self.lig = pybel.readfile(suffix, self.lig_path).next()
        print("This FixedPocketBioLipReferencedSpearmanR")


def test():
    cms = BioLipReferencedSpearmanR(
        "/work/jaydy/dat/BioLip/sml_ligand_nr/03/103m_NBN_A_1.pdb",
        "/ddnB/work/jaydy/dat/BioLip/prt/03/103mA.pdb")
    cms.calculate()

    cms = BioLipReferencedSpearmanR(
        "/work/jaydy/dat/BioLip/sml_ligand_nr/08/208l_CYS_A_1.pdb",
        "/work/jaydy/dat/BioLip/prt/08/208lA.pdb")
    cms.calculate()

    cms = FixedPocketBioLipReferencedSpearmanR(
        "/ddnB/work/jaydy/working/vina_biolip/187l_PXY_A_1.pdb/187l_PXY_A_1.pdb.pdbqt.vina.pdbqt",
        "/ddnB/work/jaydy/dat/BioLip/prt/87/187lA.pdb",
        "/ddnB/work/jaydy/dat/BioLip/sml_ligand_nr/87/187l_PXY_A_1.pdb")


if __name__ == '__main__':
    test()
    pass
