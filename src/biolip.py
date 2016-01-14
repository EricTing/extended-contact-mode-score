#!/usr/bin/env python


import os
import pybel
import shlex
import subprocess32
import tempfile

from apoc_inputs import ApocInput


class FastSearch:
    def __init__(self, lig_path,
                 index="/ddnB/work/jaydy/dat/BioLip/ligand_nr.fs",
                 tani=0.5, minimum_size=6):
        self.lig_path = lig_path
        self.index = index
        self.tani = tani
        self.minimum_size = minimum_size

    def search(self, verbose=True):
        try:
            ofn = tempfile.mkstemp(suffix='.sdf')[1]
            cmd = "babel %s %s -s %s -at%f" % (self.index,
                                               ofn,
                                               self.lig_path,
                                               self.tani)
            if verbose:
                print cmd
            args = shlex.split(cmd)
            subprocess32.call(args)
            mols = [mol for mol in pybel.readfile(format="sdf", filename=ofn)
                    if len(mol.atoms) >= self.minimum_size]
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
        self.apoc_input = ApocInput(self.lig_path,
                                    self.prt_path,
                                    threshold=7.0).input4Apoc()

    def calculate(self):
        self.pkt_path = tempfile.mkstemp()[1]
        with open(self.pkt_path, 'w') as ofs:
            ofs.write(self.apoc_input)

        ref_pkt_path = tempfile.mkstemp()[1]
        query = BioLipQuery(self.lig_path,
                            index="/ddnB/work/jaydy/dat/BioLip/ligand_nr.fs")
        for ref_lig in query.search():
            ref_prt_pdb = query.correspondingPrtPdb(ref_lig)
            if os.path.exists(ref_prt_pdb):
                ref_apoc_input = ApocInput(ref_lig,
                                           ref_prt_pdb,
                                           threshold=7.0).input4Apoc()
                with open(ref_pkt_path, 'w') as ofs:
                    ofs.write(ref_apoc_input)

                cmd = "apoc -plen 1 %s %s" % (self.pkt_path, ref_pkt_path)
                print cmd
                args = shlex.split(cmd)
                apoc_result = subprocess32.check_output(args)
                print apoc_result

        os.remove(self.pkt_path)
        os.remove(ref_pkt_path)


def test():
    cms = BioLipReferencedSpearmanR("/work/jaydy/dat/BioLip/sml_ligand_nr/03/103m_NBN_A_1.pdb",
                                    "/ddnB/work/jaydy/dat/BioLip/prt/03/103mA.pdb")
    cms.calculate()


if __name__ == '__main__':
    test()
    pass
