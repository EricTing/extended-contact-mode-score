#!/usr/bin/env python


import os
import pybel
import luigi
import json
from collections import defaultdict
from run_control_apoc import LpcApocResultTask, LpcKcombuResult
from run_control_apoc import LigandExpStructureInMol2
from run_control_apoc import LpcPocketPathTask
from run_control_apoc import ApocResultParer, PkcombuAtomMatchParser
from single_complex_cms import contactVector
from scipy.stats import spearmanr


class Path(luigi.Task):
    tname = luigi.Parameter()
    qname = luigi.Parameter(default="")
    subset = luigi.Parameter(default="")

    def mid_two(self):
        return self.tname[1:3]

    def dirname(self):
        path = os.path.join("/ddnB/work/jaydy/working/apoc_spearmanr",
                            self.subset,
                            self.mid_two())
        try:
            os.makedirs(path)
        except:
            pass
        return path


class ApocSpearmanR(Path):
    def requires(self):
        return [
            LpcApocResultTask(self.tname, self.qname, self.subset),
            LpcKcombuResult(self.tname, self.qname, self.subset),
        ]

    def output(self):
        path = os.path.join(
            self.dirname(),
            "%s__%s.json" % (self.tname, self.qname)
        )
        return luigi.LocalTarget(path)

    def run(self):

        class Complex:
            def __init__(self, name):
                self.name = name

            def readProteinCoords(self, chainid, resids):
                path = LpcPocketPathTask(self.name).output().path
                coords = defaultdict(list)
                with open(path, 'r') as ifs:
                    lines = ifs.readlines()
                for idx in resids:
                    for line in lines:
                        tokens = line.split()
                        if (len(tokens) > 11
                            and (tokens[4] == chainid
                                 and int(tokens[5]) == idx)):
                            coords[idx].append(map(float, tokens[6:9]))
                        if line.startswith("PKT"):
                            break
                return coords

            def readLigandCoords(self, name, mlist):
                path = LigandExpStructureInMol2(name).output().path
                mol = pybel.readfile("mol2", path).next()
                mol.removeh()
                all_coords = [atom.coords for atom in mol
                              if not atom.OBAtom.IsHydrogen()]
                coords = []
                for idx in mlist:
                    coords.append(all_coords[idx - 1])

                return coords

        def addVirtualPoints(pocket_alignment, t_coords, q_coords):
            t_virtual_coords, q_virtual_coords = [], []
            virtual_point = [9999.0, 9999.0, 9999.0]
            for tidx, qidx in zip(pocket_alignment.template_res,
                                  pocket_alignment.query_res):
                tc = t_coords[tidx]
                qc = q_coords[qidx]
                t_virtual_coords.extend(tc)
                q_virtual_coords.extend(qc)
                lt = len(tc)
                lq = len(qc)
                if lt > lq:
                    q_virtual_coords.extend([virtual_point
                                             for _ in range(lt - lq)])
                if lq > lt:
                    t_virtual_coords.extend([virtual_point
                                             for _ in range(lq - lt)])
            return t_virtual_coords, q_virtual_coords

        apoc_task = LpcApocResultTask(self.tname, self.qname, self.subset)
        with apoc_task.output().open('r') as ifs:
            apoc_parser = ApocResultParer(ifs.read())
        pocket_alignment = apoc_parser.queryPocket(self.tname, self.qname)

        kcombu_task = LpcKcombuResult(self.tname, self.qname, self.subset)
        kcombu_parser = PkcombuAtomMatchParser(kcombu_task.output().path)
        t_lig_list, q_lig_list = kcombu_parser.getMatchingSerialNums()

        t_complex = Complex(self.tname)
        t_prt_coords = t_complex.readProteinCoords(pocket_alignment.template_chainid,
                                                   pocket_alignment.template_res)
        t_lig_coords = t_complex.readLigandCoords(self.tname, t_lig_list)

        q_complex = Complex(self.qname)
        q_prt_coords = q_complex.readProteinCoords(pocket_alignment.query_chainid,
                                                   pocket_alignment.query_res)
        q_lig_coords = q_complex.readLigandCoords(self.qname, q_lig_list)

        t_prt_coords, q_prt_coords = addVirtualPoints(pocket_alignment,
                                                      t_prt_coords,
                                                      q_prt_coords)
        assert len(q_prt_coords) == len(t_prt_coords)

        t_vec = contactVector(t_prt_coords, t_lig_coords)
        q_vec = contactVector(q_prt_coords, q_lig_coords)

        corr = spearmanr(t_vec, q_vec)
        sparman_corr = corr.correlation
        pval = corr.pvalue
        tc = kcombu_parser.getTc()
        ps_score = pocket_alignment.ps_score

        data = {
            "spearmanr": sparman_corr,
            "pval": pval,
            "tc": tc,
            "ps_score": ps_score,
            "lig_size": len(t_lig_coords),
            "prt_size": len(t_prt_coords),
        }
        to_write = json.dumps(data, sort_keys=True,
                              indent=4, separators=(',', ': '))
        with open(self.output().path, 'w') as ofs:
            ofs.write(to_write)


def main(tname, qname):
    targets = [
        # ApocSpearmanR("3vn9_ANK_A_401", "4a2a_ATP_B_1391", "subject"),
        # ApocSpearmanR("4a9w_FAD_B_1353", "4dql_FAD_B_1101", "subject"),
        ApocSpearmanR(tname, qname, "subject"),
    ]
    luigi.build(targets, local_scheduler=True)
    pass


if __name__ == '__main__':
    import sys
    if sys.argv[1] != "tname":
        main(sys.argv[1], sys.argv[2])
