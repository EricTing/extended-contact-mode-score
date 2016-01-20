#!/usr/bin/env python

import luigi
import json
import os
from biolip import BioLipReferencedSpearmanR

OUTPUT_DIR = "/ddnB/work/jaydy/working/biolip/"


class Path(luigi.Task):
    lig_pdb = luigi.Parameter()
    prt_dir = luigi.Parameter(default="/work/jaydy/dat/BioLip/prt/")
    lig_dir = luigi.Parameter(default="/ddnB/work/jaydy/dat/BioLip/ligand_nr/")

    def prtPdb(self):
        code = os.path.basename(self.lig_pdb)
        mid_two = code[1:3]
        pdb = code.split('_')[0] + code.split('_')[2] + ".pdb"
        path = os.path.join(self.prt_dir, mid_two, pdb)
        return path

    def ligPdb(self):
        mid_two = self.lig_pdb[1:3]
        path = os.path.join(self.lig_dir, mid_two, self.lig_pdb)
        return path

    def mid_two(self):
        return self.lig_pdb[1:3]


class BioLipBioLip(Path):

    def output(self):
        path = os.path.join(OUTPUT_DIR, self.mid_two(), self.lig_pdb + ".json")
        return luigi.LocalTarget(path)

    def run(self):
        prt_pdb = self.prtPdb()
        lig_pdb = self.ligPdb()
        biolip_spearmanr = BioLipReferencedSpearmanR(lig_pdb, prt_pdb)
        inf = float('inf')
        result = biolip_spearmanr.calculate(maximum_search_results=inf,
                                            max_tani=0.9)

        to_write = json.dumps(result, sort_keys=True,
                              indent=4, separators=(',', ': '))
        ofn = self.output().path
        print(ofn)
        with open(ofn, 'w') as ofs:
            ofs.write(to_write)

    def requires(self):
        pass


def main(name):
    luigi.build([
        BioLipBioLip(name)
    ],
                local_scheduler=True
    )

if __name__ == '__main__':
    import sys
    name = sys.argv[1]
    main(name)
