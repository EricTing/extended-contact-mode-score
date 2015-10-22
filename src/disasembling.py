#!/usr/bin/env python

import luigi
import shlex
import os
import subprocess32

from run_control_apoc import PdbMol2
from urls import WORKING_DIR
import glob

# 5 fails from 200 in astex
# 2MIP 1A07 1PSO 1MCQ 1HEF, where ligands are not annotated as non-polymers

def checkAstexDisasemble():
    dirname = "/ddnB/work/jaydy/working/xcms/disasemble"
    mol2_fns = glob.glob(dirname + "/*/*/*.mol2")

    pdb_ids = set(map(lambda path: os.path.basename(path).split('_')[0],
                      mol2_fns))

    total_astex = set([_.rstrip().upper() for _ in file("../dat/astex.lst")])
    for pdb_id in total_astex - pdb_ids:
        print pdb_id


class Path(luigi.Task):

    pdb_id = luigi.Parameter()

    def mid_two(self):
        return self.pdb_id[1:3]

    def dirname(self):
        path = os.path.join(WORKING_DIR, 'disasemble',
                            self.mid_two(), self.pdb_id)
        try:
            os.makedirs(path)
        except:
            pass
        return path


class CheckAstexDisasemble(Path):

    sdf_id = luigi.Parameter()

    def astex_dir(self):
        path = os.path.join("/ddnB/work/jaydy/dat/astex",
                            self.sdf_id)
        return path

    def astex_sdf(self):
        return os.path.join(self.astex_dir(),
                            self.sdf_id + ".sdf")

    def run(self):
        disasembled_mol2_fns = glob.glob(self.dirname() + "/*mol2")
        for disasembled_mol2 in disasembled_mol2_fns:
            oam_fn = os.path.splitext(disasembled_mol2)[0] + ".oam"
            cmd = "pkcombu -A %s -B %s -oam %s" % (self.astex_sdf(),
                                                   disasembled_mol2,
                                                   oam_fn)
            print cmd
            args = shlex.split(cmd)
            subprocess32.call(args)



class DisasemblePdb(Path):

    def requires(self):
        return PdbMol2(self.pdb_id)

    def output(self):
        pass

    def run(self):
        cif_ifn = self.requires().output().path
        cmd = "mmcif_disasembler --ifn %s --dir %s" % (cif_ifn, self.dirname())
        print cmd
        args = shlex.split(cmd)
        subprocess32.call(args)


def main(pdb_id):
    luigi.build([DisasemblePdb(pdb_id)],
                local_scheduler=True)


def check(sdf_id):
    pdb_id = sdf_id[:4]
    luigi.build([CheckAstexDisasemble(pdb_id, sdf_id)],
                local_scheduler=True)
    # check the number of full-matched pairs
    # grep -e '#tanimoto 1.0000' -e '#MoleculeA' /ddnB/work/jaydy/working/astex_disasemble.o399857 | grep '#MoleculeA'|sort|uniq|wc -l


if __name__ == '__main__':
    import sys
    pdb_id = sys.argv[1]
    main(pdb_id)
    # sdf_id = sys.argv[1]
    # check(sdf_id)
