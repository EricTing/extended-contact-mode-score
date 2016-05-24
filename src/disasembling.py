#!/usr/bin/env python

import luigi
import shlex
import os
import subprocess32

from run_control_apoc import PdbMol2, PkcombuAtomMatchParser
from myurls import WORKING_DIR
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
        return luigi.LocalTarget(self.dirname())

    def run(self):
        cif_ifn = self.requires().output().path
        cmd = "mmcif_disasembler --ifn %s --dir %s" % (cif_ifn, self.dirname())
        print cmd
        args = shlex.split(cmd)
        subprocess32.call(args)


class Code:
    def __init__(self, code):
        self.code = code
        tokens = self.code.split('_')
        self.pdb_id = tokens[0]
        self.lig_id = tokens[1]
        self.chain_id = tokens[2]


class ApocRandomSet2(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()

    def my_midtwo(self):
        return self.tname[1:3]

    def my_dir(self):
        path = os.path.join(WORKING_DIR,
                            self.my_midtwo(),
                            self.tname)
        try:
            os.makedirs(path)
        except:
            pass
        return path

    def requires(self):
        tcode = Code(self.tname)
        qcode = Code(self.qname)
        return [DisasemblePdb(tcode.pdb_id),
                DisasemblePdb(qcode.pdb_id)]

    def output(self):
        path = os.path.join(self.my_dir(),
                            self.tname + '__' + self.qname + ".xcms")
        return luigi.LocalTarget(path)

    def runKcombu(self, lig1_path, lig2_path):
        oam_path = os.path.join(self.my_dir(),
                                self.tname + '__' + self.qname + ".oam")
        cmds = ["pkcombu",
                "-A", lig1_path,
                "-B", lig2_path,
                "-oam", oam_path]
        try:
            print " ".join(cmds)
            subprocess32.call(cmds)
            return oam_path
        except:
            pass

    def PrtAndLig(self, name):
        code = Code(name)
        pdb_id = code.pdb_id
        src_dir = DisasemblePdb(pdb_id).output().path
        mol_fn = "%s_%s%s_%s.mol2" % (pdb_id.upper(), pdb_id.upper(),
                                      code.chain_id,
                                      code.lig_id.upper())
        pdb_fn = "%s_%s%s.pdb" % (pdb_id.upper(),
                                  pdb_id.upper(), code.chain_id)
        return os.path.join(src_dir, pdb_fn), os.path.join(src_dir, mol_fn)

    def run(self):
        pdb1, mol1 = self.PrtAndLig(self.tname)
        pdb2, mol2 = self.PrtAndLig(self.qname)
        paths = [pdb1, pdb2, mol1, mol2]
        if not all(map(os.path.exists, paths)):
            print "not enough inputs"
        else:
            oam_path = self.runKcombu(mol1, mol2)
            parser = PkcombuAtomMatchParser(oam_path)
            output_path = self.output().path
            lml_path = os.path.splitext(output_path)[0] + '.lml'
            parser.writeMatchingSerialNums(lml_path)
            cmds = ['run_xcms',
                    '--la', mol1,
                    '--lb', mol2,
                    '--pa', pdb1,
                    '--pb', pdb2,
                    '--lml', lml_path]
            print " ".join(cmds)
            stdout = subprocess32.check_output(cmds)
            with open(output_path, 'w') as ofs:
                ofs.write(stdout)


def main(tname, qname):
    luigi.build([ApocRandomSet2(tname, qname)],
                local_scheduler=True)


def check(sdf_id):
    pdb_id = sdf_id[:4]
    luigi.build([CheckAstexDisasemble(pdb_id, sdf_id)],
                local_scheduler=True)
    # check the number of full-matched pairs
    # grep -e '#tanimoto 1.0000' -e '#MoleculeA' /ddnB/work/jaydy/working/astex_disasemble.o399857 | grep '#MoleculeA'|sort|uniq|wc -l


if __name__ == '__main__':
    import sys
    tname = sys.argv[1]
    qname = sys.argv[2]
    main(tname, qname)
    # sdf_id = sys.argv[1]
    # check(sdf_id)
