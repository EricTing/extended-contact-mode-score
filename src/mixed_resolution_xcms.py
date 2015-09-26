#!/usr/bin/env python


import luigi
import os
import subprocess32
from apoc_inputs import DecompressedPdb
from run_control_apoc import LigandExpStructureInMol2
from run_control_apoc import LigandMatchingList
from run_control_apoc import ProteinMatchingList
from urls import APOC_WORKING_DIR


class MixedResolutionXcms(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()
    subset = luigi.Parameter(default="subject")

    def _convert(self, name):
        return name.split('_')

    def my_midtwo(self):
        return self.tname[1:3]

    def mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.subset,
                                self.my_midtwo())
        return dir_path

    def requires(self):
        t_pdb = self.tname.split('_')[0]
        q_pdb = self.qname.split('_')[0]
        return [LigandExpStructureInMol2(self.tname),
                LigandExpStructureInMol2(self.qname),
                LigandMatchingList(self.tname, self.qname, self.subset),
                ProteinMatchingList(self.tname, self.qname, self.subset),
                DecompressedPdb(t_pdb),
                DecompressedPdb(q_pdb)]

    def output(self):
        path = os.path.join(self.mydir(),
                            "%s_%s.mxcms" % (self.tname, self.qname))
        return luigi.LocalTarget(path)

    def _run(self):
        tname = self.tname
        qname = self.qname

        t_pdb = tname.split('_')[0]
        q_pdb = qname.split('_')[0]

        t_lig_path = LigandExpStructureInMol2(tname).output().path
        q_lig_path = LigandExpStructureInMol2(qname).output().path

        t_prt_path = DecompressedPdb(t_pdb).output().path
        q_prt_path = DecompressedPdb(q_pdb).output().path

        lml_path = LigandMatchingList(tname, qname, self.subset).output().path
        pml_path = ProteinMatchingList(tname, qname, self.subset).output().path

        cmds = ["run_xcms",
                "--la", t_lig_path,
                "--lb", q_lig_path,
                "--pa", t_prt_path,
                "--pb", q_prt_path,
                "--lml", lml_path,
                "--pml", pml_path]

        print " ".join(cmds)
        subprocess32.call(cmds)

    def run(self):
        self._run()


def main(tname, qname):
    if tname != "tname":
        luigi.build([MixedResolutionXcms(tname, qname, "subject")],
                    local_scheduler=True)


if __name__ == '__main__':
    import sys
    tname, qname = sys.argv[1], sys.argv[2]
    main(tname, qname)
