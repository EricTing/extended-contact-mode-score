#!/usr/bin/env python


import luigi
import os
import subprocess32
import json
from apoc_inputs import DecompressedPdb
from run_control_apoc import LigandExpStructureInMol2
from run_control_apoc import LigandMatchingList
from run_control_apoc import ProteinMatchingList
from run_control_apoc import LpcKcombuResult, LpcApocResultTask
from run_control_apoc import PkcombuAtomMatchParser, ApocResultParer
from myurls import APOC_WORKING_DIR


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
                LpcKcombuResult(self.tname, self.qname, self.subset),
                LpcApocResultTask(self.tname, self.qname, self.subset),
                DecompressedPdb(t_pdb),
                DecompressedPdb(q_pdb)]

    def output(self):
        path = os.path.join(self.mydir(),
                            "%s__%s_mixedres.json" % (self.tname, self.qname))
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
        stdout = subprocess32.check_output(cmds)
        data = {}
        for line in stdout.splitlines():
            if "#atom" in line:
                data['# ligand atoms'] = int(line.split()[-1])
            if "#res" in line:
                data['# residues'] = int(line.split()[-1])
            if "#xcms" in line:
                data['xcms'] = float(line.split()[-1])
        return data

    def _kcombu_results(self):
        path = LpcKcombuResult(self.tname,
                               self.qname,
                               self.subset).output().path
        return PkcombuAtomMatchParser(path)

    def run(self):
        data = self._run()
        data['Kcombu tanimoto'] = self._kcombu_results().data.tanimoto

        with LpcApocResultTask(self.tname,
                               self.qname,
                               self.subset).output().open('r') as f:
            apoc_parser = ApocResultParer(f.read())

        data['tname'] = self.tname
        data['qname'] = self.qname
        data['Apoc result'] = LpcApocResultTask(self.tname, self.qname, self.subset).output().path
        data['Kcombu result'] = LpcKcombuResult(self.tname, self.qname, self.subset).output().path

        global_alignment = apoc_parser.queryGlobal(self.tname, self.qname)
        data['seq identity'] = global_alignment.seq_identity
        pocket_alignment = apoc_parser.queryPocket(self.tname, self.qname)
        if pocket_alignment.has_pocket_alignment:
            data['Apoc ps-score'] = pocket_alignment.ps_score
            data['Apoc p-value'] = pocket_alignment.p_value
            data['Apoc z-score'] = pocket_alignment.z_score

        f = self.output().open('w')
        to_write = json.dumps(data, sort_keys=True, indent=4, separators=(',', ': '))
        f.write(to_write)
        f.close()

        print "xcms output %s" % (self.output().path)


def main(tname, qname):
    if tname != "tname":
        luigi.build([MixedResolutionXcms(tname, qname, "subject")],
                    local_scheduler=True)


if __name__ == '__main__':
    import sys
    tname, qname = sys.argv[1], sys.argv[2]
    main(tname, qname)
