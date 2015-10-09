#!/usr/bin/env python

import os
import luigi
import glob
import subprocess32

from run_control_apoc import PkcombuAtomMatchParser


class Path(luigi.Task):

    sdf_id = luigi.Parameter()

    def astex_dir(self):
        return "/ddnB/work/jaydy/dat/astex"

    def complex_dir(self):
        return os.path.join(self.astex_dir(),
                            self.sdf_id)

    def pdb_path(self):
        return os.path.join(self.complex_dir(),
                            self.sdf_id[:5] + ".pdb")

    def sdf_path(self):
        return os.path.join(self.complex_dir(),
                            self.sdf_id + ".sdf")

    def working_dir(self):
        return os.path.join("/work/jaydy/working/astex_traces",
                            self.sdf_id)


def runPkcombu(a_path, b_path, oam_path):
    cmds = ["pkcombu",
            "-A", a_path,
            "-B", b_path,
            "-oam", oam_path]
    subprocess32.call(cmds)


class AstexXcmsBackCompatibility(Path):

    def output(self):
        path = os.path.join(self.working_dir(),
                            self.sdf_id + ".xcms.result")
        return luigi.LocalTarget(path)

    def requires(self):
        pass

    def run(self):

        def getConfNum(path):
            name = os.path.basename(path)
            return int(name.split('.')[0].split('_')[-1])

        def getXcmsVal(std_out):
            for line in std_out.splitlines():
                if "#xcms" in line:
                    return float(line.split()[-1])

        def getCmsVal(std_out):
            for line in std_out.splitlines():
                if "cms value:" in line:
                    return float(line.split()[-1])

        def getRmsdVal(std_out):
            for line in std_out.splitlines():
                if "rmsd value:" in line:
                    return float(line.split()[-1])

        def getFracVal(std_out):
            for line in std_out.splitlines():
                if "fraction value:" in line:
                    return float(line.split()[-1])

        variational_confs = glob.glob(os.path.join(self.working_dir(),
                                                   self.sdf_id + "_new*"))

        ofs = open(self.output().path, 'w')
        ofs.write("conf xcms cms rmsd fraction\n")
        for conf_path in sorted(variational_confs):
            conf_num = getConfNum(conf_path)
            oam_path = os.path.join(self.working_dir(),
                                    str(conf_num) + ".oam")
            runPkcombu(self.sdf_path(), conf_path, oam_path)
            lml_path = os.path.join(self.working_dir(),
                                    str(conf_num) + ".lml")
            PkcombuAtomMatchParser(oam_path).writeMatchingSerialNums(lml_path)

            xcms_cmds = ["run_xcms",
                         "--la", self.sdf_path(),
                         "--lb", conf_path,
                         "--pa", self.pdb_path(),
                         "--pb", self.pdb_path(),
                         "--lml", lml_path]
            std_out = subprocess32.check_output(xcms_cmds)
            xcms = getXcmsVal(std_out)

            cms_cmds = ["cms", "-frc",
                        "--lig1", self.sdf_path(),
                        "--lig2", conf_path,
                        "--prt1", self.pdb_path(),
                        "--prt2", self.pdb_path()]
            cms_out = subprocess32.check_output(cms_cmds)
            cms = getCmsVal(cms_out)
            frac = getFracVal(cms_out)
            rmsd = getRmsdVal(cms_out)

            ofs.write("%d %f %f %f %f\n" % (conf_num, xcms, cms, rmsd, frac))

        ofs.close()


def main(sdf_id):
    luigi.build([AstexXcmsBackCompatibility(sdf_id)],
                local_scheduler=True)

if __name__ == '__main__':
    import sys
    main(sys.argv[1])
