#!/usr/bin/env python

import luigi
import random
import os
import subprocess32

from glob import glob
from astex_xcms import Path
from sample_confs import SampleConf, EvenlySampleConf, EvenlySampleConfGryration


class PairwiseCms(SampleConf):

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_cms.txt")
        return luigi.LocalTarget(path)

    def requires(self):
        return SampleConf(self.sdf_id, maximum_radius=self.maximum_radius)

    def trace_bin(self):
        return "/home/jaydy/Workspace/Bitbucket/xcms/src/astex_trace_sampled.bash"

    def getSdfs(self):
        regx = os.path.join(self.dirname(),
                            self.sdf_id, "*.sdf")
        return glob(regx)

    def removeSdfs(self):
        sdfs = self.getSdfs()
        for sdf in sdfs:
            os.remove(sdf)

    def run(self):
        self.removeSdfs()
        cmds = ['bash', self.trace_bin(),
                self.sdf_id,
                self.requires().output().path]
        _ = subprocess32.check_output(cmds)

        path_task = Path(self.sdf_id)
        pdb_path = path_task.pdb_path()

        sdfs = self.getSdfs()
        with open(self.output().path, 'w') as ofs:
            if len(sdfs) > 0:
                count = 0
                while count < len(sdfs) ** 2:
                    sdf1 = random.choice(sdfs)
                    sdf2 = random.choice(sdfs)
                    cmds = ['cms', '-frc',
                            '--lig1', sdf1,
                            '--lig2', sdf2,
                            '--prt1', pdb_path,
                            '--prt2', pdb_path]
                    if sdf1 != sdf2:
                        try:
                            stdout = subprocess32.check_output(cmds)
                            ofs.write(stdout)
                            count += 1
                        except:
                            pass
            else:
                pass


class PairwiseCmsEvenlyDistributed(EvenlySampleConf):

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_even_cms.txt")
        return luigi.LocalTarget(path)

    def requires(self):
        return EvenlySampleConf(self.sdf_id,
                                maximum_radius=self.maximum_radius)

    def trace_bin(self):
        return "/home/jaydy/Workspace/Bitbucket/xcms/src/astex_trace_sampled.bash"

    def getSdfs(self):
        regx = os.path.join(self.dirname(),
                            self.sdf_id, "*.sdf")
        return glob(regx)

    def removeSdfs(self):
        sdfs = self.getSdfs()
        for sdf in sdfs:
            os.remove(sdf)

    def run(self):
        self.removeSdfs()
        cmds = ['bash', self.trace_bin(),
                self.sdf_id,
                self.requires().output().path]
        _ = subprocess32.check_output(cmds)

        path_task = Path(self.sdf_id)
        pdb_path = path_task.pdb_path()

        sdfs = self.getSdfs()
        with open(self.output().path, 'w') as ofs:
            if len(sdfs) > 0:
                count = 0
                while count < len(sdfs):
                    sdf1 = random.choice(sdfs)
                    sdf2 = random.choice(sdfs)
                    cmds = ['cms', '-frc',
                            '--lig1', sdf1,
                            '--lig2', sdf2,
                            '--prt1', pdb_path,
                            '--prt2', pdb_path]
                    if sdf1 != sdf2:
                        try:
                            stdout = subprocess32.check_output(cmds)
                            ofs.write(stdout)
                            count += 1
                        except:
                            pass
            else:
                pass


class PairwiseCmsEvenlyGryration(PairwiseCmsEvenlyDistributed):

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_even_gyration_cms.txt")
        return luigi.LocalTarget(path)

    def requires(self):
        return EvenlySampleConfGryration(self.sdf_id)


def main(sdf):
    luigi.build([
        PairwiseCms(sdf, maximum_radius=8.0),
        PairwiseCmsEvenlyDistributed(sdf, maximum_radius=8.0),
        PairwiseCmsEvenlyGryration(sdf)
    ],
                local_scheduler=True)


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
