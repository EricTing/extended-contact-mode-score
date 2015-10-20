#!/usr/bin/env python

import luigi
import random
import os
import subprocess32

from glob import glob

from astex_xcms import Path
from sample_confs import SampleConf


class PairwiseCms(SampleConf):

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_cms.txt")
        return luigi.LocalTarget(path)

    def requires(self):
        return SampleConf(self.sdf_id)

    def trace_bin(self):
        return "/home/jaydy/Workspace/Bitbucket/xcms/src/astex_trace_sampled.bash"

    def run(self):
        cmds = ['bash', self.trace_bin(),
                self.sdf_id]
        _ = subprocess32.check_output(cmds)

        regx = os.path.join(self.dirname(),
                            self.sdf_id, "*.sdf")

        path_task = Path(self.sdf_id)
        pdb_path = path_task.pdb_path()

        sdfs = glob(regx)
        with open(self.output().path, 'w') as ofs:
            count = 0
            while count < 200:
                sdf1 = random.choice(sdfs)
                sdf2 = random.choice(sdfs)
                cmds = ['cms', '-frc',
                        '--lig1', sdf1,
                        '--lig2', sdf2,
                        '--prt1', pdb_path,
                        '--prt2', pdb_path]
                try:
                    stdout = subprocess32.check_output(cmds)
                    ofs.write(stdout)
                    count += 1
                except:
                    pass


def main(sdf):
    luigi.build([PairwiseCms(sdf)],
                local_scheduler=True)


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
