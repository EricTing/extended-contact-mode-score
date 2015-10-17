#!/usr/bin/env python

import luigi
import pandas as pd
from cms_vals_for_p_val import PairwiseCms


class CollectGridCms(luigi.Task):

    def requires(self):
        sdfs = [_.rstrip() for _ in file("../dat/astex_sdf.lst")]
        return [PairwiseCms(sdf) for sdf in sdfs]

    def output(self):
        path = "../dat/astex_grid_cms.csv"
        return luigi.LocalTarget(path)

    def run(self):
        dset = {}
        for job in self.requires():
            with job.output().open('r') as ifs:
                lines = ifs.readlines()
                cms_lines = filter(lambda line: "cms value" in line, lines)
                cms_valus = map(lambda line: float(line.split()[-1]),
                                cms_lines)
                sdf_id = job.sdf_id
                dset[sdf_id] = cms_valus

        dset = pd.DataFrame(dset)
        dset.to_csv(self.output().path)


def main():
    luigi.build([CollectGridCms()],
                local_scheduler=True)
    pass


if __name__ == '__main__':
    main()
