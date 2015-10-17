#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import luigi


class SampleConf(luigi.Task):

    sdf_id = luigi.Parameter()

    def dirname(self):
        return "/ddnB/work/jaydy/working/astex_weak"

    def requires(self):
        pass

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_grid.csv")
        return luigi.LocalTarget(path)

    def run(self):
        ifn = os.path.join(self.dirname(),
                           self.sdf_id,
                           self.sdf_id + ".csv")

        dset = pd.read_csv(ifn, sep=' ', index_col=False)
        total_samples = 20000
        if len(dset > total_samples):
            rows = np.random.choice(dset.index.values, total_samples)
            dset = dset.ix[rows]

        xg = np.linspace(dset.t0.min(), dset.t0.max(), 10)
        yg = np.linspace(dset.t1.min(), dset.t1.max(), 10)
        zg = np.linspace(dset.t2.min(), dset.t2.max(), 10)

        dx = xg[1] - xg[0]
        dy = yg[1] - yg[0]
        dz = zg[1] - zg[0]

        grid_pts_with_confs = pd.DataFrame()
        for x in xg:
            for y in yg:
                for z in zg:
                    mydset = dset[(dset.t0 > x)
                                  & (dset.t0 < x + dx)
                                  & (dset.t1 > y)
                                  & (dset.t1 < y + dy)
                                  & (dset.t2 > z)
                                  & (dset.t2 < z + dz)]
                    if len(mydset) > 0:
                        mydset.index = range(len(mydset))
                        grid_pts_with_confs = pd.concat([grid_pts_with_confs, mydset.ix[[0]]])

        num_clusters = 100
        if len(grid_pts_with_confs) > num_clusters:
            rows = np.random.choice(grid_pts_with_confs.index.values, num_clusters)
            grid_pts_with_confs = grid_pts_with_confs.ix[rows]

        cols = dset.columns
        grid_pts_with_confs = grid_pts_with_confs[cols]
        ofn = self.output().path
        grid_pts_with_confs.to_csv(ofn, sep=' ', index=False)


def main(sdf):
    luigi.build([SampleConf(sdf)],
                local_scheduler=True)


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
