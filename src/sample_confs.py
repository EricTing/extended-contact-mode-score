#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import luigi
import random


class SampleConf(luigi.Task):

    sdf_id = luigi.Parameter()
    radius = luigi.Parameter(default=8.0)

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

        def calculateDist2Center():
                vals = dset[['t0', 't1', 't2']].values
                return np.sqrt(np.sum(np.multiply(vals, vals), axis=1))

        dset['dst'] = calculateDist2Center()
        dset = dset[dset['dst'] < self.radius]

        dx, dy, dz = [0.1] * 3
        grid_pts_with_confs = pd.DataFrame()
        for x in np.arange(dset.t0.min(), dset.t0.max(), dx):
            for y in np.arange(dset.t1.min(), dset.t1.max(), dx):
                for z in np.arange(dset.t2.min(), dset.t2.max(), dx):
                    mydset = dset[(dset.t0 > x)
                                  & (dset.t0 < x + dx)
                                  & (dset.t1 > y)
                                  & (dset.t1 < y + dy)
                                  & (dset.t2 > z)
                                  & (dset.t2 < z + dz)]
                    if len(mydset) > 0:
                        mydset.index = range(len(mydset))
                        grid_pts_with_confs = pd.concat([grid_pts_with_confs,
                                                         mydset.ix[[random.choice(mydset.index)]]])

        num_clusters = 500
        if len(grid_pts_with_confs) > num_clusters:
            rows = np.random.choice(grid_pts_with_confs.index.values,
                                    num_clusters)
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
