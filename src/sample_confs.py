#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import luigi
import random


class SampleConf(luigi.Task):

    sdf_id = luigi.Parameter()
    maximum_radius = luigi.Parameter(default=10.0)
    minimum_radius = luigi.Parameter(default=00.0)

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

        def resort(angle):
            while abs(angle) > 3.14:
                if angle > 0:
                    angle -= 6.28
                else:
                    angle += 6.28
            return angle

        dset[['r0', 'r1', 'r2']] = dset[['r0', 'r1', 'r2']].applymap(resort)

        if not dset['dst'].max() > self.minimum_radius:
            ofs = open(self.output().path, 'w')
            ofs.close()
        else:
            dset = dset[dset['dst'] < self.maximum_radius]

            dx, dy, dz = [2.0] * 3
            d_angle = 1.256        # 72 degree
            grid_cell_dset = pd.DataFrame(columns=dset.columns)

            x_grids = np.arange(dset.t0.min() - 0.1, dset.t0.max() + 0.1, dx)
            y_grids = np.arange(dset.t1.min() - 0.1, dset.t1.max() + 0.1, dy)
            z_grids = np.arange(dset.t2.min() - 0.1, dset.t2.max() + 0.1, dz)
            r0_grids = np.arange(dset.r0.min() - 0.1, dset.r0.max() + 0.1, d_angle)
            r1_grids = np.arange(dset.r1.min() - 0.1, dset.r1.max() + 0.1, d_angle)
            r2_grids = np.arange(dset.r2.min() - 0.1, dset.r2.max() + 0.1, d_angle)

            for x in x_grids:
                for y in y_grids:
                    for z in z_grids:
                        for r0 in r0_grids:
                            for r1 in r1_grids:
                                for r2 in r2_grids:
                                    mydset = dset[(dset.t0 > x)
                                                  & (dset.t0 < x + dx)
                                                  & (dset.t1 > y)
                                                  & (dset.t1 < y + dy)
                                                  & (dset.t2 > z)
                                                  & (dset.t2 < z + dz)
                                                  & (dset.r0 > r0)
                                                  & (dset.r0 < r0 + d_angle)
                                                  & (dset.r1 > r1)
                                                  & (dset.r1 < r1 + d_angle)
                                                  & (dset.r2 > r2)
                                                  & (dset.r2 < r2 + d_angle)]
                                    if len(mydset) > 0:
                                        mydset.index = range(len(mydset))
                                        grid_cell_dset = pd.concat([grid_cell_dset,
                                                                    mydset.ix[[random.choice(mydset.index)]]])

            grid_cell_dset.index = range(len(grid_cell_dset))
            ofn = self.output().path
            grid_cell_dset.to_csv(ofn, sep=' ', index=False)


def main(sdf):
    luigi.build([SampleConf(sdf)],
                local_scheduler=True)


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
