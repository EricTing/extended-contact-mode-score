#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import luigi
import random
import subprocess32


class Path(luigi.Task):

    sdf_id = luigi.Parameter()

    def pdb_id(self):
        return self.sdf_id[:5]

    def astex_dir(self):
        return os.path.join("/work/jaydy/dat/astex/",
                            self.sdf_id)

    def astex_sdf(self):
        return os.path.join(self.astex_dir(),
                            self.sdf_id + ".sdf")

    def astex_pdb(self):
        return os.path.join(self.astex_dir(),
                            self.sdf_id[:5] + ".pdb")

    def astex_pdb_conf_num(self):
        with open(self.astex_pdb(), 'r') as ifs:
            lines = ifs.readlines()
            return len([line for line in lines if 'MODEL' in line])

    def pdbqt(self, conf_num=0):
        return os.path.join("/work/jaydy/dat/vinaq",
                            self.sdf_id,
                            "%s.%d.pdbqt" % (self.pdb_id(), conf_num))

    def pdb(self, conf_num=0):
        return os.path.join("/work/jaydy/dat/splited_pdbs",
                            self.sdf_id,
                            self.pdb_id() + ".%d" % conf_num)

    def vina_bin(self):
        return "/project/michal/apps/autodock_vina_1_1_2_linux_x86/bin/vina"

    def ehalf_box_size_bin(self):
        return "/home/jaydy/Workspace/Bitbucket/ConfSpace/src/box-gyration.pl"

    def prepare_receptor4(self):
        return "/home/jaydy/Workspace/Bitbucket/ConfSpace/src/\
        prepare_receptor4.pl"

    def work_dir(self):
        path = os.path.join("/work/jaydy/working/vina",
                            self.sdf_id)
        try:
            os.makedirs(path)
        except:
            pass
        return path


class EvenlySampleConf(luigi.Task):

    sdf_id = luigi.Parameter()
    maximum_radius = luigi.Parameter(default=8.0)

    def dirname(self):
        return "/ddnB/work/jaydy/working/astex_weak"

    def requires(self):
        pass

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_even_grid.csv")
        return luigi.LocalTarget(path)

    def half_box_size(self):
        return self.maximum_radius / 2.0

    def run(self):
        ifn = os.path.join(self.dirname(),
                           self.sdf_id,
                           self.sdf_id + ".csv")

        dset = pd.read_csv(ifn, sep=' ', index_col=False)
        size = self.half_box_size()
        tra_pts = np.linspace(-size, size, 11)
        rot_pts = np.linspace(-3.14, 3.14, 7)

        data = []
        for x in tra_pts:
            for y in tra_pts:
                for z in tra_pts:
                    for r0 in rot_pts:
                        for r1 in rot_pts:
                            for r2 in rot_pts:
                                vector = [0, 0] + [x, y, z, r0, r1, r2]
                                data.append(vector)
        np.random.shuffle(data)
        total_samples = 100
        grid_cell_dset = pd.DataFrame(data[:total_samples],
                                      columns=dset.columns[:8])
        grid_cell_dset.to_csv(self.output().path, sep=' ', index=False)


class EvenlySampleConfGryration(EvenlySampleConf):

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_even_gyration_grid.csv")
        return luigi.LocalTarget(path)

    def half_box_size(self):
        path = Path(self.sdf_id)
        cmds = ['perl', path.ehalf_box_size_bin(), path.astex_sdf()]
        stdout = subprocess32.check_output(cmds)
        box_size, x, y, z = stdout.split()
        return float(box_size) / 2


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
            grid_cell_dset = pd.DataFrame(columns=dset.columns)

            x_grids = np.arange(dset.t0.min() - 0.1, dset.t0.max() + 0.1, dx)
            y_grids = np.arange(dset.t1.min() - 0.1, dset.t1.max() + 0.1, dy)
            z_grids = np.arange(dset.t2.min() - 0.1, dset.t2.max() + 0.1, dz)

            for x in x_grids:
                for y in y_grids:
                    for z in z_grids:
                        mydset = dset[(dset.t0 > x) &
                                      (dset.t0 < x + dx) &
                                      (dset.t1 > y) &
                                      (dset.t1 < y + dy) &
                                      (dset.t2 > z) &
                                      (dset.t2 < z + dz)]
                        if len(mydset) > 0:
                            mydset.index = range(len(mydset))
                            rnd_idx = random.choice(mydset.index)
                            grid_cell_dset = pd.concat([grid_cell_dset,
                                                        mydset.ix[[rnd_idx]]])

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
