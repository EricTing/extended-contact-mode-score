#!/usr/bin/env python

import luigi
import math
import pybel
import subprocess32
import os
import pandas as pd
import numpy as np
import json
import random
from scipy.spatial.distance import euclidean
from collections import defaultdict
from sample_confs import Path
from scipy.stats import spearmanr, pearsonr
from glob import glob


class MyPath(Path):
    def work_dir(self):
        path = os.path.join("/work/jaydy/working/single_cms",
                            self.sdf_id)
        try:
            os.makedirs(path)
        except:
            pass
        return path

    def trace_bin(self):
        return "/home/jaydy/Workspace/Bitbucket/xcms/src/astex_trace_sampled.bash"


class VariousBoxSizeCms(luigi.Task):
    sdf_id = luigi.Parameter()
    box_size = luigi.Parameter()
    max_rot_angle = luigi.Parameter(default=3.14)
    num_samples = luigi.Parameter(default=100)
    num_cms_vals = luigi.Parameter(default=100)

    def work_dir(self):
        my_path = MyPath(self.sdf_id)
        path = os.path.join(my_path.work_dir(),
                            str(self.box_size))
        try:
            os.makedirs(path)
        except:
            pass
        return path

    def requires(self):
        pass

    def output(self):
        path = os.path.join(self.work_dir(),
                            'rnd_grid.csv')
        return luigi.LocalTarget(path)

    def half_box_size(self):
        return self.box_size / 2.0

    def run(self):
        data = []
        for _ in range(self.num_samples):
            tra_vec = np.random.uniform(-self.half_box_size(),
                                        self.half_box_size(),
                                        3)
            rot_vec = np.random.uniform(-self.max_rot_angle,
                                        self.max_rot_angle,
                                        3)
            conf_nums = [0, 0]
            data.append(np.hstack((conf_nums, tra_vec, rot_vec)))

        cols = ['lig', 'prt', 't0', 't1', 't2', 'r0', 'r1', 'r2']
        dset = pd.DataFrame(data, columns=cols)
        dset.to_csv(self.output().path, sep=' ', index=False)


class PairwiseDisSimilarity(VariousBoxSizeCms):

    def requires(self):
        return VariousBoxSizeCms(self.sdf_id,
                                 self.box_size,
                                 num_cms_vals=self.num_cms_vals,
                                 max_rot_angle=self.max_rot_angle,
                                 num_samples=self.num_samples)

    def output(self):
        path = os.path.join(self.work_dir(),
                            'rnd_grid.json')
        return luigi.LocalTarget(path)

    def getSdfs(self):
        regx = os.path.join(self.work_dir(), '*.sdf')
        return glob(regx)

    def removeSdfs(self):
        for sdf in self.getSdfs():
            os.remove(sdf)

    def readOutput(self):
        return json.loads(open(self.output().path, 'r').read())

    def run(self):
        self.removeSdfs()
        my_path = MyPath(self.sdf_id)
        cmds = ['bash', my_path.trace_bin(),
                self.sdf_id,
                self.requires().output().path,
                self.work_dir()]
        subprocess32.call(cmds)

        sdfs = self.getSdfs()
        data = defaultdict(list)
        pdb = my_path.astex_pdb()
        native_sdf = my_path.astex_sdf()
        for time in range(self.num_cms_vals):
            var_sdf = random.choice(sdfs)
            cmds = ['cms', '-frc',
                    '--lig1', native_sdf,
                    '--lig2', var_sdf,
                    '--prt1', pdb,
                    '--prt2', pdb]
            stdout = subprocess32.check_output(cmds)
            for line in stdout.splitlines():
                if 'rmsd' in line:
                    rmsd = float(line.split()[-1])
                    data['rmsd'].append(rmsd)
                if 'cms' in line:
                    cms = float(line.split()[-1])
                    data['cms'].append(cms)
                if 'fraction' in line:
                    fraction = float(line.split()[-1])
                    data['fraction'].append(fraction)
        to_write = json.dumps(data, sort_keys=True,
                              indent=4, separators=(',', ': '))
        with open(self.output().path, 'w') as ofs:
            ofs.write(to_write)


def squared_distance(coordsA, coordsB):
    """Find the squared distance between two 3-tuples"""
    sqrdist = sum((a-b)**2 for a, b in zip(coordsA, coordsB))
    return sqrdist


def rmsd(allcoordsA, allcoordsB):
    """Find the RMSD between two lists of 3-tuples"""
    deviation = sum(squared_distance(atomA, atomB) for
                    (atomA, atomB) in zip(allcoordsA, allcoordsB))
    return math.sqrt(deviation / float(len(allcoordsA)))


def contactVector(lig_coords, prt_coords):
    vec = [math.sqrt(squared_distance(c1, c2)) for c1 in lig_coords
           for c2 in prt_coords]
    return vec


def readLigCoords(path, lig_format="sdf"):
    mol = pybel.readfile(lig_format, path).next()
    return [atom.coords for atom in mol if not atom.OBAtom.IsHydrogen()]


def readPrtCoords(path, prt_format="pdb"):
    mol = pybel.readfile(prt_format, path).next()
    return [atom.coords for atom in mol if not atom.OBAtom.IsHydrogen()]


def compare(lig1_path, lig2_path, prt1_path, prt2_path,
            lig_format="sdf", prt_format="pdb"):
    lig1 = readLigCoords(lig1_path)
    lig2 = readLigCoords(lig2_path)

    rmsd_val = rmsd(lig1, lig2)

    prt1 = readPrtCoords(prt1_path)
    prt2 = readPrtCoords(prt2_path)

    vec1 = contactVector(lig1, prt1)
    vec2 = contactVector(lig2, prt2)

    masked_vec1 = [0 if val > 4.8 else 1 for val in vec1]
    masked_vec2 = [0 if val > 4.8 else 1 for val in vec2]

    contact_rmsd = euclidean(vec1, vec2) / math.sqrt(len(vec1))

    neighbor_vec1, neighbor_vec2 = [], []
    for d1, d2 in zip(vec1, vec2):
        if d1 < 7.5 and d2 < 7.5:
            neighbor_vec1.append(d1)
            neighbor_vec2.append(d2)

    return {
        "rmsd": rmsd_val,
        "spearmanr": spearmanr(vec1, vec2),
        "pearsonr": pearsonr(vec1, vec2),
        "masked_pearsonr": pearsonr(masked_vec1, masked_vec2),
        "masked_spearmanr": spearmanr(masked_vec1, masked_vec2),
        "contact_rmsd": contact_rmsd,
        "neighbor_pearsonr": pearsonr(neighbor_vec1, neighbor_vec2),
        "neighbor_spearmanr": spearmanr(neighbor_vec1, neighbor_vec2),
    }


class PairwiseDisSimilaritySpearman(PairwiseDisSimilarity):
    def output(self):
        path = os.path.join(self.work_dir(),
                            'spearmanr.json')
        return luigi.LocalTarget(path)

    def run(self):
        self.removeSdfs()
        my_path = MyPath(self.sdf_id)
        cmds = ['bash', my_path.trace_bin(),
                self.sdf_id,
                self.requires().output().path,
                self.work_dir()]
        subprocess32.call(cmds)

        sdfs = self.getSdfs()
        data = defaultdict(list)
        pdb = my_path.astex_pdb()
        native_sdf = my_path.astex_sdf()

        ligs = [readLigCoords(sdf) for sdf in sdfs]
        prt = readPrtCoords(pdb)
        native_lig = readLigCoords(native_sdf)
        native_vec = contactVector(native_lig, prt)
        for lig in ligs:
            vec = contactVector(lig, prt)
            max_pmf = 7.38 + 0.1
            neighboring_vec, neighboring_native_vec = [], []
            for dist, native_dist in zip(vec, native_vec):
                if dist < max_pmf or native_dist < max_pmf:
                    neighboring_vec.append(dist)
                    neighboring_native_vec.append(native_dist)

            contact_pairs = len(neighboring_native_vec)
            if contact_pairs > 3:
                rmsd_val = rmsd(lig, native_lig)
                spearman_corr = spearmanr(vec, native_vec)
                pearson_corr = pearsonr(vec, native_vec)

                n_spearman_corr = spearmanr(neighboring_vec,
                                            neighboring_native_vec)
                n_pearson_corr = pearsonr(neighboring_vec, neighboring_native_vec)

                data['rmsd_val'].append(rmsd_val)
                data['pearsonr'].append(pearson_corr[0])
                data['pearsonr_p_val'].append(pearson_corr[1])
                data['spearmanr'].append(spearman_corr.correlation)
                data['spearmanr_p_val'].append(spearman_corr.pvalue)
                data['neighboring_pearsonr'].append(n_pearson_corr[0])
                data['neighboring_pearsonr_p_val'].append(n_pearson_corr[1])
                data['neighboring_spearmanr'].append(n_spearman_corr.correlation)
                data['neighboring_spearmanr_p_val'].append(n_spearman_corr.pvalue)
                data['contact_pairs'].append(len(neighboring_vec))

        to_write = json.dumps(data, sort_keys=True,
                              indent=4, separators=(',', ': '))
        with open(self.output().path, 'w') as ofs:
            ofs.write(to_write)


def main(sdf):
    luigi.build([
        # PairwiseDisSimilarity(sdf, 5.0),
        # PairwiseDisSimilarity(sdf, 10.0),
        PairwiseDisSimilarity(sdf, 10.68, num_cms_vals=1000, num_samples=500),
        PairwiseDisSimilaritySpearman(sdf, 10.68, num_samples=100),
        # PairwiseDisSimilarity(sdf, 15.0),
        # PairwiseDisSimilarity(sdf, 20.0)
    ],
                local_scheduler=True)

    pass


if __name__ == '__main__':
    import sys
    main(sys.argv[1])

