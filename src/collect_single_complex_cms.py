#!/usr/bin/env python

import luigi
import subprocess32
import os
import json
import pandas as pd
from single_complex_cms import PairwiseDisSimilarity, MyPath
from single_complex_cms import PairwiseDisSimilaritySpearman


class GyrationSize(luigi.Task):
    def getSdfIds(self):
        return [_.rstrip() for _ in file("../dat/astex_sdf.lst")]

    def requires(self):
        pass

    def output(self):
        path = os.path.join("../dat", "gyration.txt")
        return luigi.LocalTarget(path)

    def run(self):
        data = []
        for sdf_id in self.getSdfIds():
            my_path = MyPath(sdf_id)
            cmds = ['perl', my_path.ehalf_box_size_bin(), my_path.astex_sdf()]
            stdout = subprocess32.check_output(cmds)
            data.append([sdf_id] + stdout.split())
        dset = pd.DataFrame(data, columns=['sdf', 'box_size', 'x', 'y', 'z'])
        dset.to_csv(self.output().path)


class CollectPairwiseDisSimilarity(luigi.Task):

    box_size = luigi.Parameter()

    def getSdfIds(self):
        return [_.rstrip() for _ in file("../dat/astex_sdf.lst")]

    def requires(self):
        return [PairwiseDisSimilarity(sdf_id, self.box_size)
                for sdf_id in self.getSdfIds()]

    def output(self):
        path = os.path.join("../dat", "rnd_%f.json" % self.box_size)
        return luigi.LocalTarget(path)

    def run(self):
        data = {}
        for pair_wise_job in self.requires():
            data[pair_wise_job.sdf_id] = pair_wise_job.readOutput()

        to_write = json.dumps(data, sort_keys=True,
                              indent=4, separators=(',', ': '))
        with open(self.output().path, 'w') as ofs:
            ofs.write(to_write)


class CollectPairwiseDisSimilaritySpearman(CollectPairwiseDisSimilarity):
    def requires(self):
        return [PairwiseDisSimilaritySpearman(sdf, 10.68, num_samples=100)
                for sdf in self.getSdfIds()]

    def output(self):
        path = os.path.join("../dat/spearmanr.csv")
        return luigi.LocalTarget(path)

    def run(self):
        all_dset = pd.DataFrame()
        for task in self.requires():
            dset = pd.read_json(task.output().path)
            dset['sdf'] = [task.sdf_id] * len(dset)
            all_dset = all_dset.append(dset)

        all_dset.to_csv(self.output().path)


def main():
    luigi.build([
        CollectPairwiseDisSimilarity(5.0),
        CollectPairwiseDisSimilarity(10.0),
        CollectPairwiseDisSimilarity(10.68),
        CollectPairwiseDisSimilarity(15.0),
        CollectPairwiseDisSimilarity(20.0),
        GyrationSize(),
        CollectPairwiseDisSimilaritySpearman(10.68),
    ],
                local_scheduler=True)


if __name__ == '__main__':
    main()
