#!/usr/bin/env python

import luigi
import itertools
import os
import pandas as pd
import run_control_apoc
from urls import WORKING_DIR
from collect_xcms import AtomicXcmsTable
from sklearn.cluster import DBSCAN
import numpy as np


class SuccessfulPocketList(luigi.Task):

    subset = luigi.Parameter()

    def requires(self):
        return AtomicXcmsTable(self.subset)

    def output(self):
        path = os.path.join(WORKING_DIR, self.subset + "_xcms_success.lst")
        return luigi.LocalTarget(path)

    def readLst(self):
        with self.output().open('r') as inputObj:
            return inputObj.read().splitlines()

    def run(self):
        ifn = self.requires().output().path
        dset = pd.read_csv(ifn, index_col=0)
        good_pockets = set(dset.tname)

        with self.output().open('w') as outputObj:
            for tname in good_pockets:
                outputObj.write("%s\n" % (tname))


class Curate:

    class LpcApocResultTask(run_control_apoc.LpcApocResultTask):

        def _my_midtwo(self):
            return self.qname[1:3]

        def _mydir(self):
            return os.path.join(WORKING_DIR,
                                self.subset,
                                self.tname[1:3],
                                self.tname,
                                self._my_midtwo())

    class PairWisePsScore(luigi.Task):

        tname = luigi.Parameter()
        subset = luigi.Parameter()

        def requires(self):
            list_path = SuccessfulPocketList(self.subset).output().path
            assert(os.path.exists(list_path))
            with open(list_path, 'r') as inputObj:
                names = inputObj.read().splitlines()

            my_idx = names.index(self.tname)
            qnames = []
            for idx, name in enumerate(names):
                if idx > my_idx:
                    qnames.append(name)

            return [Curate.LpcApocResultTask(self.tname, qname, self.subset)
                    for qname in qnames]

        def output(self):
            path = os.path.join(WORKING_DIR,
                                self.subset,
                                self.tname[1:3],
                                self.tname + "_ps_scores.txt")
            return luigi.LocalTarget(path)

        def _run(self):
            with self.output().open('w') as outputObj:
                for apoc_task in self.requires():
                    with apoc_task.output().open('r') as inputObj:
                        parser = run_control_apoc.ApocResultParer(inputObj.read())
                        qname = apoc_task.qname
                        pocket = parser.queryPocket(self.tname, qname)
                        ps_score = pocket.ps_score
                        outputObj.write("%s %s %f\n" % (self.tname,
                                                        qname,
                                                        ps_score))

        def run(self):
            self._run()

    class PocketSimilarityMatrix(luigi.Task):

        subset = luigi.Parameter()

        def requires(self):
            list_path = SuccessfulPocketList(self.subset).output().path
            assert(os.path.exists(list_path))
            with open(list_path, 'r') as inputObj:
                names = inputObj.read().splitlines()

            return [Curate.PairWisePsScore(tname, self.subset)
                    for tname in names]

        def output(self):
            path = os.path.join(WORKING_DIR, self.subset,
                                self.subset + "_ps_score_matrix.txt")
            return luigi.LocalTarget(path)

        def _readPsScores(self, path):
            with open(path, 'r') as inputObj:
                return [t.split()[-1] for t in inputObj.read().splitlines()]

        def run(self):
            ps_score_results = self.requires()
            size = len(ps_score_results)
            matrix = [[1.0 for x in range(size)]
                      for y in range(size)]

            for idx, task in enumerate(ps_score_results):
                ifn = task.output().path
                ps_scores = self._readPsScores(ifn)
                matrix[idx][(idx + 1):] = ps_scores[:(size - idx - 1)]

            # transpose since it is a similarity matrix
            for idx in range(size):
                for idy in range(size):
                    if idy < idx:
                        matrix[idx][idy] = matrix[idy][idx]

            with self.output().open('w') as outputObj:
                for row in matrix:
                    outputObj.write(" ".join(map(str, row)))
                    outputObj.write("\n")

    class ReSample(luigi.Task):

        subset = luigi.Parameter()
        eps = luigi.Parameter(default=1.15)

        def requires(self):
            return [Curate.PocketSimilarityMatrix(self.subset),
                    AtomicXcmsTable(self.subset),
                    SuccessfulPocketList(self.subset)]

        def output(self):
            matrix_path = AtomicXcmsTable(self.subset).output().path
            path = os.path.join(os.path.dirname(matrix_path),
                                self.subset + "_curated_xcms.csv")
            return luigi.LocalTarget(path)

        def run(self):
            ifn = Curate.PocketSimilarityMatrix(self.subset).output().path
            mat = np.loadtxt(ifn)
            db = DBSCAN(eps=self.eps, min_samples=1).fit(mat)
            labels = db.labels_
            print len(set(labels))
            names = SuccessfulPocketList(self.subset).readLst()
            assert(len(names) == len(labels))

            groups = itertools.groupby(zip(labels, names), key=lambda t: t[0])
            sampled_names = []
            for _, group in groups:
                tnames = [t[-1] for t in list(group)]
                sampled = np.random.choice(tnames)
                sampled_names.append(sampled)

            dset = AtomicXcmsTable(self.subset).readDset()
            sampled = dset[dset['tname'].isin(sampled_names)]
            sampled.to_csv(self.output().path)


def test():
    luigi.build([SuccessfulPocketList(subset="subject"),
                 Curate.LpcApocResultTask(subset="subject",
                                          tname="3oz2_OZ2_A_502",
                                          qname="3ocj_PLM_A_305"),
                 Curate.PocketSimilarityMatrix(subset="subject"),
                 Curate.PairWisePsScore(tname="3oz2_OZ2_A_502",
                                        subset="subject")],
                local_scheduler=True)


def main():
    luigi.build([SuccessfulPocketList(subset="subject"),
                 Curate.ReSample(subset="subject"),
                 Curate.PocketSimilarityMatrix(subset="subject")],
                local_scheduler=True)


if __name__ == '__main__':
    # test()
    main()
