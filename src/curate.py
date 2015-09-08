#!/usr/bin/env python

import luigi
import os
import pandas as pd
import run_control_apoc
from urls import WORKING_DIR
from collect_xcms import AtomicXcmsTable


class SuccessfulPocketList(luigi.Task):

    subset = luigi.Parameter()

    def requires(self):
        return AtomicXcmsTable(self.subset)

    def output(self):
        path = os.path.join(WORKING_DIR, self.subset + "_xcms_success.lst")
        return luigi.LocalTarget(path)

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


def test():
    luigi.build([SuccessfulPocketList(subset="subject"),
                 Curate.LpcApocResultTask(subset="subject",
                                          tname="3oz2_OZ2_A_502",
                                          qname="3ocj_PLM_A_305"),
                 Curate.PairWisePsScore(tname="3oz2_OZ2_A_502",
                                        subset="subject")],
                local_scheduler=True)


if __name__ == '__main__':
    test()
