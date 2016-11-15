#!/usr/bin/env python


import luigi
import os
import pandas as pd
from collections import defaultdict
from dude_vina import SpearmanR, BioLipSpearmanR


class CollectDudeVinaSpearmanR(luigi.Task):
    def getTaskNames(self):
        task_names = [_.rstrip() for _
                      in file("../dat/crystal-active-jobs.txt")]
        return task_names

    def outDir(self):
        return "/ddnB/work/jaydy/working/output-crystal-active-spearmr"

    def run(self):
        task_names = self.getTaskNames()
        results = defaultdict(list)
        for name in task_names:
            task = SpearmanR(name, "output-crystal-active-mod-optimal")
            if task.complete():
                dfs = task.result2Dataframe()
                for key, df in dfs.iteritems():
                    results[key].append(df)

                print(name, "done")
            else:
                print(name, "not done")

        total_df = pd.DataFrame()
        for key, dfs in results.iteritems():
            df = pd.concat(dfs, axis=0)
            total_df = total_df.append(df)
            ofn = os.path.join(self.outDir(), key + '.csv')
            df.to_csv(ofn)

        ofn = os.path.join(self.outDir(), 'active.csv')
        total_df.to_csv(ofn)


class CollectDudeVinaBioLipSpearmanR(CollectDudeVinaSpearmanR):
    def run(self):
        pass

    def check(self):
        for task_name in self.getTaskNames():
            try:
                task = BioLipSpearmanR(task_name,
                                       'output-crystal-active-mod-optimal')
                if not task.complete():
                    print task_name, "TO-DO"
                else:
                    print task_name, "Done"
            except Exception as detail:
                print detail


def main():
    luigi.build([
        # CollectDudeVinaSpearmanR()
    ],
                local_scheduler=True)

    collect_biolip_spearmanr = CollectDudeVinaBioLipSpearmanR()
    collect_biolip_spearmanr.check()


if __name__ == '__main__':
    main()
