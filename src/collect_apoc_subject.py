#!/usr/bin/env python


import luigi
import json
import pandas as pd
from apoc_subject_spearmanr import ApocSpearmanR


class CollectSubject(luigi.Task):
    def output(self):
        path = "../dat/apoc_subject_spearmanr.csv"
        return luigi.LocalTarget(path)

    def getNames(self):
        for line in file("../dat/subject.lst"):
            if not line.startswith("tname"):
                yield line.split()

    def run(self):
        big_dset = pd.DataFrame()
        count, total = 0, 0
        for tname, qname in self.getNames():
            total += 1
            task = ApocSpearmanR(tname, qname, "subject")
            if task.complete():
                count += 1
                json_path = task.output().path
                data = json.loads(open(json_path, 'r').read())
                df = pd.DataFrame(data.items())
                df.set_index(0, inplace=True)
                big_dset = big_dset.append(df.T)

        print("%d completes out of%d" % (count, total))
        big_dset.to_csv(self.output().path)


def main():
    luigi.build(
    [
        CollectSubject(),
    ],
        local_scheduler=True
    )
    pass

if __name__ == '__main__':
    main()
