#!/usr/bin/env python

import luigi
import pickle
import json
import os
import pandas as pd
from urls import WORKING_DIR
from xcms import LpcApocXcms
from apoc_inputs import ApocListPathTask


class AtomicXcmsCollection(luigi.Task):

    subset = luigi.Parameter(default="subject")

    def requires(self):
        return ApocListPathTask(self.subset)

    def output(self):
        collected = os.path.join(WORKING_DIR, self.subset + "_atmic_xcms.pkl")
        missed = os.path.join(WORKING_DIR, self.subset + "_atmic_xcms.missed")
        return luigi.LocalTarget(collected), luigi.LocalTarget(missed)

    def run(self):
        collected, missed = self.output()
        collected_content = {}
        missed_content = []
        with self.requires().output().open('r') as f:
            lines = f.read().splitlines()
            for line in lines:
                tname, qname = line.split()
                key = tname + " " + qname
                lpc_result = LpcApocXcms(tname, qname).output()
                if lpc_result.exists():
                    with lpc_result.open('r') as f:
                        collected_content[key] = json.loads(f.read())
                else:
                    missed_content.append(key)

        with collected.open('w') as f:
            pickle.dump(collected_content, f)

        with missed.open('w') as f:
            for key in missed_content:
                f.write(key + "\n")


class AtomicXcmsTable(luigi.Task):

    def requires(self):
        return AtomicXcmsCollection()

    def output(self):
        csv_path = os.path.join(WORKING_DIR, "atmic_xcms.csv")
        return luigi.LocalTarget(csv_path)

    def run(self):
        collected, missed = self.requires().output()
        with collected.open('r') as f:
            dicts = pickle.load(f)

        cols = ['Apoc ps-score', 'Kcombu tanimoto',
                'xcms', '# ligand atoms', '# residue atoms']
        datas = []
        for key, data in dicts.iteritems():
            tname, qname = key.split()
            try:
                mydata = [data[_] for _ in cols]
                mydata.append(tname)
                mydata.append(qname)
                datas.append(mydata)
            except:
                print tname, qname

        cols = cols + ['tname', 'qname']
        dset = pd.DataFrame(datas, columns=cols)
        dset.to_csv(self.output().path)


def main():
    luigi.build([AtomicXcmsCollection(),
                 AtomicXcmsTable()],
                local_scheduler=True)


if __name__ == '__main__':
    main()
