#!/usr/bin/env python

import luigi
import pickle
import json
import os
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
                output = LpcApocXcms(tname, qname).output()
                if output.exists():
                    with output.open('r') as f:
                        collected_content[key] = json.loads(f.read())
                else:
                    missed_content.append(key)

        with collected.output().open('w') as f:
            pickle.dump(collected_content, f)

        with missed.output().open('w') as f:
            for key in missed_content:
                f.write(key + "\n")


def main():
    luigi.build([AtomicXcmsCollection()],
                local_scheduler=True)


if __name__ == '__main__':
    main()
