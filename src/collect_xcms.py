#!/usr/bin/env python

import luigi
import pickle
import json
import os
import pandas as pd
from myurls import WORKING_DIR
from xcms import LpcApocXcms
from apoc_inputs import ApocListPathTask
from run_control_apoc import LpcApocResultTask, ApocResultParer
from mixed_resolution_xcms import MixedResolutionXcms
from disasembling import ApocRandomSet2


class ApocPocketsResultsCollection(luigi.Task):

    subset = luigi.Parameter(default="rs2")

    def requires(self):
        return ApocListPathTask(self.subset)

    def output(self):
        path = os.path.join(WORKING_DIR, self.subset + "_ps_score.pkl")
        return luigi.LocalTarget(path)

    def run(self):
        results = {}
        with self.requires().output().open('r') as f:
            for line in f.read().splitlines():
                tname, qname = line.split()
                key = tname + " " + qname
                apoc_result = LpcApocResultTask(tname, qname, self.subset)
                try:
                    with apoc_result.output().open('r') as apoc_f:
                        parser = ApocResultParer(apoc_f.read())
                        pocket_alignment = parser.queryPocket(tname, qname)
                        global_alignment = parser.queryGlobal(tname, qname)
                        tm_score = global_alignment.tm_score
                        ps_score = pocket_alignment.ps_score
                        results[key] = {"ps_score": ps_score,
                                        "tm_score": tm_score}
                except:
                    results[key] = None

        to_write = json.dumps(results, sort_keys=True,
                              indent=4, separators=(',', ': '))
        with self.output().open('w') as f:
            pickle.dump(to_write, f)


class ApocPocketsResultsTable(luigi.Task):

    subset = luigi.Parameter()

    def requires(self):
        return ApocPocketsResultsCollection(self.subset)

    def output(self):
        path = os.path.join(WORKING_DIR, self.subset + "_apoc_scores.csv")
        return luigi.LocalTarget(path)

    def _run(self):
        with self.requires().output().open('r') as inputObj:
            data = json.loads(pickle.load(inputObj))
            dset = pd.DataFrame(data)
            dset = dset.T
            dset.to_csv(self.output().path)

    def run(self):
        self._run()


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
                try:
                    tname, qname = line.split()
                    key = tname + " " + qname
                    lpc_result = LpcApocXcms(tname, qname, self.subset).output()
                    if lpc_result.exists():
                        with lpc_result.open('r') as f:
                            collected_content[key] = json.loads(f.read())
                    else:
                        missed_content.append(key)
                except ValueError:
                    print "ERROR!", line
                    pass

        with collected.open('w') as f:
            pickle.dump(collected_content, f)

        with missed.open('w') as f:
            for key in missed_content:
                f.write(key + "\n")


class AtomicXcmsTable(luigi.Task):

    subset = luigi.Parameter()

    def requires(self):
        return AtomicXcmsCollection(self.subset)

    def output(self):
        csv_path = os.path.join(WORKING_DIR, self.subset + "_atmic_xcms.csv")
        return luigi.LocalTarget(csv_path)

    def run(self):
        collected, missed = self.requires().output()
        with collected.open('r') as f:
            dicts = pickle.load(f)

        cols = ['Apoc ps-score', 'Kcombu tanimoto',
                'xcms', '# ligand atoms', '# residues']
        datas = []
        # TODO: calculate the fraction of contact
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

    def readDset(self):
        with self.output().open('r') as inputObj:
            return pd.read_csv(inputObj, index_col=0)


class MixedResolutionXcmsCollection(AtomicXcmsCollection):

    def output(self):
        collected = os.path.join(WORKING_DIR, self.subset + "_mixed_xcms.pkl")
        missed = os.path.join(WORKING_DIR, self.subset + "_mixed_xcms.missed")
        return luigi.LocalTarget(collected), luigi.LocalTarget(missed)

    def run(self):
        collected, missed = self.output()
        collected_content = {}
        missed_content = []
        with self.requires().output().open('r') as f:
            lines = f.read().splitlines()
            count = 0
            for line in lines:
                count += 1
                print("%s %f%%\n" % (line, count  * 100. / len(lines)))
                try:
                    tname, qname = line.split()
                    key = tname + " " + qname
                    lpc_result = MixedResolutionXcms(tname, qname, self.subset).output()
                    if lpc_result.exists():
                        with lpc_result.open('r') as f:
                            collected_content[key] = json.loads(f.read())
                    else:
                        missed_content.append(key)
                except ValueError:
                    print "ERROR!", line
                    pass

        with collected.open('w') as f:
            pickle.dump(collected_content, f)

        with missed.open('w') as f:
            for key in missed_content:
                f.write(key + "\n")


class MixedResolutionXcmsTable(AtomicXcmsTable):

    def requires(self):
        return MixedResolutionXcmsCollection(self.subset)

    def output(self):
        csv_path = os.path.join(WORKING_DIR, self.subset + "_mixed_xcms.csv")
        return luigi.LocalTarget(csv_path)


class CollectXcmsApocRS2(luigi.Task):

    def requires(self):
        return ApocListPathTask("rs2")

    def run(self):
        rs2_lst = self.requires().output().path
        data = []
        for tname, qname in [line.split() for line in file(rs2_lst)]:
            rs2_task = ApocRandomSet2(tname, qname)
            result_path = rs2_task.output().path
            if os.path.exists(result_path):
                with open(result_path, 'r') as ifs:
                    cms_line = next((l for l in ifs.readlines()
                                     if 'Extended CMS' in l),
                                    None)
                    if cms_line is not None:
                        key = tname + "__" + qname
                        val = float(cms_line.split()[-1])
                        data.append((key, val))
        dset = pd.DataFrame(data)
        dset.to_csv(self.output().path)

    def output(self):
        path = "../dat/rs2_xcms.csv"
        return luigi.LocalTarget(path)


def main():
    luigi.build([CollectXcmsApocRS2()],
                local_scheduler=True)


if __name__ == '__main__':
    main()
