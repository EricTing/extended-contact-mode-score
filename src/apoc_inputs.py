#!/usr/bin/env python

import luigi
import os
import urllib


class PdbPathTask(luigi.Task):

    pdb_id = luigi.Parameter()

    def mypath(self):
        mid_two = self.pdb_id[1:3]
        pdb_path = os.path.join("/ddnB/work/jaydy/dat/pdb",
                                mid_two,
                                "pdb%s.ent.gz" % self.pdb_id)
        return pdb_path

    def output(self):
        return luigi.LocalTarget(self.mypath())

    def run(self):
        print "downloading %s from protein data bank" % self.pdb_id
        url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % self.pdb_id
        urllib.urlretrieve(url, self.mypath())


class ApocListPathTask(luigi.Task):

    subset = luigi.Parameter()

    def output(self):
        subset_lst_path = "../dat/%s.lst" % self.subset
        return luigi.LocalTarget(subset_lst_path)


class BuildInputs4ApocPipeline(luigi.Task):
    subset = luigi.Parameter(default="control")

    def pdb_ids(self):
        ids = []
        with ApocListPathTask(self.subset).output().open('r') as lst_f:
            for line in lst_f:
                if 'tname' not in line:
                    tname, qname = line.split()
                    tpdb = tname.split('_')[0]
                    qpdb = qname.split('_')[0]
                    ids.append(tpdb)
                    ids.append(qpdb)
        return ids

    def output(self):
        return luigi.LocalTarget("../dat/%s_apoc_inputs.txt" % self.subset)

    def requires(self):
        return ApocListPathTask(self.subset), \
            [PdbPathTask(pdb_id) for pdb_id in self.pdb_ids()]

    def run(self):
        with ApocListPathTask(self.subset).output().open('r') as lst_f:
            with self.output().open('w') as apoc_inputs_f:
                for line in lst_f:
                    if 'tname' not in line:
                        tname, qname = line.split()
                        tpdb = tname.split('_')[0]
                        qpdb = qname.split('_')[0]
                        tpdb_path = PdbPathTask(tpdb).output().path
                        qpdb_path = PdbPathTask(qpdb).output().path
                        apoc_inputs_f.write("%s %s\n" % (tpdb_path, qpdb_path))


if __name__ == '__main__':
    luigi.build([BuildInputs4ApocPipeline("control"),
                 BuildInputs4ApocPipeline("subject")],
                local_scheduler=True)
