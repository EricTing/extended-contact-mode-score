#!/usr/bin/env python

import luigi
import os
import urllib

WORKING_DIR = "/ddnB/work/jaydy/working/xcms"
DAT_DIR = "/home/jaydy/Workspace/Bitbucket/xcms/dat"


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
        subset_lst_path = "%s/%s.lst" % (DAT_DIR, self.subset)
        return luigi.LocalTarget(subset_lst_path)


class ApocPdbPathsTask(luigi.Task):

    def pdb_ids(self):

        ids = []
        with ApocListPathTask("subject").output().open('r') as lst_f:
            for line in lst_f:
                if 'tname' not in line:
                    tname, qname = line.split()
                    tpdb = tname.split('_')[0]
                    qpdb = qname.split('_')[0]
                    ids.append(tpdb)
                    ids.append(qpdb)
        with ApocListPathTask("control").output().open('r') as lst_f:
            for line in lst_f:
                if 'tname' not in line:
                    tname, qname = line.split()
                    tpdb = tname.split('_')[0]
                    qpdb = qname.split('_')[0]
                    ids.append(tpdb)
                    ids.append(qpdb)
        ids = set(ids)
        return list(ids)

    def output(self):
        output_path = os.path.join(WORKING_DIR,
                                   "apoc_pdbs_paths.txt")
        return luigi.LocalTarget(output_path)

    def requires(self):
        return ApocListPathTask("subject"), ApocListPathTask("control"),\
            [PdbPathTask(pdb_id) for pdb_id in self.pdb_ids()]

    def run(self):
        paths = [PdbPathTask(pdb_id).output().path
                 for pdb_id in self.pdb_ids()]
        with self.output().open('w') as f:
            for path in paths:
                f.write(path + "\n")


if __name__ == '__main__':
    luigi.build([ApocPdbPathsTask()],
                local_scheduler=True)
