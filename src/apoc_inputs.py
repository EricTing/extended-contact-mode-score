#!/usr/bin/env python

import luigi
import gzip
import os
import urllib

from urls import WORKING_DIR, DAT_DIR, APOC_WORKING_DIR


class PdbPathTask(luigi.Task):

    pdb_id = luigi.Parameter()

    def _mypath(self):
        mid_two = self.pdb_id[1:3]
        pdb_path = os.path.join("/ddnB/work/jaydy/dat/pdb",
                                mid_two,
                                "pdb%s.ent.gz" % self.pdb_id)
        return pdb_path

    def output(self):
        return luigi.LocalTarget(self._mypath())

    def run(self):
        print "downloading %s from protein data bank" % self.pdb_id
        url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % self.pdb_id
        urllib.urlretrieve(url, self._mypath())


class DecompressedPdb(luigi.Task):

    pdb_id = luigi.Parameter()

    def requires(self):
        return PdbPathTask(self.pdb_id)

    def output(self):
        mid_two = self.pdb_id[1:3]
        pdb_path = os.path.join("/ddnB/work/jaydy/dat/pdb",
                                mid_two,
                                "pdb%s.ent" % self.pdb_id)
        return luigi.LocalTarget(pdb_path)

    def run(self):
        with gzip.open(self.requires().output().path, 'rb') as f:
            file_content = f.read()
            with self.output().open('w') as of:
                of.write(file_content)


class LigandExpStructureInPdb(luigi.Task):

    pdb_id = luigi.Parameter()
    lig_id = luigi.Parameter()
    chain_id = luigi.Parameter()

    def output(self):
        mid_two = self.pdb_id[1:3]
        path = os.path.join(APOC_WORKING_DIR,
                            mid_two,
                            "%s_%s_%s.pdb" % (self.pdb_id,
                                              self.lig_id,
                                              self.chain_id))
        return luigi.LocalTarget(path)

    def requires(self):
        return PdbPathTask(self.pdb_id)

    def run(self):
        with gzip.open(self.requires().output().path, 'rb') as f:
            content = f.read()
        ligand_line_mark = self.lig_id + " " + self.chain_id
        with self.output().open('w') as f:
            for line in content.splitlines():
                if ligand_line_mark in line:
                    f.write(line + "\n")


class ApocListPathTask(luigi.Task):

    subset = luigi.Parameter()

    def run(self):
        if self.subset == "subject" or self.subset == "control":
            subset_lst_path = "%s/%s.lst" % (DAT_DIR, self.subset)
            assert(os.path.exists(subset_lst_path))
            pass
        elif self.subset == "rs2":
            subset_lst_path = "/ddnB/work/jaydy/dat/apoc/RS2.lst"
            assert(os.path.exists(subset_lst_path))

            valid = self.output().open('w')
            invalid = open("/ddnB/work/jaydy/dat/apoc/RS2_invalid.lst", 'w')
            with open(subset_lst_path, 'r') as f:
                while True:
                    line = f.readline()
                    if line:
                        tokens = line.split()
                        if len(tokens) == 2:
                            valid.write(line)
                        else:
                            invalid.write(line)
                    else:
                        break
            valid.close()
            invalid.close()
        else:
            raise KeyError("cannot find the list for %s" % self.subset)

    def output(self):
        if self.subset == "subject" or self.subset == "control":
            subset_lst_path = "%s/%s.lst" % (DAT_DIR, self.subset)
            return luigi.LocalTarget(subset_lst_path)
        elif self.subset == "rs2":
            subset_lst_path = "/ddnB/work/jaydy/dat/apoc/RS2_valid.lst"
            return luigi.LocalTarget(subset_lst_path)
        else:
            raise KeyError("cannot find the list for %s" % self.subset)


class ApocPdbPathsTask(luigi.Task):

    def _pdb_ids(self):

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
            [PdbPathTask(pdb_id) for pdb_id in self._pdb_ids()]

    def run(self):
        paths = [PdbPathTask(pdb_id).output().path
                 for pdb_id in self._pdb_ids()]
        with self.output().open('w') as f:
            for path in paths:
                f.write(path + "\n")


if __name__ == '__main__':
    # luigi.build([ApocPdbPathsTask()],
    #             local_scheduler=True)
    luigi.build([LigandExpStructureInPdb('104m', 'NBN', 'A'),
                 LigandExpStructureInPdb('1m8e', '7NI', 'A')],
                local_scheduler=True)
