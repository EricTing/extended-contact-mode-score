#!/usr/bin/env python

import luigi
import subprocess32
import os
from urls import LPC_DIR, APOC_WORKING_DIR, APOC_BIN, KCOMBU_BIN
from apoc_inputs import LigandExpStructureInPdb


def run_cmd(cmd):
    p = subprocess32.Popen(cmd, shell=False,
                           universal_newlines=True, stdout=subprocess32.PIPE)
    p.wait()
    output = p.communicate()[0]
    return output


class LpcPocketPathTask(luigi.Task):

    target_id = luigi.Parameter()

    def mypath(self):
        mid_two = self.target_id[1:3]
        pdb_path = os.path.join(LPC_DIR, mid_two,
                                self.target_id[:-1], self.target_id + '.pdb')
        return pdb_path

    def output(self):
        return luigi.LocalTarget(self.mypath())


class Data:
    pass


class ApocResultParer:

    def __init__(self, result):
        '''
        result: str of list of strs
        '''
        self.result = result.splitlines(False)
        self.global_property = Data()
        self.pocket_property = Data()
        self.matching_list = []

        global_idx, pocket_idx = 0, 0
        for idx, line in enumerate(self.result):
            if "Global alignment" in line:
                global_idx = idx
            if "Pocket alignment" in line:
                pocket_idx = idx

        for idx, line in enumerate(self.result):
            if idx > global_idx and idx < pocket_idx:
                if "TM-score" in line:
                    self.global_property.tm_score = float(line.split()[-1])
                if "RMSD" in line and "Seq identity" in line:
                    self.global_property.rmsd = float(line.split(',')[0].split()[-1])
                    self.global_property.seq_identity  = float(line.split(',')[-1].split()[-1])
            if idx > pocket_idx:
                if "TM-score" in line:
                    self.pocket_property.tm_score = float(line.split()[-1])
                if "RMSD" in line and "Seq identity" in line:
                    self.pocket_property.rmsd = float(line.split(',')[0].split()[-1])
                    self.pocket_property.seq_identity  = float(line.split(',')[-1].split()[-1])
                if "PS-score" in line:
                    self.pocket_property.ps_score = float(line.split(',')[0].split()[-1])
                if "*" in line and "*" == line[-1] and "******" not in line:
                    tokens = line.split()
                    self.matching_list.append([[tokens[1], int(tokens[2])],
                                               [tokens[4], int(tokens[5])]])


class LpcKcombuResult(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()

    def my_midtwo(self):
        return self.tname[1:3]

    def mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.my_midtwo())
        return dir_path

    def mypath(self):
        path = os.path.join(self.mydir(),
                            self.tname + "__" + self.qname + ".oam")
        return path

    def output(self):
        return luigi.LocalTarget(self.mypath())

    def convert(self, name):
        tokens = name.split('_')
        return tokens[:3]

    def requires(self):
        t_tokens = self.convert(self.tname)
        q_tokens = self.convert(self.qname)
        return LigandExpStructureInPdb(t_tokens[0], t_tokens[1], t_tokens[2]),\
            LigandExpStructureInPdb(q_tokens[0], q_tokens[1], q_tokens[2])

    def _run_kcombu(self):
        inputs = [_.output().path for _ in self.requires()]
        args = [KCOMBU_BIN,
                "-A", inputs[0],
                "-B", inputs[1],
                "-oam", self.output().path]
        output = "No Matches!\n"
        try:
            output = subprocess32.check_output(args)
        except:
            with self.output().open('w') as f:
                f.write(output)

    def run(self):
        self._run_kcombu()


class LpcApocResultTask(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()

    def my_midtwo(self):
        return self.tname[1:3]

    def mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.my_midtwo())
        return dir_path

    def mypath(self):
        path = os.path.join(self.mydir(),
                            self.tname + "__" + self.qname)
        return path

    def output(self):
        return luigi.LocalTarget(self.mypath())

    def convert(self, name):
        tokens = name.split('_')
        return tokens[0] + tokens[2]

    def requires(self):
        return [LpcPocketPathTask(self.convert(self.tname)),
                LpcPocketPathTask(self.convert(self.qname))]

    def run_apoc(self):
        subprocess32.call(["mkdir", "-p", self.mydir()])
        paths = [_.output().path for _ in self.requires()]
        cmd = [APOC_BIN] + paths
        print " ".join(cmd)
        stdout = subprocess32.check_output(cmd)
        return stdout

    def run(self):
        stdout = self.run_apoc()
        with self.output().open('w') as f:
            f.write(stdout)


if __name__ == "__main__":
    import sys
    tname, qname = sys.argv[1], sys.argv[2]
    if tname != "tname":
        luigi.build([LpcApocResultTask(tname, qname),
                     LpcKcombuResult(tname, qname)],
                    local_scheduler=True)
