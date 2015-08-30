#!/usr/bin/env python

import luigi
import subprocess32
import os
import re
from Bio.PDB import PDBParser
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
        tokens = self.target_id.split('_')
        pdb_id = tokens[0]
        mid_two = pdb_id[1:3]
        pdb_path = os.path.join(LPC_DIR, mid_two,
                                pdb_id, pdb_id + tokens[-2] + '.pdb')
        assert os.path.exists(pdb_path), "%s does not exist" % pdb_path
        return pdb_path

    def output(self):
        return luigi.LocalTarget(self.mypath())


class Data:
    pass


class ApocResultParer:

    def __init__(self, content):
        self._content = content
        self.has_pocket_alignment = False

        self.global_sections = re.findall(r'''>>>>>>>>>>>>>>>>>>>>>>>>>   Global alignment   <<<<<<<<<<<<<<<<<<<<<<<<<<(.*?)seconds''',
                                          content,
                                          re.DOTALL)

        self.pocket_sections = re.findall(r'''>>>>>>>>>>>>>>>>>>>>>>>>>   Pocket alignment   <<<<<<<<<<<<<<<<<<<<<<<<<<(.*?)seconds''',
                                          content,
                                          re.DOTALL)

        if len(self.pocket_sections) > 0:
            self.has_pocket_alignment = True

        self._read_global()
        self._read_pocket()

    def getContent(self):
        return self._content

    def _read_global(self):
        self.global_alignments = {}
        for string in self.global_sections:
            data = Data()
            for line in string.splitlines():
                if "Structure 1" in line:
                    data.structure_1 = line.split()[2].split('.')[0]
                if "Structure 2" in line:
                    data.structure_2 = line.split()[2].split('.')[0]
                if "TM-score" in line:
                    data.tm_score = float(line.split()[-1])
                if "RMSD" in line and "Seq identity" in line:
                    data.rmsd = float(line.split(',')[0].split()[-1])
                    data.seq_identity = float(line.split(',')[-1].split()[-1])
            key = data.structure_1 + " " + data.structure_2
            self.global_alignments[key] = data

    def _read_pocket(self):
        self.pocket_alignmets = {}
        for string in self.pocket_sections:
            data = Data()
            data.has_pocket_alignment = False
            data.template_res, data.query_res = [], []
            for line in string.splitlines():
                if "Structure 1" in line and "Pocket:" in line:
                    data.tname = line.split()[-1].split(':')[-1]
                if "Structure 2" in line and "Pocket:" in line:
                    data.qname = line.split()[-1].split(':')[-1]
                if "PS-score" in line:
                    data.ps_score = float(line.split(',')[0].split()[-1])
                if "RMSD" in line and "Seq identity" in line:
                    data.rmsd = float(line.split(',')[0].split()[-1])
                    data.seq_identity = float(line.split(',')[-1].split()[-1])
                if "*" in line and "*" == line[-1] and "******" not in line:
                    data.has_pocket_alignment = True
                    tokens = line.split()
                    data.template_chainid = tokens[1]
                    data.query_chainid = tokens[4]
                    data.template_res.append(int(tokens[2]))
                    data.query_res.append(int(tokens[5]))
            key = data.tname + " " + data.qname
            self.pocket_alignmets[key] = data

    def queryGlobal(self, structure_1, structure_2):
        return self.global_alignments[structure_1 + " " + structure_2]

    def queryPocket(self, tname, qname):
        try:
            key = tname + " " + qname
            return self.pocket_alignmets[key]
        except:
            raise KeyError, "Can not find pockets for %s and %s" % (tname, qname)


class PkcombuAtomMatchParser:
    def __init__(self, oam_fn):
        f = open(oam_fn, 'r')
        self.content = f.read()
        f.close()

        self.data = Data()

    def _readPdbMatchingSerialNums(self):
        lines = self.content.splitlines()
        list_a, list_b = [], []
        for line in lines:
            if line[0].isdigit():
                tokens = line.split()
                list_a.append(int(tokens[1]))
                list_b.append(int(tokens[5]))

        return list_a, list_b

    def getMatchingSerialNums(self):
        return self._readPdbMatchingSerialNums()


def getPdbAtomsBySerialNum(pdb_fn, serial_nums):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('x', pdb_fn)
    atoms = {atom.serial_number : atom for atom in structure.get_atoms()}
    re_ordered = []
    for num in serial_nums:
        re_ordered.append(atoms[num])

    return re_ordered


class LpcKcombuResult(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()
    subset = luigi.Parameter(default="subject")

    def my_midtwo(self):
        return self.tname[1:3]

    def mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.subset,
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
        cmds = [KCOMBU_BIN,
                "-A", inputs[0],
                "-B", inputs[1],
                "-oam", self.output().path]
        try:
            subprocess32.call(["mkdir", "-p", self.mydir()])
            subprocess32.call(cmds)
        except:
            pass

    def run(self):
        self._run_kcombu()


class LpcApocResultTask(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()
    subset = luigi.Parameter(default="subject")

    def my_midtwo(self):
        return self.tname[1:3]

    def mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.subset,
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
        return [LpcPocketPathTask(self.tname),
                LpcPocketPathTask(self.qname)]

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
        luigi.build([LpcApocResultTask(tname, qname, "control"),
                     LpcKcombuResult(qname, tname, "control"),
                     LpcKcombuResult(tname, qname, "control")],
                    local_scheduler=True)
