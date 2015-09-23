#!/usr/bin/env python

import luigi
import subprocess32
import os
import re
import urllib
from Bio.PDB import PDBParser
from urls import LPC_DIR, APOC_WORKING_DIR, APOC_BIN, KCOMBU_BIN, PREPARE_MOL2_BIN
from apoc_inputs import LigandExpStructureInPdb


def run_cmd(cmd):
    p = subprocess32.Popen(cmd, shell=False,
                           universal_newlines=True, stdout=subprocess32.PIPE)
    p.wait()
    output = p.communicate()[0]
    return output


class LpcPocketPathTask(luigi.Task):

    target_id = luigi.Parameter()

    def _mypath(self):
        tokens = self.target_id.split('_')
        pdb_id = tokens[0]
        mid_two = pdb_id[1:3]
        pdb_path = os.path.join(LPC_DIR, mid_two,
                                pdb_id, pdb_id + tokens[-2] + '.pdb')
        assert os.path.exists(pdb_path), "%s does not exist" % pdb_path
        return pdb_path

    def output(self):
        return luigi.LocalTarget(self._mypath())


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
                    data.p_value = float(line.split(',')[1].split()[-1])
                    data.z_score = float(line.split(',')[2].split()[-1])
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

    def queryGlobal(self, tname, qname):
        structure_1 = tname.split('_')[0] + tname.split('_')[2]
        structure_2 = qname.split('_')[0] + qname.split('_')[2]
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

        for line in self.content.splitlines():
            if "tanimoto" in line:
                self.data.tanimoto = float(line.split()[-1])

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


class PdbMol2(luigi.Task):
    pdb_id = luigi.Parameter()

    def _mypath(self):
        mid_two = self.pdb_id[1:3]
        dirname = os.path.join("/ddnB/work/jaydy/dat/pdb_mol2", mid_two)
        os.system("mkdir -p %s" % (dirname))
        pdb_path = os.path.join(dirname, "%s.cif" % self.pdb_id)
        return pdb_path

    def output(self):
        return luigi.LocalTarget(self._mypath())

    def run(self):
        print "downloading %s from protein data bank" % self.pdb_id
        url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=cif&compression=NO&structureId=%s" % self.pdb_id
        urllib.urlretrieve(url, self._mypath())


class LigandExpStructureInMol2(luigi.Task):

    ligand_code = luigi.Parameter()

    def _my_pdbid(self):
        return self.ligand_code.split('_')[0]

    def _my_resn(self):
        return self.ligand_code.split('_')[1]

    def _my_chain(self):
        return self.ligand_code.split('_')[2]

    def _my_num(self):
        return self.ligand_code.split('_')[3]

    def _my_midtwo(self):
        return self._my_pdbid()[1:3]

    def _my_dir(self):
        return os.path.join(APOC_WORKING_DIR,
                            "mol2",
                            self._my_midtwo(),
                            self._my_pdbid())

    def requires(self):
        return PdbMol2(self._my_pdbid())

    def output(self):
        path = os.path.join(self._my_dir(),
                            self.ligand_code + '.mol2')
        return luigi.LocalTarget(path)

    def run(self):
        try:
            os.makedirs(self._my_dir())
        except:
            pass

        ifn = PdbMol2(self._my_pdbid()).output().path
        ofn = self.output().path
        cmds = [PREPARE_MOL2_BIN,
                "--ifn", ifn,
                "--ofn", ofn,
                "--code", self.ligand_code]
        subprocess32.call(cmds)


class LpcKcombuResult(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()
    subset = luigi.Parameter(default="subject")

    def _my_midtwo(self):
        return self.tname[1:3]

    def _mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.subset,
                                self._my_midtwo())
        return dir_path

    def _mypath(self):
        path = os.path.join(self._mydir(),
                            self.tname + "__" + self.qname + ".oam")
        return path

    def output(self):
        return luigi.LocalTarget(self._mypath())

    def _convert(self, name):
        tokens = name.split('_')
        return tokens[:3]

    def requires(self):
        return [LigandExpStructureInMol2(self.tname),
                LigandExpStructureInMol2(self.qname)]

    def _run_kcombu(self):
        try:
            os.makedirs(self._mydir())
        except:
            pass
        inputs = [_.output().path for _ in self.requires()]
        cmds = [KCOMBU_BIN,
                "-A", inputs[0],
                "-B", inputs[1],
                "-oam", self.output().path]
        try:
            subprocess32.call(["mkdir", "-p", self._mydir()])
            subprocess32.call(cmds)
        except:
            pass

    def run(self):
        self._run_kcombu()


class LpcApocResultTask(luigi.Task):

    tname = luigi.Parameter()
    qname = luigi.Parameter()
    subset = luigi.Parameter(default="subject")

    def _my_midtwo(self):
        return self.tname[1:3]

    def _mydir(self):
        dir_path = os.path.join(APOC_WORKING_DIR,
                                self.subset,
                                self._my_midtwo())
        return dir_path

    def _mypath(self):
        path = os.path.join(self._mydir(),
                            self.tname + "__" + self.qname)
        return path

    def output(self):
        return luigi.LocalTarget(self._mypath())

    def _convert(self, name):
        tokens = name.split('_')
        return tokens[0] + tokens[2]

    def requires(self):
        return [LpcPocketPathTask(self.tname),
                LpcPocketPathTask(self.qname)]

    def run_apoc(self):
        subprocess32.call(["mkdir", "-p", self._mydir()])
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
