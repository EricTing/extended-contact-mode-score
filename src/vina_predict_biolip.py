#!/usr/bin/env python

import os
import shlex
import json
import luigi
import pybel
import sample_confs
import biolip_query_biolip
import subprocess32
import dockedpose

from biolip import BioLipReferencedSpearmanR
from biolip import FixedPocketBioLipReferencedSpearmanR
from restore_order import getPredConfs, getOrderedSdf

OUTPUT_DIR = "/ddnB/work/jaydy/working/vina_biolip/"


class VinaPredictBiolipStructure(biolip_query_biolip.Path):
    def workdir(self):
        path = os.path.join(OUTPUT_DIR, self.lig_pdb)
        try:
            os.makedirs(path)
        except:
            pass
        return path

    @property
    def lig_pdbqt(self):
        return os.path.join(self.workdir(), self.lig_pdb + '.pdbqt')

    @property
    def lig_sdf(self):
        return os.path.join(self.workdir(), self.lig_pdb + '.sdf')

    def prepareInputFiles(self):
        if (not os.path.exists(self.lig_pdbqt)) or (
                not os.path.exists(self.lig_sdf)):
            lig = pybel.readfile('pdb', self.ligPdb()).next()
            lig.write('pdbqt', self.lig_pdbqt, overwrite=True)
            lig.write('sdf', self.lig_sdf, overwrite=True)

    def output(self):
        self.prepareInputFiles()
        path = self.lig_pdbqt + '.vina.pdbqt'
        return luigi.LocalTarget(path)

    def run(self):
        self.prepareInputFiles()

        bin_path = sample_confs.Path('')
        cmds = ['perl', bin_path.ehalf_box_size_bin(), self.lig_sdf]
        stdout = subprocess32.check_output(cmds)
        box_size, x, y, z = stdout.split()

        cmd = '''%s --receptor %s --ligand %s \
        --center_x %s\
        --center_y %s\
        --center_z %s\
        --size_x %s\
        --size_y %s\
        --size_z %s\
        --cpu 1\
        --out %s
        ''' % (bin_path.vina_bin(), self.prtPdbqt, self.lig_pdbqt, x, y, z,
               box_size, box_size, box_size, self.output().path)
        print cmd
        vina_out = subprocess32.check_output(shlex.split(cmd))
        ofn = self.output().path + ".txt"
        with open(ofn, 'w') as ofs:
            ofs.write(vina_out)


class VinaRandomizeBiolipStructure(VinaPredictBiolipStructure):
    def requires(self):
        return VinaPredictBiolipStructure(self.lig_pdb)

    def output(self):
        path = self.lig_pdbqt + ".rnd.pdbqt"
        return luigi.LocalTarget(path)

    def run(self):
        self.prepareInputFiles()

        bin_path = sample_confs.Path('')
        cmds = ['perl', bin_path.ehalf_box_size_bin(), self.lig_sdf]
        stdout = subprocess32.check_output(cmds)
        box_size, x, y, z = stdout.split()

        cmd = '''%s --receptor %s --ligand %s \
        --center_x %s\
        --center_y %s\
        --center_z %s\
        --size_x %s\
        --size_y %s\
        --size_z %s\
        --randomize_only\
        --cpu 1\
        --out %s
        ''' % (bin_path.vina_bin(), self.prtPdbqt, self.lig_pdbqt, x, y, z,
               box_size, box_size, box_size, self.output().path)
        print cmd
        vina_out = subprocess32.check_output(shlex.split(cmd))
        ofn = self.output().path + ".rnd.txt"
        with open(ofn, 'w') as ofs:
            ofs.write(vina_out)


class QueryVinaResultOnBioLip(VinaPredictBiolipStructure):
    def requires(self):
        return VinaPredictBiolipStructure(self.lig_pdb)

    def __run(self, vina_task):
        prt_pdb = vina_task.prtPdb
        lig_pdbqt = vina_task.output().path
        biolip_spearmanr = BioLipReferencedSpearmanR(lig_pdbqt, prt_pdb)
        inf = float('inf')
        # maximum_search_results = 100
        result = biolip_spearmanr.calculate(maximum_search_results=inf,
                                            max_tani=0.9)
        to_write = json.dumps(result,
                              sort_keys=True,
                              indent=4,
                              separators=(',', ': '))
        return to_write

    def run(self):
        vina_task = self.requires()
        to_write = self.__run(vina_task)
        with open(self.output().path, 'w') as ofs:
            ofs.write(to_write)

    def output(self):
        path = self.requires().output().path + '.1000.json'
        return luigi.LocalTarget(path)


class QueryVinaResultOnBioLipFixedPocket(QueryVinaResultOnBioLip):
    def output(self):
        path = self.requires().output().path + '.fixed.1000.json'
        return luigi.LocalTarget(path)

    def runFixed(self, vina_task):
        prt_pdb = vina_task.prtPdb
        native_lig_path = vina_task.lig_sdf
        lig_pdbqt = vina_task.output().path
        biolip_spearmanr = FixedPocketBioLipReferencedSpearmanR(
            lig_pdbqt, prt_pdb, native_lig_path)
        inf = float('inf')
        # maximum_search_results = 100
        result = biolip_spearmanr.calculate(maximum_search_results=inf,
                                            max_tani=0.9)
        to_write = json.dumps(result,
                              sort_keys=True,
                              indent=4,
                              separators=(',', ': '))
        return to_write

    def run(self):
        vina_task = self.requires()
        to_write = self.runFixed(vina_task)
        with open(self.output().path, 'w') as ofs:
            ofs.write(to_write)


class QueryVinaRandomResultOnBioLipFixedPocket(
        QueryVinaResultOnBioLipFixedPocket):
    def requires(self):
        return VinaRandomizeBiolipStructure(self.lig_pdb)


class VinaResultAccuracy(VinaPredictBiolipStructure):
    def requires(self):
        return VinaPredictBiolipStructure(self.lig_pdb)

    def output(self):
        path = self.requires().output().path + '.actual.json'
        return luigi.LocalTarget(path)

    def caculateRMSD(self):
        vina_task = self.requires()
        predicted_pdbqt = vina_task.output().path
        predicted_mol = pybel.readfile('pdbqt', predicted_pdbqt).next()
        crystal_pdbqt = vina_task.lig_pdbqt
        crystal_mol = pybel.readfile('pdbqt', crystal_pdbqt).next()

        def rmsd(m1, m2):
            c1 = [a.coords for a in m1 if not a.OBAtom.IsHydrogen()]
            c2 = [a.coords for a in m2 if not a.OBAtom.IsHydrogen()]
            return dockedpose.rmsd(c1, c2)

        return rmsd(predicted_mol, crystal_mol)

    @property
    def predicted_sdf(self):
        vina_task = self.requires()
        predicted_pdbqt = vina_task.output().path
        predicted_sdf = os.path.splitext(predicted_pdbqt)[0] + '.sdf'
        return predicted_sdf


    def CMS(self):
        vina_task = self.requires()
        predicted_pdbqt = vina_task.output().path
        crystal_sdf = vina_task.lig_sdf
        predicted_sdf = self.predicted_sdf

        init_sdf = [l.rstrip() for l in file(crystal_sdf)]
        vina_in = [l.rstrip() for l in file(vina_task.lig_pdbqt)]
        vina_out = [l.rstrip() for l in file(predicted_pdbqt)]
        all_atom_zones = getPredConfs(vina_out)

        ordered_sdfs = []
        for atom_zone in all_atom_zones:
            ordered_sdfs.append(getOrderedSdf(init_sdf, atom_zone, vina_in))

        with open(predicted_sdf, 'w') as ofs:
            for line in ordered_sdfs[0]:
                ofs.write(line)
                ofs.write("\n")

        prt_pdb = vina_task.prtPdb
        cmd = "cms -frc --lig1 %s --prt1 %s --lig2 %s --prt2 %s" % (
            predicted_sdf, prt_pdb, crystal_sdf, prt_pdb)
        print(cmd)
        cmds = shlex.split(cmd)
        output = subprocess32.check_output(cmds)
        for line in output.splitlines():
            if "cms value" in line:
                cms = float(line.split()[-1])
            if "fraction value" in line:
                fraction = float(line.split()[-1])

        return cms, fraction

    def run(self):
        rmsd = self.caculateRMSD()
        cms, fraction = self.CMS()
        result = {'rmsd': rmsd, "cms": cms, "fraction": fraction}
        to_write = json.dumps(result,
                              sort_keys=True,
                              indent=4,
                              separators=(',', ': '))
        with self.output().open('w') as ofs:
            ofs.write(to_write)


class VinaRandomAccuracy(VinaResultAccuracy):
    def requires(self):
        return VinaRandomizeBiolipStructure(self.lig_pdb)

    def output(self):
        path = self.requires().output().path + '.actual.json'
        return luigi.LocalTarget(path)


def test():
    task = VinaResultAccuracy("3ofl_JHM_B_1.pdb")
    luigi.build([task], local_scheduler=True)


def main(name):
    luigi.build(
        [
            # VinaPredictBiolipStructure(name),
            VinaResultAccuracy(name),
            # QueryVinaResultOnBioLip(name),
            # QueryVinaResultOnBioLipFixedPocket(name),
            # VinaRandomizeBiolipStructure(name),
            # QueryVinaRandomResultOnBioLipFixedPocket(name),
            VinaRandomAccuracy(name),
        ],
        local_scheduler=True)
    pass


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
