#!/usr/bin/env python

import luigi
import pybel
import os
import sample_confs
import biolip_query_biolip
import subprocess32
import shlex
import json
import dockedpose

from biolip import BioLipReferencedSpearmanR
from biolip import FixedPocketBioLipReferencedSpearmanR

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

    def __prepareInputFiles(self):
        if (not os.path.exists(self.lig_pdbqt)) or (
                not os.path.exists(self.lig_sdf)):
            lig = pybel.readfile('pdb', self.ligPdb()).next()
            lig.write('pdbqt', self.lig_pdbqt, overwrite=True)
            lig.write('sdf', self.lig_sdf, overwrite=True)

    def output(self):
        self.__prepareInputFiles()
        path = self.lig_pdbqt + '.vina.pdbqt'
        return luigi.LocalTarget(path)

    def run(self):
        self.__prepareInputFiles()

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
        path = self.requires().output().path + '.json'
        return luigi.LocalTarget(path)


class QueryVinaResultOnBioLipFixedPocket(QueryVinaResultOnBioLip):
    def output(self):
        path = self.requires().output().path + '.fixed.json'
        return luigi.LocalTarget(path)

    def __runFixed(self, vina_task):
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
        to_write = self.__runFixed(vina_task)
        with open(self.output().path, 'w') as ofs:
            ofs.write(to_write)


class VinaResultAccuracy(VinaPredictBiolipStructure):
    def requires(self):
        return VinaPredictBiolipStructure(self.lig_pdb)

    def output(self):
        path = self.requires().output().path + '.actual.json'
        return luigi.LocalTarget(path)

    def __calculateRmsd(self):
        vina_task = VinaPredictBiolipStructure(self.lig_pdb)
        predicted_pdbqt = vina_task.output().path
        predicted_mol = pybel.readfile('pdbqt', predicted_pdbqt).next()
        crystal_pdbqt = vina_task.lig_pdbqt
        crystal_mol = pybel.readfile('pdbqt', crystal_pdbqt).next()

        def __rmsd(m1, m2):
            c1 = [a.coords for a in m1 if not a.OBAtom.IsHydrogen()]
            c2 = [a.coords for a in m2 if not a.OBAtom.IsHydrogen()]
            return dockedpose.rmsd(c1, c2)

        return __rmsd(predicted_mol, crystal_mol)

    def __cms(self):
        vina_task = VinaPredictBiolipStructure(self.lig_pdb)
        predicted_pdbqt = vina_task.output().path
        crystal_sdf = vina_task.lig_sdf
        predicted_sdf = os.path.splitext(predicted_pdbqt)[0] + '.sdf'
        if not os.path.exists(predicted_sdf):
            # use the first predicted model
            predicted_mol = next(pybel.readfile('pdbqt', predicted_pdbqt))
            predicted_mol.write('sdf', predicted_sdf)

        prt_pdb = vina_task.prtPdb
        cmds = shlex.split("cms -frc --lig1 %s --prt1 %s --lig2 %s --prt2 %s" %
                           (predicted_sdf, prt_pdb, crystal_sdf, prt_pdb))
        print(cmds)
        subprocess32.call(cmds)

    def run(self):
        rmsd = self.__calculateRmsd()
        result = {'rmsd': rmsd}
        to_write = json.dumps(result,
                              sort_keys=True,
                              indent=4,
                              separators=(',', ': '))
        with self.output().open('w') as ofs:
            ofs.write(to_write)


def test():
    task = VinaResultAccuracy("3ofl_JHM_B_1.pdb")
    luigi.build([task], local_scheduler=True)


def main(name):
    luigi.build(
        [
            # VinaPredictBiolipStructure(name),
            VinaResultAccuracy(name),
            QueryVinaResultOnBioLip(name),
            QueryVinaResultOnBioLipFixedPocket(name),
        ],
        local_scheduler=True)
    pass


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
