#!/usr/bin/env python


import luigi
import pybel
import os
import sample_confs
import biolip_query_biolip
import subprocess32
import shlex


OUTPUT_DIR = "/ddnB/work/jaydy/working/vina_biolip/"


class VinaPredictBiolipStructure(biolip_query_biolip.Path):

    def workdir(self):
        path = os.path.join(OUTPUT_DIR, self.lig_pdb)
        try:
            os.makedirs(path)
        except:
            pass
        return path

    def __prepareInputFiles(self):
        # ligand pdbqt
        self.lig_pdbqt = os.path.join(self.workdir(), self.lig_pdb + '.pdbqt')
        self.lig_sdf = os.path.join(self.workdir(), self.lig_pdb + '.sdf')
        lig = pybel.readfile('pdb', self.ligPdb()).next()
        if not os.path.exists(self.lig_pdbqt):
            lig.write('pdbqt', self.lig_pdbqt, overwrite=True)
        if not os.path.exists(self.lig_sdf):
            lig.write('sdf', self.lig_sdf, overwrite=True)

    def output(self):
        self.__prepareInputFiles()
        path = self.lig_pdbqt + '.vina.pdbqt'
        return luigi.LocalTarget(path)

    def run(self):
        self.__prepareInputFiles()

        bin_path = sample_confs.Path('')
        cmds = ['perl',
                bin_path.ehalf_box_size_bin(),
                self.lig_sdf]
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
        ''' % (bin_path.vina_bin(),
               self.prtPdbqt,
               self.lig_pdbqt,
               x, y, z,
               box_size,
               box_size,
               box_size,
               self.output().path
               )
        print cmd
        vina_out = subprocess32.check_output(shlex.split(cmd))
        ofn = self.output().path + ".txt"
        with open(ofn, 'w') as ofs:
            ofs.write(vina_out)


def main():
    luigi.build([
        VinaPredictBiolipStructure('3e3t_I3C_A_3.pdb'),
        VinaPredictBiolipStructure('1i9m_INW_A_2.pdb'),
    ], local_scheduler=True
    )
    pass

if __name__ == '__main__':
    main()
