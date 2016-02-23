#!/usr/bin/env python

from __future__ import print_function

import os
import re
import luigi
import pybel
from vina_predict_biolip import VinaResultAccuracy
from vina_predict_biolip import QueryVinaResultOnBioLipFixedPocket


def addLig2Prt(lig_sdf, prt_pdb, concated_ofn):
    """append the ligand to the protein
    """
    lig = pybel.readfile('sdf', lig_sdf).next()
    pattern = r"HETATM.*\n|ATOM  .*\n"
    atom_lines = re.findall(pattern, lig.write(format='pdb'))
    HETATM_lines = [line.replace("ATOM  ", "HETATM") for line in atom_lines]

    pdb_lines = [_ for _ in file(prt_pdb)]
    outputs = pdb_lines[:-2] + HETATM_lines + pdb_lines[-2:]
    with open(concated_ofn, 'w') as f:
        f.writelines(outputs)


class SuperPose(VinaResultAccuracy):
    """superpose the predicted protein-ligand complex onto the template's protein
    """

    def requires(self):
        """
        VinaResultAccuracy provides predicted SDF file,
        QueryVinaResultOnBioLipFixedPocket provides aligned residue numbers and ligand numbers
        """
        return [VinaResultAccuracy(self.lig_pdb),
                QueryVinaResultOnBioLipFixedPocket(self.lig_pdb)]

    def output(self):
        pass

    def append_ligand(self):
        """applend the predicted ligand to the native protein
        """
        accuracy_task = VinaResultAccuracy(self.lig_pdb)
        prt_pdb = accuracy_task.prtPdb
        predicted_sdf = accuracy_task.predicted_sdf
        concated_pdb = os.path.join(accuracy_task.workdir(),
                                    self.lig_pdb + '.concated.pdb')
        addLig2Prt(predicted_sdf, prt_pdb, concated_pdb)

    def run(self):
        self.append_ligand()


def test():
    """testing
    """
    luigi.build([SuperPose('2yiw_YIW_A_1.pdb')], local_scheduler=True)


def main():
    """launch the job
    """
    pass


if __name__ == '__main__':
    # main()
    test()
