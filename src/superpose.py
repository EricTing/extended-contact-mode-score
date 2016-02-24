#!/usr/bin/env python

from __future__ import print_function

import os
import re
import json
import shlex
import shutil
import subprocess32
import luigi
import pybel
from vina_predict_biolip import VinaResultAccuracy
from vina_predict_biolip import QueryVinaResultOnBioLipFixedPocket
from biolip_query_biolip import Path
from readme import read2Df, similarPocketsLigands, clean


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

    superpose_perl = luigi.Parameter(
        default=
        "/home/jaydy/Workspace/Bitbucket/GeauxFindSite/src/superpose.pl")

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

        return the path of the concated pdb
        """
        accuracy_task = VinaResultAccuracy(self.lig_pdb)
        prt_pdb = accuracy_task.prtPdb
        predicted_sdf = accuracy_task.predicted_sdf
        concated_pdb = os.path.join(accuracy_task.workdir(),
                                    self.lig_pdb + '.concated.pdb')
        addLig2Prt(predicted_sdf, prt_pdb, concated_pdb)
        return concated_pdb

    def superpose(self):
        """superpose and copy the referenced protein and ligand
        """
        result = json.loads(QueryVinaResultOnBioLipFixedPocket(
            self.lig_pdb).output().open('r').read())
        dset = read2Df(result, self.lig_pdb)
        dset = similarPocketsLigands(clean(dset), minimum_Tc=0.4)
        work_dir = os.path.join(self.workdir(), 'superpose')
        try:
            os.makedirs(work_dir)
        except Exception as detail:
            print(detail)

        mob_pdb = self.append_ligand()
        for template_pdb in dset.index:
            ref_pdb = Path(template_pdb).prtPdb
            lig_code = Path(template_pdb).lig_code
            ref_lig = Path(lig_code).ligPdb() + '.pdb'
            sup_pdb = os.path.join(work_dir, lig_code + '.sup.pdb')
            cmd = shlex.split("perl %s all %s %s %s" %
                              (self.superpose_perl, ref_pdb, mob_pdb, sup_pdb))
            subprocess32.call(cmd)

            shutil.copy(ref_pdb, work_dir)
            shutil.copy(ref_lig, work_dir)

    def run(self):
        self.append_ligand()
        self.superpose()


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
