#!/usr/bin/env python

from __future__ import print_function

import os
import pybel
import json
import luigi
import subprocess32
from collections import defaultdict

from run_control_apoc import PkcombuAtomMatchParser
from astex_xcms import runPkcombu


class CheckNames(luigi.Task):
    def output(self):
        path = "../dat/dude_jobs.txt"
        return luigi.LocalTarget(path)

    def run(self):
        names = []
        with open("../dat/protein-name-match.dat", 'rb') as ifs:
            for line in ifs:
                lig, prt = line.split()
                try:
                    assert lig == prt[:-1]
                    names.append(prt + "\n")
                except Exception as inst:
                    print(inst)

        with open(self.output().path, 'w') as ofs:
            ofs.writelines(names)


class Path(luigi.Task):

    prt_id = luigi.Parameter()

    def getPrtId(self):
        return self.prt_id

    def getLigId(self):
        return self.prt_id[:-1]

    def getCrystalLigDir(self):
        return "/ddnB/work/jaydy/dat/edud_vina/EDUD-PDBQT/crystal-ligands-pdbqt"

    def getCrystalPrtDir(self):
        return "/ddnB/work/jaydy/dat/edud_vina/EDUD-PDBQT/crystal-proteins-pdbqt"

    def getCrystalLigPdbqt(self):
        return os.path.join(self.getCrystalLigDir(),
                            self.getLigId() + ".pdbqt")

    def getCrystalPrtPdbqt(self):
        return os.path.join(self.getCrystalPrtDir(),
                            self.getPrtId() + ".pdbqt")


class PrtPdb(Path):

    def datDir(self):
        return "/ddnB/work/jaydy/dat/edud_vina/EDUD-PDBQT/crystal-proteins-pdb"

    def output(self):
        path = os.path.join(self.datDir(),
                            self.getPrtId() + ".pdb")
        return luigi.LocalTarget(path)

    def getPdbPath(self):
        if not os.path.exists(self.output().path):
            self.run()
        return self.output().path

    def run(self):
        mol = pybel.readfile("pdbqt", self.getCrystalPrtPdbqt()).next()
        mol.write(format="pdb", filename=self.output().path)


class LigSdf(Path):

    def datDir(self):
        return "/ddnB/work/jaydy/dat/edud_vina/EDUD-PDBQT/crystal-ligands-pdbqt"

    def output(self):
        path = os.path.splitext(self.getCrystalLigPdbqt())[0] + '.sdf'
        return luigi.LocalTarget(path)

    def getSdfPath(self):
        if not os.path.exists(self.output().path):
            self.run()
        return self.output().path

    def run(self):
        mol = pybel.readfile("pdbqt", self.getCrystalLigPdbqt()).next()
        mol.removeh()
        mol.write(format="sdf", filename=self.output().path, overwrite=True)


class Calculate(luigi.Task):
    tar_name = luigi.Parameter()
    subset = luigi.Parameter()

    def output(self):
        try:
            os.makedirs(os.path.join("/work/jaydy/working", self.subset))
        except:
            pass
        path = os.path.join("/work/jaydy/working",
                            self.subset, self.tar_name + '.json')
        return luigi.LocalTarget(path)

    def run(self):
        def untar():
            my_dir = '/var/scratch/' + os.getenv('USER')
            try:
                os.makedirs(my_dir)
            except:
                pass

            targz_path = os.path.join("/work/jaydy/dat/edud_vina",
                                      self.subset, self.tar_name + '.tar.gz')
            cmds = [
                'tar', 'xf', targz_path,
                '-C', my_dir
            ]
            subprocess32.call(cmds)
            return os.path.join(my_dir, self.tar_name)

        def getPdbqtFiles(root_dir):
            paths = []
            for (dirname, subshere, fileshere) in os.walk(root_dir):
                for fname in fileshere:
                    if fname.endswith('pdbqt'):
                        path = os.path.join(dirname, fname)
                        paths.append(path)
            return paths

        results = defaultdict(dict)

        untared_dir = untar()
        pdbqt_fns = getPdbqtFiles(untared_dir)

        def runXcms():
            for fn in pdbqt_fns:
                lig_sdf = os.path.splitext(fn)[0] + '.sdf'

                if not os.path.exists(lig_sdf):
                    lig = pybel.readfile("pdbqt", fn).next()
                    lig.removeh()
                    lig.write(format="sdf", filename=lig_sdf)

                prt_id = os.path.split(fn)[-1].split('-')[0]
                crystal_lig_sdf = LigSdf(prt_id).getSdfPath()

                oam_path = os.path.splitext(fn)[0] + '.oam'

                runPkcombu(lig_sdf, crystal_lig_sdf, oam_path)

                try:
                    lml_path = os.path.splitext(fn)[0] + '.lml'
                    kcombu = PkcombuAtomMatchParser(oam_path)
                    kcombu.writeMatchingSerialNums(lml_path)
                    tc = kcombu.getTc()
                    list_a, list_b = kcombu.getMatchingSerialNums()

                    crystal_prt_pdb = PrtPdb(prt_id).getPdbPath()
                    cmds = ['run_xcms',
                            '--la', lig_sdf,
                            '--lb', crystal_lig_sdf,
                            '--pa', crystal_prt_pdb,
                            '--pb', crystal_prt_pdb,
                            '--lml', lml_path]

                    lig_id = os.path.basename(lig_sdf)
                    results[prt_id][lig_id] = {
                        'Tc': tc,
                        'dude_ligand_indices': list_a,
                        'crystal_ligand_indices': list_b,
                    }

                    try:
                        xcms_result = subprocess32.check_output(cmds)
                        for line in xcms_result.splitlines():
                            if 'PS-score' in line:
                                ps_score = float(line.split(':')[-1])
                                results[prt_id][lig_id].update(
                                    {'ps-score': ps_score})
                            if 'residue list of the first protein' in line:
                                prt_contact_indices = map(
                                    int, line.split(':')[-1].split())
                                results[prt_id][lig_id].update(
                                    {'dude_contact_indices': prt_contact_indices})
                            if 'Extended CMS' in line:
                                xcms_val = float(line.split()[-1])
                                results[prt_id][lig_id].update({'xcms': xcms_val})
                    except:
                        print(" ".join(cmds))
                        results[prt_id][lig_id].update({'xcms': None})
                        pass

                except Exception as inst:
                    print(inst)
                    pass

        runXcms()
        to_write = json.dumps(results, sort_keys=True,
                              indent=4, separators=(',', ': '))
        with open(self.output().path, 'a') as ofs:
            ofs.write(to_write)


def main(task_name):
    luigi.build([
        Calculate(task_name,
                  'output-crystal-active-mod-optimal'),
    ],
                local_scheduler=True
    )
    pass


if __name__ == '__main__':
    import sys
    main(sys.argv[1])
