#!/usr/bin/env python

import luigi
import random
import pybel
import json
import os
import scipy.cluster.hierarchy as hr
import numpy as np
from biolip import BioLipReferencedSpearmanR
from collections import defaultdict

BIOLIP_LIGS_SDF = "/ddnB/work/jaydy/dat/BioLip/ligand_nr.sdf"
OUTPUT_DIR = "/ddnB/work/jaydy/working/biolipbiolip/"


class Path(luigi.Task):
    lig_pdb = luigi.Parameter()
    prt_dir = luigi.Parameter(default="/work/jaydy/dat/BioLip/prt/")
    lig_dir = luigi.Parameter(default="/ddnB/work/jaydy/dat/BioLip/ligand_nr/")

    @property
    def pdbqt_dir(self):
        return os.path.join(
            "/ddnB/work/jaydy/dat/BioLip/prt_pdbqt/",
            self.mid_two,
        )

    @property
    def prtPdbqt(self):
        pdbqt = self.lig_pdb.split(
            '_')[0] + self.lig_pdb.split('_')[2] + ".pdbqt"
        path = os.path.join(
            "/ddnB/work/jaydy/dat/BioLip/prt_pdbqt/",
            self.mid_two,
            pdbqt
        )
        return path

    @property
    def prtPdb(self):
        code = os.path.basename(self.lig_pdb)
        mid_two = code[1:3]
        pdb = code.split('_')[0] + code.split('_')[2] + ".pdb"
        path = os.path.join(self.prt_dir, mid_two, pdb)
        return path

    def ligPdb(self):
        mid_two = self.lig_pdb[1:3]
        path = os.path.join(self.lig_dir, mid_two, self.lig_pdb)
        return path

    @property
    def mid_two(self):
        return self.lig_pdb[1:3]


class Cluster(luigi.Task):
    def requires(self):
        pass

    def run(self):
        mols = list(pybel.readfile('sdf', BIOLIP_LIGS_SDF))
        # mols = list(pybel.readfile('sdf',
        #                            "/ddnB/work/jaydy/dat/BioLip/sml_ligand_nr/a5.sdf"))

        def calMat(mols):
            fps = [mol.calcfp('maccs') for mol in mols]
            tanimoto_vals = [[myfp | otherfp for otherfp in fps]
                             for myfp in fps]
            return np.asmatrix(tanimoto_vals)

        def cluster(mols, cutoff=0.95):
            tanimoto_vals = calMat(mols)
            z = hr.linkage(tanimoto_vals, method='ward', metric='euclidean')
            cluster_ids = hr.fcluster(z, cutoff, criterion='distance')
            titles = [mol.title for mol in mols]
            clusters = defaultdict(list)
            for title, cluster_id in zip(titles, cluster_ids):
                clusters[str(cluster_id)].append(title)
            return clusters

        random.shuffle(mols)
        total_chunks = 10
        chunk_size = len(mols) / total_chunks
        chunks = [mols[i:i + chunk_size]
                  for i in xrange(0, len(mols), chunk_size)]
        results = []
        for chunk in chunks:
            results.append(cluster(chunk, cutoff=0.95))

        to_write = json.dumps(results, sort_keys=True,
                              indent=4, separators=(',', ': '))
        ofn = self.output().path
        print(ofn)
        with open(ofn, 'w') as ofs:
            ofs.write(to_write)

    def output(self):
        path = "../dat/biolip_cluster.json"
        return luigi.LocalTarget(path)


class SampleClusters(luigi.Task):
    path = luigi.Parameter(default="../dat/biolipbiolip_sampled.txt")
    total = luigi.Parameter(default=2000)

    def requires(self):
        return Cluster()

    def output(self):
        return luigi.LocalTarget(self.path)

    def run(self):
        representatives = []
        with self.requires().output().open('r') as ifs:
            n_clusters = json.load(ifs)
            for clusters in n_clusters:
                for members in clusters.values():
                    representatives.append(random.choice(members))

        if len(representatives) > self.total:
            representatives = random.sample(representatives, self.total)

        names = map(lambda name: name.split('/')[-1],
                    representatives)
        with open(self.output().path, 'w') as ofs:
            for name in names:
                ofs.write(name + "\n")


class BioLipBioLip(Path):

    def output(self):
        outdir = os.path.join(OUTPUT_DIR, self.mid_two)
        try:
            os.makedirs(outdir)
        except:
            pass
        path = os.path.join(outdir, self.lig_pdb + ".1000.json")
        return luigi.LocalTarget(path)

    def run(self):
        prt_pdb = self.prtPdb
        lig_pdb = self.ligPdb()
        biolip_spearmanr = BioLipReferencedSpearmanR(lig_pdb, prt_pdb)
        inf = float('inf')
        result = biolip_spearmanr.calculate(maximum_search_results=inf,
                                            max_tani=0.9)

        to_write = json.dumps(result, sort_keys=True,
                              indent=4, separators=(',', ': '))
        ofn = self.output().path
        print(ofn)
        with open(ofn, 'w') as ofs:
            ofs.write(to_write)

    def requires(self):
        pass


def main(name):
    luigi.build([
        BioLipBioLip(name)
    ], local_scheduler=True
    )

if __name__ == '__main__':
    import sys
    name = sys.argv[1]
    main(name)
    # luigi.build([
    #     Cluster(),
    #     SampleClusters()
    # ], local_scheduler=True
    # )
