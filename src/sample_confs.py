#!/usr/bin/env python

import os
import pandas as pd
from sklearn.cluster import KMeans
import numpy as np
from sklearn import metrics
import luigi


class SampleConf(luigi.Task):

    sdf_id = luigi.Parameter()

    def dirname(self):
        return "/ddnB/work/jaydy/working/astex_weak"

    def requires(self):
        pass

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_kmeans.csv")
        return luigi.LocalTarget(path)

    def run(self):
        ifn = os.path.join(self.dirname(),
                           self.sdf_id,
                           self.sdf_id + ".csv")

        dset = pd.read_csv(ifn, sep=' ', index_col=False)
        total_samples = 20000
        if len(dset > total_samples):
            rows = np.random.choice(dset.index.values, total_samples)
            dset = dset.ix[rows]
        conf_vecs = dset[['t0', 't1', 't2', 'r0', 'r1', 'r2']].values.copy()

        num_clusters = 100
        kmeans = KMeans(n_clusters=num_clusters, n_jobs=16)
        kmeans.fit(conf_vecs)
        centers = kmeans.cluster_centers_

        ofn = self.output().path
        with open(ofn, 'w') as ofs:
            ofs.write("lig prt t0 t1 t2 r0 r1 r2\n")
            for center in centers:
                ofs.write("0 0 " + " ".join(map(str, center)) + "\n")


def main(sdf):
    luigi.build([SampleConf(sdf)],
                local_scheduler=True)


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
