#!/usr/bin/env python

import os
import pandas as pd
from sklearn.cluster import DBSCAN
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
                            self.sdf_id + ".csv")
        return luigi.LocalTarget(path)

    def run(self):
        ifn = os.path.join(self.dirname(),
                           self.sdf_id,
                           self.sdf_id + ".csv")

        dset = pd.read_csv(ifn, sep=' ', index_col=False)
        # conf_vecs = dset[['t0', 't1', 't2', 'r0', 'r1', 'r2']].values.copy()
        conf_vecs = dset[['t0', 't1', 't2']].values.copy()
        total = 10000
        if len(conf_vecs) > total:
            conf_vecs = conf_vecs[:total]

        ofn = os.path.join(self.dirname(),
                           self.sdf_id,
                           self.sdf_id + ".db_scan")

        with open(ofn, 'w') as ofs:
            ofs.write("min_samples eps #cluster ss\n")
            for num_samples in np.arange(1, len(conf_vecs), 50):
                for my_eps in [0.001 * 2 ** i for i in range(1, 14)]:
                    db = DBSCAN(eps=my_eps, min_samples=num_samples).fit(conf_vecs)
                    try:
                        ss = metrics.silhouette_score(conf_vecs, db.labels_)
                        ofs.write("%d %f %d %f\n" % (num_samples,
                                                     my_eps,
                                                     len(db.core_sample_indices_),
                                                     ss))
                    except:
                        pass


def main(sdf):
    luigi.build([SampleConf(sdf)],
                local_scheduler=True)


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
