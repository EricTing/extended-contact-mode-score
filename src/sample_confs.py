#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics


def main(sdf):
    ifn = os.path.join("/ddnB/work/jaydy/working/astex_weak",
                       sdf,
                       sdf + ".csv")

    dset = pd.read_csv(ifn, sep=' ', index_col=False)
    # conf_vecs = dset[['t0', 't1', 't2', 'r0', 'r1', 'r2']].values.copy()
    conf_vecs = dset[['t0', 't1', 't2']].values.copy()
    total = 10000
    if len(conf_vecs) > total:
        conf_vecs = conf_vecs[:total]

    ofn = os.path.join("/ddnB/work/jaydy/working/astex_weak",
                       sdf,
                       sdf + ".db_scan")

    with open(ofn, 'w') as ofs:
        ofs.write("min_samples eps #cluster ss\n")
        for num_samples in range(1, 20, 2):
            for my_eps in np.linspace(0.01, 5, 50, endpoint=True):
                db = DBSCAN(eps=my_eps, min_samples=num_samples).fit(conf_vecs)
                try:
                    ss = metrics.silhouette_score(conf_vecs, db.labels_)
                    ofs.write("%d %f %d %f\n" % (num_samples,
                                                 my_eps,
                                                 len(db.core_sample_indices_),
                                                 ss))
                except:
                    pass


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
