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
                            self.sdf_id,
                            self.sdf_id + ".db_scan")
        return luigi.LocalTarget(path)

    def get_data(self):
        ifn = os.path.join(self.dirname(),
                           self.sdf_id,
                           self.sdf_id + ".csv")

        dset = pd.read_csv(ifn, sep=' ', index_col=False)
        # conf_vecs = dset[['t0', 't1', 't2', 'r0', 'r1', 'r2']].values.copy()
        conf_vecs = dset[['t0', 't1', 't2']].values.copy()
        total = 10000
        if len(conf_vecs) > total:
            conf_vecs = conf_vecs[:total]

        return conf_vecs

    def run(self):
        conf_vecs = self.get_data()
        ofn = self.output().path

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


class ClusterConfs(SampleConf):

    def requires(self):
        return SampleConf(self.sdf_id)

    def output(self):
        path = os.path.join(self.dirname(),
                            self.sdf_id,
                            self.sdf_id + "_clustered.csv")
        return luigi.LocalTarget(path)

    def get_cluster_paras(self):
        ifn = SampleConf(self.sdf_id).output().path
        dset = pd.read_csv(ifn, sep=' ', index_col=False)
        dset = dset[dset['ss'] > 0.0]
        target_cluster_num = 100
        dset['Dist2Hundred'] = dset['#cluster'].apply(lambda x: abs(x - target_cluster_num))
        dset['deviedByss'] = dset['Dist2Hundred'] / dset['ss']
        sorted_dset = dset.sort(['deviedByss'])
        min_samples, eps = sorted_dset.values[0][:2]
        return min_samples, eps

    def run(self):
        min_samples, eps = self.get_cluster_paras()
        conf_vecs = self.get_data()
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(conf_vecs)

        ifn = os.path.join(self.dirname(),
                           self.sdf_id,
                           self.sdf_id + ".csv")

        dset = pd.read_csv(ifn, sep=' ', index_col=False)
        dset.ix[db.core_sample_indices_].to_csv(self.output().path,
                                                sep=' ',
                                                index=False)


def main(sdf):
    luigi.build([SampleConf(sdf), ClusterConfs(sdf)],
                local_scheduler=True)


if __name__ == '__main__':
    import sys
    sdf = sys.argv[1]
    main(sdf)
