import pandas as pd

import matplotlib
matplotlib.use('Agg')

import json
import os
from sklearn import metrics
from pprint import pprint
import minepy

import luigi

import biolip_query_biolip
import vina_predict_biolip

from vina_predict_biolip import QueryVinaResultOnBioLipFixedPocket
from vina_predict_biolip import QueryVinaRandomResultOnBioLipFixedPocket

sampled_list = "../dat/biolipbiolip_sampled_2.txt"
sampled_names = [_.rstrip() for _ in file(sampled_list)]


def read2Df(result, task_name):
    df = pd.DataFrame(dict(result)).T
    df['query'] = [task_name] * len(df)
    df['template'] = df.index
    return df


class Read(luigi.Task):
    def output(self):
        path = "/ddnB/work/jaydy/working/biolip_sampled.csv"
        return luigi.LocalTarget(path)

    def run(self):
        sampled_df = pd.DataFrame()
        for name in sampled_names:
            task = biolip_query_biolip.BioLipBioLip(name)
            if task.complete():
                with task.output().open('r') as ifs:
                    result = json.loads(ifs.read())
                    df = read2Df(result, name)
                    sampled_df = sampled_df.append(df, ignore_index=True)

        sampled_df.to_csv(self.output().path, ignore_index=True)


def cutRedundantTemplates(df):
    my_df = df.copy()
    my_df['prt_name'] = my_df['template'].apply(lambda t: t[3:7])
    my_df = my_df.groupby('prt_name').apply(lambda x: x.sample(1))
    my_df.drop('prt_name', inplace=True, axis=1)
    # my_df['lig_name'] = my_df['template'].apply(lambda t: t.split('_')[1])
    # TODO: rank instead of random sampling
    # my_df = my_df.groupby('lig_name').apply(lambda x: x.sample(1))
    return my_df


class CutRedundancy(luigi.Task):
    def requires(self):
        return Read()

    def output(self):
        path = "../dat/biolip_sampled_cutted.csv"
        return luigi.LocalTarget(path)

    def run(self):
        df = pd.read_csv(Read().output().path, index_col=0)
        cutted = df.groupby('query', as_index=False)\
                   .apply(cutRedundantTemplates)
        cutted.set_index('query', inplace=True)
        cutted.to_csv(self.output().path, ignore_index=True)


class CuttedVinaPredictBioLip(luigi.Task):
    def requires(self):
        pass

    def output(self):
        path = "/work/jaydy/working/vina_biolip_sampled.csv"
        return luigi.LocalTarget(path)

    def check(self):
        completes, incompletes = [], []
        for name in sampled_names:
            task = vina_predict_biolip.QueryVinaResultOnBioLip(name)
            if task.complete():
                completes.append(name)
            else:
                incompletes.append(name)
        print("{} completes and {} incompletes".format(
            len(completes), len(incompletes)))

        return completes

    def run(self):
        sampled_df = pd.DataFrame()
        for name in self.check():
            task = vina_predict_biolip.QueryVinaResultOnBioLip(name)
            with task.output().open('r') as ifs:
                result = json.loads(ifs.read())
                df = read2Df(result, name)
                sampled_df = sampled_df.append(df, ignore_index=True)
        cutted = cutRedundantTemplates(sampled_df)
        cutted.set_index('query', inplace=True)
        cutted.to_csv(self.output().path, ignore_index=True)


class VinaPredictBioLipFixed(luigi.Task):
    def check(self):
        completes, incompletes = [], []
        for name in sampled_names:
            task = QueryVinaResultOnBioLipFixedPocket(name)
            if task.complete():
                completes.append(name)
            else:
                incompletes.append(name)
        print("{} completes and {} incompletes".format(
            len(completes), len(incompletes)))
        return completes

    def output(self):
        path = "/work/jaydy/working/vina_biolip_sampled_fixed.csv"
        return luigi.LocalTarget(path)

    def run(self):
        sampled_df = pd.DataFrame()
        for name in self.check():
            task = QueryVinaResultOnBioLipFixedPocket(name)
            with task.output().open('r') as ifs:
                result = json.loads(ifs.read())
                df = read2Df(result, name)
                sampled_df = sampled_df.append(df, ignore_index=True)
        sampled_df.to_csv(self.output().path, ignore_index=True)


class VinaRandomizedBioLipFixed(VinaPredictBioLipFixed):
    def check(self):
        completes, incompletes = [], []
        for name in sampled_names:
            task = QueryVinaRandomResultOnBioLipFixedPocket(name)
            if task.complete():
                completes.append(name)
            else:
                incompletes.append(name)
        print("{} completes and {} incompletes".format(
            len(completes), len(incompletes)))
        return completes

    def output(self):
        path = "/work/jaydy/working/vina_biolip_sampled_rnd.csv"
        return luigi.LocalTarget(path)


class CuttedVinaPredictBioLipFixed(luigi.Task):
    def requires(self):
        return VinaPredictBioLipFixed()

    def output(self):
        path = os.path.splitext(self.requires().output().path)[
            0] + '.cutted.csv'
        return luigi.LocalTarget(path)

    def run(self):
        df = pd.read_csv(self.requires().output().path, index_col=0)
        cutted = df.groupby('query', as_index=False)\
                   .apply(cutRedundantTemplates)
        cutted.set_index('query', inplace=True)
        cutted.to_csv(self.output().path, ignore_index=True)


class UnCuttedVinaPredictBioLip(CuttedVinaPredictBioLip):
    def output(self):
        path = "/work/jaydy/working/uncutted_vina_biolip_sampled.csv"
        return luigi.LocalTarget(path)

    def run(self):
        sampled_df = pd.DataFrame()
        for name in self.check():
            task = vina_predict_biolip.QueryVinaResultOnBioLip(name)
            with task.output().open('r') as ifs:
                result = json.loads(ifs.read())
                df = read2Df(result, name)
                sampled_df = sampled_df.append(df, ignore_index=True)
        sampled_df.to_csv(self.output().path, ignore_index=True)


class CheckVinaResultAccuracy(luigi.Task):
    def output(self):
        path = "../dat/vina_biolip_actual_accuracy.csv"
        return luigi.LocalTarget(path)

    def check(self):
        completes, incompletes = [], []
        for name in sampled_names:
            task = vina_predict_biolip.VinaResultAccuracy(name)
            if task.complete():
                completes.append(name)
            else:
                incompletes.append(name)
        print("{} completes and {} incompletes".format(
            len(completes), len(incompletes)))

        return completes

    def run(self):
        index, data = [], []
        for name in self.check():
            task = vina_predict_biolip.VinaResultAccuracy(name)
            with task.output().open('r') as ifs:
                result = json.loads(ifs.read())
                data.append(result['rmsd'])
                index.append(name)

        df = pd.DataFrame({'query': index, 'rmsd': data})
        df.to_csv(self.output().path, ignore_index=True)


def clean(df):
    print("{} queries and {} records in the original dataset".format(df[
        'query'].unique().size, df.shape[0]))
    # regular ps_score
    df = df[(df.ps_score < 1.001) & (df.ps_score) > 0.00]
    print("{} queries and {} records after removing wield ps-score".format(df[
        'query'].unique().size, df.shape[0]))
    # different systems
    df = df[(df.seq_identity < 0.9) & (df.Tc < 0.9)]
    print("{} queries after filtering same systems".format(df['query'].unique(
    ).size))
    # drop nan
    df = df.dropna()
    print("{} queries after droping nan".format(df['query'].unique().size))

    return df


def similarPocketsLigands(df):
    queries_before = df['query'].unique().size
    df = df[(df.Tc > 0.5) & (df.ps_score > 0.4)]
    print("{} queries after filtering dissimilar pockets and ligands".format(
        df['query'].unique().size))
    queries_after = df['query'].unique().size
    print("%.3f%% survive" % (float(queries_after) / queries_before * 100))
    return df


def correlation(df):
    tc_times_ps = df.Tc * df.ps_score
    print("Maximal information-based stats")
    pprint(minepy.minestats(tc_times_ps, df['spearmanr']))


def ratio(array, filter_fn=lambda x: x > 0):
    filtered = [_ for _ in array if filter_fn(_)]
    r = float(len(filtered)) / len(array)
    return r


def preprocess(df):
    df = clean(df.copy())
    # correlation(df)
    df = similarPocketsLigands(df.copy())
    # pprint(df.describe())
    print("%.3f of the spearmanr > 0" % ratio(df.spearmanr))
    print("%.3f of the p-value < 0.05" % ratio(df[df.spearmanr > 0]['pval'],
                                               filter_fn=lambda x: x < 0.05))
    return df.copy()


def analysis():

    print("native structures in biolip querying biolip")
    df = pd.read_csv(Read().output().path, index_col=0)

    print
    print("predicted structures in biolip querying biolip")
    fixed_df = pd.read_csv(VinaPredictBioLipFixed().output().path, index_col=0)

    df_queries = set(df['query'].tolist())
    fixed_df_queries = set(fixed_df['query'].tolist())
    shared_queries = df_queries & fixed_df_queries

    df = df[df['query'].isin(shared_queries)]
    fixed_df = fixed_df[fixed_df['query'].isin(shared_queries)]

    # def checkSameTemplates(df1, df2):
    #     shared_queries = set(df1['query']) & set(df2['query'])
    #     same, different = [], []
    #     for query in shared_queries:
    #         tmp1 = df1[df1['query'] == query]["template"]
    #         tmp2 = df2[df2['query'] == query]["template"]
    #         if set(tmp1) == set(tmp2):
    #             same.append(query)
    #         else:
    #             different.append(query)
    #     return same, different

    # same, different = checkSameTemplates(fixed_df, df)

    # TODO: why same query found different templates, such as 2h3q_PBU_A_1.pdb
    # TODO: because of different settings for the binding site volumne

    df = preprocess(df.copy())

    print
    fixed_df = preprocess(fixed_df.copy())

    ifn = CheckVinaResultAccuracy().output().path
    actual_rmsd_df = pd.read_csv(ifn, index_col=0)

    fixed_accuracies = fixed_df[["query", "spearmanr", "pval", "tc_times_ps",
                                 "TM-score", "template"]]

    fixed_rmsd_accuracies = pd.merge(fixed_accuracies, actual_rmsd_df)

    # fixed_rmsd_accuracies.plot(kind='scatter', x='rmsd', y="spearmanr")

    # plt.savefig('fixed_rmsd_spearmanr.png')

    # print minepy.minestats(fixed_rmsd_accuracies.rmsd,
    #                        fixed_rmsd_accuracies.spearmanr)

    def comparedWithRmsd(fixed_rmsd_accuracies):
        def calculateAUC(scores):
            scores["native_like"] = scores.rmsd.apply(
                lambda r: 1 if r < 3 else 0)
            fpr, tpr, thresholds = metrics.roc_curve(scores.native_like,
                                                     scores.spearmanr,
                                                     pos_label=1)
            return metrics.auc(fpr, tpr)

        '''
        # averaged spearmanr
        print(
            "################################################################################")
        ave_spearmanr = fixed_rmsd_accuracies.groupby("query").apply(
            lambda g: g[["spearmanr", "rmsd"]].mean())
        print("AUC:", calculateAUC(ave_spearmanr))
        pprint(ave_spearmanr.corr())
        pprint(pd.Series(minepy.minestats(ave_spearmanr.rmsd,
                                          ave_spearmanr.spearmanr)))

        # ave_spearmanr.plot(kind='scatter', x='rmsd', y="spearmanr")

        # plt.savefig('rmsd_ave_spearmanr.png')

        # tc_times_ps weighted spearmanr
        def weighted(g):
            score = g["tc_times_ps"].dot(g["spearmanr"]) / g[
                "tc_times_ps"].sum()
            return pd.DataFrame(
                {"rmsd": g["rmsd"].tolist()[0],
                 "spearmanr": score},
                index=[0])

        # tc_times_ps weighted spearmanr
        print(
            "################################################################################")
        weighted_spearmanr = fixed_rmsd_accuracies.groupby("query").apply(
            weighted)
        print("AUC:", calculateAUC(weighted_spearmanr))
        pprint(weighted_spearmanr.corr())
        pprint(pd.Series(minepy.minestats(weighted_spearmanr.rmsd,
                                          weighted_spearmanr.spearmanr)))
        # weighted_spearmanr.plot(kind='scatter', x='rmsd', y="spearmanr")
        # plt.savefig('rmsd_weighted_spearmanr.png')
        '''

        # tc_times_ps ranked spearmanr
        print(
            "################################################################################"
        )
        ranked_spearmanr = fixed_rmsd_accuracies.groupby("query").apply(
            lambda g: g.sort_values("tc_times_ps", ascending=False).iloc[0][["spearmanr", "rmsd"]])
        print "AUC:", calculateAUC(ranked_spearmanr)
        print
        pprint(ranked_spearmanr.corr())
        print
        pprint(pd.Series(minepy.minestats(ranked_spearmanr.rmsd,
                                          ranked_spearmanr.spearmanr)))
        # ranked_spearmanr.plot(kind='scatter', x='rmsd', y="spearmanr")
        # plt.savefig('rmsd_ranked_spearmanr.png')

    comparedWithRmsd(fixed_rmsd_accuracies[fixed_rmsd_accuracies["TM-score"] >
                                           0.5])
    comparedWithRmsd(fixed_rmsd_accuracies[fixed_rmsd_accuracies["TM-score"] <
                                           0.5])


if __name__ == "__main__":
    luigi.build(
        [
            Read(),
            # CutRedundancy(),
            CuttedVinaPredictBioLip(),
            UnCuttedVinaPredictBioLip(),
            CheckVinaResultAccuracy(),
            VinaPredictBioLipFixed(),
            VinaRandomizedBioLipFixed(),
            # CuttedVinaPredictBioLipFixed(),
        ],
        local_scheduler=True)
