import pandas as pd

import json
import os
import random
from pprint import pprint

import luigi

import biolip_query_biolip
import vina_predict_biolip

from vina_predict_biolip import QueryVinaResultOnBioLipFixedPocket

sampled_list = "../dat/biolipbiolip_sampled_2.txt"
sampled_names = [_.rstrip() for _ in file(sampled_list)]

unfinished = []
finished = []
for name in sampled_names:
    task = biolip_query_biolip.BioLipBioLip(name)
    if not task.complete():
        unfinished.append(name)
    else:
        finished.append(name)

len(unfinished)

path = biolip_query_biolip.Path(unfinished[0]).prtPdb
if not os.path.exists(path):
    print '''%s fails because %s does not exist
    ''' % (unfinished[0], path)

missing_proteins_list = []
for name in unfinished:
    task = biolip_query_biolip.BioLipBioLip(name)
    if not os.path.exists(task.prtPdb):
        missing_proteins_list.append(name)

len(missing_proteins_list)

set(unfinished) - set(missing_proteins_list)

task = biolip_query_biolip.BioLipBioLip(finished[0])

result = []
with task.output().open('r') as ifs:
    result = json.loads(ifs.read())
    pprint(result[:2])


def read2Df(result, task_name):
    df = pd.DataFrame(dict(result)).T
    df['query'] = [task_name] * len(df)
    df['template'] = df.index
    return df


class Read(luigi.Task):
    def output(self):
        path = "../dat/biolip_sampled.csv"
        return luigi.LocalTarget(path)

    def run(self):
        sampled_finished = random.sample(finished, 2000)

        sampled_df = pd.DataFrame()
        for name in sampled_finished:
            task = biolip_query_biolip.BioLipBioLip(name)
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
        sampled_df = pd.read_csv(Read().output().path, index_col=0)
        cutted = cutRedundantTemplates(sampled_df)
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


class CuttedVinaPredictBioLipFixed(luigi.Task):
    def requires(self):
        return VinaPredictBioLipFixed()

    def output(self):
        path = os.path.splitext(self.requires().output().path)[
            0] + '.cutted.csv'
        return luigi.LocalTarget(path)

    def run(self):
        df = pd.read_csv(self.requires().output().path, index_col=0)
        cutted = cutRedundantTemplates(df)
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


if __name__ == "__main__":
    luigi.build(
        [
            Read(),
            CutRedundancy(),
            CuttedVinaPredictBioLip(),
            UnCuttedVinaPredictBioLip(),
            CheckVinaResultAccuracy(),
            VinaPredictBioLipFixed(),
            CuttedVinaPredictBioLipFixed(),
        ],
        local_scheduler=True)
