#!/usr/bin/env python
"""paper figures
"""
from __future__ import print_function

import re
import pandas as pd
import numpy as np
import analysis as readme

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from pprint import pprint
from analysis import analysis, loadData, clean

sns.set_style("white", {'ytick.major.size': 10.0})
sns.set_context("poster", font_scale=1.2)
matplotlib.rcParams.update({'font.size': 24, 'font.family': 'serif'})


def main():
    native_df, fixed_df, rnd_df = loadData()

    native_df.head()

    native_df = readme.clean(native_df)
    fixed_df = readme.clean(fixed_df)
    rnd_df = readme.clean(rnd_df)

    fixed_df[(fixed_df.num_binding_res < 30) & (
        fixed_df.num_lig_atoms < 20)].sort_values(
            ['num_binding_res', 'num_lig_atoms'],
            ascending=False).head()

    def normalizeRanges(s):
        lo, hi = re.findall(r'[+-]?\d+\.*\d*', s)
        return "(%.3f, %.3f]" % (float(lo), float(hi))

    plt.figure()
    df = native_df
    df = df[(df.pval > 0) & (df.spearmanr > 0)]

    data = df.groupby([pd.cut(df['Tc'], 20), pd.cut(df['ps_score'], 20)])[
        'pval'].apply(lambda pvals: np.mean(np.log(pvals))).unstack()
    data = data.reindex(index=data.index[::-1])
    sns.heatmap(data,
                square=True,
                xticklabels=map(normalizeRanges, data.index[::-1]),
                yticklabels=map(normalizeRanges, data.columns[::-1]))
    plt.xlabel("Pocket similarity", fontsize=24)
    plt.ylabel("Ligand similarity", fontsize=24)
    plt.tight_layout()
    plt.savefig("/work/jaydy/working/xcms_plot/p_val.tiff", dpi=200)

    plt.figure()
    df = native_df
    data = df.groupby([pd.cut(df['Tc'], 20), pd.cut(df['ps_score'], 20)])[
        'spearmanr'].mean().unstack()
    data = data.reindex(index=data.index[::-1])
    sns.heatmap(data,
                square=True,
                xticklabels=map(normalizeRanges, data.index[::-1]),
                yticklabels=map(normalizeRanges, data.columns[::-1]),
                vmin=-1,
                vmax=1)
    plt.xlabel("Pocket similarity", fontsize=24)
    plt.ylabel("Ligand similarity", fontsize=24)
    plt.tight_layout()
    plt.savefig("/work/jaydy/working/xcms_plot/native_xcms.tiff", dpi=200)

    plt.figure()
    df = fixed_df
    data = df.groupby([pd.cut(df['Tc'], 20), pd.cut(df['ps_score'], 20)])[
        'spearmanr'].mean().unstack()
    data = data.reindex(index=data.index[::-1])
    ax = sns.heatmap(data,
                     square=True,
                     xticklabels=map(normalizeRanges, data.index[::-1]),
                     yticklabels=map(normalizeRanges, data.columns[::-1]),
                     vmin=-1,
                     vmax=1)
    plt.xlabel("Pocket similarity")
    plt.ylabel("Ligand similarity")
    plt.tight_layout()
    plt.savefig("/work/jaydy/working/xcms_plot/predicted_xcms.tiff", dpi=200)

    plt.figure()
    df = rnd_df
    data = df.groupby([pd.cut(df['Tc'], 20), pd.cut(df['ps_score'], 20)])[
        'spearmanr'].mean().unstack()
    data = data.reindex(index=data.index[::-1])
    ax = sns.heatmap(data,
                     square=True,
                     xticklabels=map(normalizeRanges, data.index[::-1]),
                     yticklabels=map(normalizeRanges, data.columns[::-1]),
                     vmin=-1,
                     vmax=1)
    plt.xlabel("Pocket similarity")
    plt.ylabel("Ligand similarity")
    plt.tight_layout()
    plt.savefig("/work/jaydy/working/xcms_plot/random_xcms.tiff", dpi=200)

    predicted_rmsd = pd.read_csv(
        readme.CheckVinaResultAccuracy().output().path,
        index_col=0)

    plt.figure()
    predicted_rmsd.rmsd.hist(bins=30, color='darkgrey', alpha=0.9)
    plt.xlabel("RMSD [$\AA$]")
    plt.savefig("/work/jaydy/working/xcms_plot/predict_rmsd_hist.tiff",
                dpi=200)

    plt.figure()
    predicted_rmsd.cms.hist(bins=30, color='darkgrey', alpha=0.9)
    plt.xlabel("CMS")
    plt.xlim((0, 1))
    plt.savefig("/work/jaydy/working/xcms_plot/predict_cms_hist.tiff", dpi=200)

    rnd_rmsd = pd.read_csv(readme.CheckVinaRandomRmsd().output().path,
                           index_col=0)

    plt.figure()
    rnd_rmsd.rmsd.hist(bins=30, color='darkgrey', alpha=0.9)
    plt.xlabel("RMSD [$\AA$]")
    plt.savefig("/work/jaydy/working/xcms_plot/random_rmsd_hist.tiff", dpi=200)

    plt.figure()
    rnd_rmsd.cms.hist(bins=30, color='darkgrey', alpha=0.9)
    plt.xlabel("CMS")
    plt.xlim((0, 1))
    plt.savefig("/work/jaydy/working/xcms_plot/random_cms_hist.tiff", dpi=200)

    native_df = readme.similarPocketsLigands(native_df)
    fixed_df = readme.similarPocketsLigands(fixed_df)
    rnd_df = readme.similarPocketsLigands(rnd_df)

    fixed_df.groupby("query").apply(lambda g: g.shape[0]).hist(
        bins=30,
        color='darkgrey',
        alpha=0.9)

    fiexed_xcms = fixed_df[["query", "spearmanr", "pval", "tc_times_ps",
                            "TM-score", "template"]]

    fixed_rmsd_xcms = pd.merge(fiexed_xcms, predicted_rmsd)

    fixed_spearmanr = fixed_rmsd_xcms.groupby("query").apply(
        lambda g: g.sort_values("tc_times_ps", ascending=False).iloc[0][["spearmanr", "rmsd", "TM-score", "cms"]])

    plt.figure()
    fixed_spearmanr.spearmanr.hist(bins=30, color='darkgrey', alpha=0.9)
    plt.xlabel("XCMS")
    plt.savefig("/work/jaydy/working/xcms_plot/predict_xcms_hist.tiff",
                dpi=200)

    rnd_xcms = rnd_df[["query", "spearmanr", "pval", "tc_times_ps", "TM-score",
                       "template"]]

    rnd_rmsd_xcms = pd.merge(rnd_xcms, rnd_rmsd)

    rnd_spearmanr = rnd_rmsd_xcms.groupby("query").apply(
        lambda g: g.sort_values("tc_times_ps", ascending=False).iloc[0][["spearmanr", "rmsd", "cms"]])

    plt.figure()
    rnd_spearmanr.spearmanr.hist(bins=30, color='darkgrey', alpha=0.9)
    plt.xlabel("XCMS")
    plt.savefig("/work/jaydy/working/xcms_plot/random_xcms_hist.tiff", dpi=200)

    fig = plt.figure(figsize=(7, 6))
    # fixed_spearmanr.plot(kind="scatter",
    #                      x='cms',
    #                      y='rmsd',
    #                      color='k',
    #                      alpha=0.3)
    ax = fig.add_subplot(1,1,1)
    ax.scatter(fixed_spearmanr['cms'],
                fixed_spearmanr['rmsd'],
                alpha=0.3,
                c='k')

    # z = np.polyfit(fixed_spearmanr.cms, fixed_spearmanr.rmsd, 3)
    # p = np.poly1d(z)
    # x = np.linspace(0, 1, num=50, endpoint=True)
    # fitted_line = p(x)
    # plt.plot(x, fitted_line, 'r--')

    # z = np.polyfit(fixed_spearmanr.cms, fixed_spearmanr.rmsd, 2)
    # p = np.poly1d(z)
    # x = np.linspace(0, 1, num=50, endpoint=True)
    # fitted_line = p(x)
    # plt.plot(x, fitted_line, 'b--')

    # z = np.polyfit(fixed_spearmanr.cms, fixed_spearmanr.rmsd, 1)
    # p = np.poly1d(z)
    # x = np.linspace(0, 1, num=50, endpoint=True)
    # fitted_line = p(x)
    # plt.plot(x, fitted_line, 'r--')
    # fixed_spearmanr[['cms', 'rmsd']].to_csv("/work/jaydy/working/cms_rmsd_scatter.csv")

    ax.set_xlabel("CMS", fontsize=24)
    ax.set_ylabel("RMSD [$\AA$]", fontsize=24)
    ax.set_xlim((0, 1))
    ax.set_ylim((0, 15))
    fig.tight_layout()
    fig.savefig("/work/jaydy/working/xcms_plot/cms_rmsd_scatter.tiff", dpi=200)

    high_tm = fixed_spearmanr[fixed_spearmanr["TM-score"] > 0.5]
    low_tm = fixed_spearmanr[fixed_spearmanr["TM-score"] < 0.5]

    high_tm.to_csv("/work/jaydy/working/xcms_plot/high_tm.csv",
                   ignore_index=True)
    low_tm.to_csv("/work/jaydy/working/xcms_plot/low_tm.csv",
                  ignore_index=True)

    plt.figure()
    # high_tm.plot(kind="scatter", x="cms", y="spearmanr", color='g', alpha=0.5, label='High global similarity')
    plt.scatter(high_tm.cms,
                high_tm.spearmanr,
                color='g',
                alpha=0.5,
                label='High global similarity')
    # plt.xlabel("CMS")
    # plt.ylabel("XCMS")
    # plt.savefig("/work/jaydy/working/xcms_plot/cms_xcms_high_scatter.tiff", dpi=200)

    # plt.figure()
    # low_tm.plot(kind="scatter", x="cms", y="spearmanr", color='r', alpha=0.5, label='Low global similarity')
    plt.scatter(low_tm.cms,
                low_tm.spearmanr,
                color='r',
                alpha=0.5,
                label='Low global similarity')

    plt.xlim((0, 1))
    plt.ylim((-0.6, 1))
    # plt.xlim((0, 1))
    # plt.ylim((-0.6, 1))
    # plt.legend(loc='upper left')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
               loc=3,
               ncol=2,
               mode="expand",
               borderaxespad=0.)
    plt.xlabel("CMS")
    plt.ylabel("XCMS")
    plt.savefig("/work/jaydy/working/xcms_plot/cms_xcms_high_low_scatter.tiff",
                dpi=200)

    # plt.figure()
    # high_tm.plot(kind="scatter", x="rmsd", y="spearmanr", color='k', alpha=0.5)
    # plt.xlim((0, 18))
    # plt.ylim((-0.8, 1))
    # plt.xlabel("RMSD [$\AA$]")
    # plt.ylabel("XCMS")
    # plt.savefig("/work/jaydy/working/xcms_plot/rmsd_xcms_high_scatter.tiff", dpi=200)

    # plt.figure()
    # low_tm.plot(kind="scatter", x="rmsd", y="spearmanr", color='k', alpha=0.5)
    # plt.xlim((0, 18))
    # plt.ylim((-0.8, 1))
    # plt.xlabel("RMSD [$\AA$]")
    # plt.ylabel("XCMS")
    # plt.savefig("/work/jaydy/working/xcms_plot/rmsd_xcms_low_scatter.tiff", dpi=200)

    plt.figure()
    from sklearn import metrics
    high_tm['native_like'] = high_tm.rmsd.apply(lambda r: 1 if r < 3 else 0)
    low_tm['native_like'] = low_tm.rmsd.apply(lambda r: 1 if r < 3 else 0)
    hfpr, htpr, _ = metrics.roc_curve(high_tm.native_like,
                                      high_tm.spearmanr,
                                      pos_label=1)
    lfpr, ltpr, _ = metrics.roc_curve(low_tm.native_like,
                                      low_tm.spearmanr,
                                      pos_label=1)

    plt.plot(hfpr, htpr, 'k-', label="high global similarity")
    plt.plot(lfpr, ltpr, 'k--', label="low global similarity")
    plt.plot([0, 1], [0, 1], 'k:', label="random guess")
    plt.legend(loc='lower right')

    print("High global similarity AUC: %f" % metrics.auc(hfpr, htpr))
    print("Low global similarity AUC: %f" % metrics.auc(lfpr, ltpr))
    plt.savefig("/work/jaydy/working/xcms_plot/xcms_roc.tiff", dpi=200)

    plt.figure()
    rnd_spearmanr.plot(kind="scatter",
                       x="rmsd",
                       y="spearmanr",
                       color='k',
                       alpha=0.5)
    plt.xlim((0, 30))
    plt.ylim((-1, 1))
    plt.xlabel("RMSD [$\AA$]")
    plt.ylabel("XCMS")

    xcms_back_df = pd.read_csv(
        readme.CheckVinaResultAgainstIdenticalSystems().output().path,
        index_col=0)

    xcms_back_df['query'] = xcms_back_df.index.values

    predicted_df = pd.merge(predicted_rmsd,
                            xcms_back_df[['query', 'spearmanr']])

    fig = plt.figure(figsize=(7, 6))

    sampled_predicted_df = predicted_df.sample(len(fixed_spearmanr))
    # sampled_predicted_df.plot(kind='scatter',
    #                           x='cms',
    #                           y='spearmanr',
    #                           color='k',
    #                           alpha=0.3)
    ax = fig.add_subplot(1,1,1)
    ax.scatter(sampled_predicted_df['cms'],
                sampled_predicted_df['spearmanr'],
                c='k',
                alpha=0.3)

    # cleaned_predicted_df = predicted_df.dropna()

    # z = np.polyfit(cleaned_predicted_df.cms, cleaned_predicted_df.spearmanr, 3)
    # p = np.poly1d(z)
    # x = np.linspace(0, 1, num=50, endpoint=True)
    # fitted_line = p(x)
    # plt.plot(x, fitted_line, 'r--')

    # z = np.polyfit(cleaned_predicted_df.cms, cleaned_predicted_df.spearmanr, 2)
    # p = np.poly1d(z)
    # x = np.linspace(0, 1, num=50, endpoint=True)
    # fitted_line = p(x)
    # plt.plot(x, fitted_line, 'b--')

    # z = np.polyfit(cleaned_predicted_df.cms, cleaned_predicted_df.spearmanr, 1)
    # p = np.poly1d(z)
    # x = np.linspace(0, 1, num=50, endpoint=True)
    # fitted_line = p(x)
    # plt.plot(x, fitted_line, 'r--')

    # cleaned_predicted_df[['cms', 'spearmanr']].to_csv("/work/jaydy/working/cms_xcms_scatter.csv")

    ax.set_xlabel('CMS', fontsize=24)
    ax.set_ylabel("XCMS", fontsize=24)
    ax.set_xlim((0, 1))
    ax.set_ylim((-1, 1))
    fig.tight_layout()
    fig.savefig("/work/jaydy/working/xcms_plot/cms_xcms_scatter.tiff", dpi=200)

    plt.figure()
    plt.hist([rnd_spearmanr.spearmanr, fixed_spearmanr.spearmanr],
             30,
             color=['lightgrey', 'dimgrey'],
             label=['Random', 'AutoDock Vina'])
    plt.xlabel('XCMS', fontsize=24)
    plt.ylabel('Number of complexes', fontsize=24)
    # plt.xlim((-1, 1))
    plt.legend(loc='best')
    plt.savefig("/work/jaydy/working/xcms_plot/rnd_pred_hist_xcms.tiff",
                dpi=200)

    print(rnd_spearmanr.spearmanr.describe())
    print(fixed_spearmanr.spearmanr.describe())

    plt.figure()
    plt.hist([rnd_rmsd.rmsd, predicted_rmsd.rmsd],
             30,
             color=['lightgrey', 'dimgrey'],
             label=['Random', 'AutoDock Vina'])
    plt.xlabel('RMSD [$\AA$]', fontsize=24)
    plt.ylabel('Number of complexes', fontsize=24)
    plt.legend(loc='best')
    plt.savefig("/work/jaydy/working/xcms_plot/rnd_pred_hist_rmsd.tiff",
                dpi=200)
    print(rnd_rmsd.rmsd.describe())
    print(predicted_rmsd.rmsd.describe())

    plt.figure()
    plt.hist([rnd_rmsd.cms, predicted_rmsd.cms],
             30,
             color=['lightgrey', 'dimgrey'],
             label=['Random', 'AutoDock Vina'])
    plt.xlabel('CMS', fontsize=24)
    plt.ylabel('Number of complexes', fontsize=24)
    plt.xlim((0, 1))
    plt.legend(loc='best')
    plt.savefig("/work/jaydy/working/xcms_plot/rnd_pred_hist_cms.tiff",
                dpi=200)
    print(rnd_rmsd.cms.describe())
    print(predicted_rmsd.cms.describe())


if __name__ == '__main__':
    main()
