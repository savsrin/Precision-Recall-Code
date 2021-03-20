#!/usr/bin/env python3

"""

scatterplotting intron lengths between mouse and human

"""
import argparse
from datetime import datetime
import gzip
from matplotlib import use; use('pdf')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import pandas as pd
from scipy import stats
import seaborn as sns; sns.set(color_codes=True)


def scatter_with_size(kdf_init, dbdf_init, cancer, can_count, out_path,
                      no_zeros=False):
    """

    :param kdf_init:
    :param dbdf_init:
    :param cancer:
    :param can_count:
    :param out_path:
    :param no_zeros:
    :return:
    """
    merge_df = pd.merge(kdf_init, dbdf_init, on=['jx', 'cancer'], how='outer')
    merge_df.drop(['cancer'], axis=1, inplace=True)
    if no_zeros:
        merge_df.dropna(inplace=True)
    else:
        mutual_jxs = len(merge_df.dropna())
        merge_df.fillna(0, inplace=True)

    num_jxs = len(merge_df)

    groups = merge_df.groupby(['rail_prev', 'k_prev'])
    scatter_dict = {'rail_prev': [], 'k_prev': [], 'count': []}
    for set, group in groups:
        scatter_dict['rail_prev'].append(set[0])
        scatter_dict['k_prev'].append(set[1])
        scatter_dict['count'].append(group.jx.count())

    scatter_df = pd.DataFrame(scatter_dict)

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 5.0, 5.0
    sns.set_context("paper")
    sns.set_style("whitegrid")

    # ax = sns.scatterplot(
    #     x='k_prev', y='rail_prev', size='count', data=scatter_df,
    #     facecolor='xkcd:bright blue', edgecolor='None', s=4
    # )
    # ax = sns.scatterplot(
    #     x='k_prev', y='rail_prev', size='count', data=scatter_df,
    #     facecolor='xkcd:yellow orange', edgecolor='xkcd:pumpkin'
    # )
    ax = sns.scatterplot(
        x='k_prev', y='rail_prev', size='count', data=scatter_df,
        facecolor='xkcd:light grey blue', edgecolor='xkcd:bluish'
    )
    # fig_name = '{}_rail_vs_k_neojx_prev.pdf'.format(cancer)
    if no_zeros:
        fig_name = '{}_rail_vs_k_neojx_prev_no_zeros.pdf'.format(cancer)
        title = (
            '{}\n{} TCGA samples, {} unique mutual junctions'
            ''.format(cancer, can_count, num_jxs)
        )
    else:
        fig_name = '{}_rail_vs_k_neojx_prev.pdf'.format(cancer)
        title = (
            '{}\n{} TCGA samples, {} total unique jxs, {} mutual jxs'
            ''.format(cancer, can_count, num_jxs, mutual_jxs)
        )

    plt.title(title)
    plt.xlabel('Kahles et al neojx prevalence')
    # plt.xticks(x_locs, x_labels)
    # plt.setp(ax.get_xticklabels(), rotation=90, fontsize=6)
    # plt.yscale('log')
    # plt.yscale('symlog', linthresy=2500)
    # ax.set_ylim(ymin=0)
    ax.set_ylim(ymin=0, ymax=1.0)
    ax.set_xlim(xmin=0, xmax=1.0)
    # ax.xaxis.grid(False)
    # ax.tick_params(axis='both', which='major', direction='out', length=3, width=2)
    plt.ylabel('rail-aligned neojx prevalence')
    fig = plt.gcf()
    fig_file = os.path.join(out_path, fig_name)
    print(fig_file)
    fig.savefig(fig_file)
    plt.clf()
    return


def joint_plot(kdf_init, dbdf_init, cancer, can_count, out_path, kind='reg',
               no_zeros=True):
    """

    :param kdf_init:
    :param dbdf_init:
    :param cancer:
    :param can_count:
    :param out_path:
    :param kind:
    :param no_zeros:
    :return:
    """
    scatter_df = pd.merge(
        kdf_init, dbdf_init, on=['jx', 'cancer'], how='outer'
    )
    scatter_df.drop(['cancer'], axis=1, inplace=True)
    if no_zeros:
        scatter_df.dropna(inplace=True)
        mutual_jxs = len(scatter_df)
    else:
        mutual_jxs = len(scatter_df.dropna())
        scatter_df.fillna(0, inplace=True)

    num_jxs = len(scatter_df)

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 5.0, 5.0
    sns.set_context("paper")
    sns.set_style("whitegrid")

    # ax = sns.scatterplot(
    #     x='k_prev', y='rail_prev', size='count', data=scatter_df,
    #     facecolor='xkcd:bright blue', edgecolor='None', s=4
    # )
    # ax = sns.scatterplot(
    #     x='k_prev', y='rail_prev', size='count', data=scatter_df,
    #     facecolor='xkcd:yellow orange', edgecolor='xkcd:pumpkin'
    # )
    ax = sns.jointplot(
        x='k_prev', y='rail_prev', data=scatter_df, kind=kind,
        color='xkcd:emerald green', xlim=(-0.005, 1.0), ylim=(-0.005, 1.0)
    )
    # ax = sns.jointplot(
    #     x='k_prev', y='rail_prev', data=scatter_df, kind=kind,
    #     color='xkcd:greyblue', xlim=(0, 1.0), ylim=(0, 1.0)
    # )
    if no_zeros:
        fig_name = (
            '{}_rail_vs_k_neojx_prev_joint_{}_no_zeros.pdf'
            ''.format(cancer, kind)
        )
        title = (
            '{}\n{} TCGA samples, {} unique mutual junctions'
            ''.format(cancer, can_count, num_jxs)
        )
    else:
        fig_name = '{}_rail_vs_k_neojx_prev_joint_{}.pdf'.format(cancer, kind)
        title = (
            '{}\n{} TCGA samples, {} total unique jxs, {} mutual jxs'
            ''.format(cancer, can_count, num_jxs, mutual_jxs)
        )
    try:
        abbr = _TCGA_ABBR[cancer]
    except KeyError:
        abbr = 'TCGA'
    title = '{}: {} samples, {} mutual jxs'.format(abbr, can_count, mutual_jxs)
    plt.title(title)
    plt.xlabel('Kahles et al neojx prevalence')
    # plt.xticks(x_locs, x_labels)
    # plt.setp(ax.get_xticklabels(), rotation=90, fontsize=6)
    # plt.yscale('log')
    # plt.yscale('symlog', linthresy=2500)
    # ax.set_ylim(ymin=0)
    # ax.set_ylim(ymin=0, ymax=1.0)
    # ax.set_xlim(xmin=0, xmax=1.0)
    # ax.xaxis.grid(False)
    # ax.tick_params(axis='both', which='major', direction='out', length=3, width=2)
    plt.ylabel('rail-aligned neojx prevalence')
    # plt.ylabel(
    #     '{} rail-aligned neojx prevalence\n{} samples, {} mutual jxs'
    #     ''.format(abbr, num_jxs, mutual_jxs))
    fig = plt.gcf()
    fig_file = os.path.join(out_path, fig_name)
    print(fig_file)
    fig.savefig(fig_file)
    plt.clf()
    return


def scatter_with_regression(x_vals, y_vals, flag, out_path, timestamp):
    """

    :param kdf_init:
    :param dbdf_init:
    :param cancer:
    :param can_count:
    :param out_path:
    :return:
    """
    x = x_vals
    y = y_vals
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line = [slope * x + intercept for x in x_vals]
    print('linear regression stats:')
    print(slope, intercept, r_value, p_value, std_err)
    if intercept > 0:
        line_label = (
            'y = {}x + {}\nR = {}'
            ''.format(round(slope, 3), round(intercept, 3), round(r_value, 2))
        )
    else:
        line_label = (
            'y = {}x - {}\nR = {}'
            ''.format(round(slope, 3), abs(round(intercept, 3)), round(r_value, 2))
        )
    ax_max = max(max(x_vals), max(y_vals)) * 1.02
    # ax_min = -0.1 * ax_max
    ax_min = -0.02 * ax_max
    # point_dict = {'fake_x': [], 'fake_y': [], 'color': []}
    # for color in scatter_dict['color']:
    #     point_dict['fake_x'].append(-ax_max * 10)
    #     point_dict['fake_y'].append(-ax_max * 10)
    #     point_dict['color'].append(color)
    # fake_df = pd.DataFrame(point_dict)

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 5.0, 5.0
    mpl.style.use('seaborn-whitegrid')
    edgecolor = '#047495'
    facecolor = '#d0fefe'
    edgecolor = '#040273'
    facecolor = '#0485d1'

    plt.scatter(
        x=x_vals, y=y_vals, c=facecolor, edgecolors=edgecolor, s=10,
        linewidths=0.05, label=None
    )
    # plt.scatter(
    #     x=x_vals, y=y_vals, linewidths=0.05, label=None
    # )
    # plt.scatter(
    #     x='fake_x', y='fake_y', c='color', data=fake_df, s=40,
    #     edgecolors='xkcd:light grey', linewidth=0.01, label='shared jx'
    # )
    plt.plot(
        x, line, c='xkcd:dark grey', linewidth=0.8, label=line_label
    )
    # plt.plot(
    #     x, line, c='xkcd:grey', linewidth=1, label=line_label
    # )
    # plt.legend(fontsize='small', markerscale=.4, frameon=True);
    plt.legend(fontsize='x-small', frameon=True)
    # lgnd = plt.legend(fontsize='small', scatterpoints=1, frameon=True)
    # lgnd.legendHandles[0]._legmarker.set_markersize(6)

    # lgnd = plt.legend(loc="lower left", scatterpoints=1, fontsize=10)
    # lgnd.legendHandles[0]._sizes = [30]
    # # lgnd.legendHandles[1]._sizes = [30]

    fig_name = 'm_v_h_intronlength_scatter_{}_{}.pdf'.format(flag, timestamp)
    plt.xlabel('human intron length', fontsize=8)
    plt.ylabel('mouse intron length', fontsize=8)
    # plt.xticks(x_locs, x_labels)

    # plt.setp(ax.get_xticklabels(), rotation=90, fontsize=6)
    # plt.yscale('log')
    # plt.yscale('symlog', linthresy=2500)
    ax = plt.gca()
    # ax.set_ylim(ymin=0)
    # ax.set_ylim(ymin=ax_min, ymax=ax_max)
    # ax.set_xlim(xmin=ax_min, xmax=ax_max)
    ax.set_xlim(xmin=0, xmax = 100000)
    ax.set_ylim(ymin=0, ymax = 100000)
    # plt.axis('equal')
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(6)
        tick.label.set_rotation('horizontal')
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)
        tick.label.set_rotation('vertical')
    # # ax.set_facecolor('w')
    # ax.set_aspect('equal', 'box')

    # ax.spines['bottom'].set_color('k')
    # ax.spines['top'].set_color('k')
    # ax.xaxis.label.set_color('k')
    # ax.tick_params(axis='x', colors='k')

    # ax.xaxis.grid(False)
    # ax.tick_params(axis='both', which='major', direction='out', length=3, width=2)
    fig = plt.gcf()
    fig_file = os.path.join(out_path, fig_name)
    fig.savefig(fig_file)
    plt.clf()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Check TCGA junctions for developmental tissue evidence.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='give path for output files: sampled spots, aligned junction '
             'reads, and SRA numbers with their p-values.'
    )
    parser.add_argument(
        '--junction-file', '-j',
        help='Merged h2m and m2h liftover junctions.'
    )

    args = parser.parse_args()
    out_path = args.output_path
    # log_mode = args.log_level
    jx_file = args.junction_file

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    human_not_shared = []
    human_shared = []
    mouse_not_shared = []
    mouse_shared = []

    with gzip.open(jx_file) as jx_info:
        for line in jx_info:
            items = line.decode('utf-8').strip().split('\t')
            h_len = int(items[4])
            m_len = int(items[21])
            if items[0] == '*' or items[-1] == '*':
                human_not_shared.append(h_len)
                mouse_not_shared.append(m_len)
            else:
                human_shared.append(h_len)
                mouse_shared.append(m_len)

    # print(
    #     'number of human_only and mouse only: {}, {}'
    #     ''.format(len(human_not_shared), len(mouse_not_shared))
    # )
    # print(
    #     'number of human shared and mouse shared: {}, {}'
    #     ''.format(len(human_shared), len(mouse_shared))
    # )
    # print(
    #     'number of total human and mouse: {}, {}'
    #     ''.format(len(human_shared + human_not_shared),
    #               len(mouse_shared + mouse_not_shared))
    # )
    # exit()
    # scatter_with_regression(
    #     human_shared, mouse_shared, 'jx_in_h&m', out_path, now
    # )
    num_longs = 0
    num_uneven = 0
    human_longer = 0
    # for i, (h, m) in enumerate(zip(human_not_shared, mouse_not_shared)):
    #     # if i == 20:
    #     #     exit()
    #     # print(h, m)
    #     if h > 100000 or m > 100000:
    #         num_longs += 1
    #     if (h - m > 90000) or (m - h > 90000):
    #         print(h, m)
    #         num_uneven += 1
    #     if h > m:
    #         human_longer += 1
    #
    # print('num longs = {}'.format(num_longs))
    # print('num unevens = {}'.format(num_uneven))
    # print(
    #     'humans longer in {}% of cases'.format(
    #         100 * human_longer / (len(human_shared + human_not_shared))
    #     )
    # )
    human_longer = 0
    mouse_longer = 0
    equal = 0
    mostly_equal = 0
    h_slightly_longer = 0
    m_slightly_longer = 0
    for i, (h, m) in enumerate(zip(human_shared, mouse_shared)):
        # # if i == 20:
        # #     exit()
        # # print(h, m)
        # if h > 100000 or m > 100000:
        #     num_longs += 1
        # if (h - m > 90000) or (m - h > 90000):
        #     print(h, m)
        #     num_uneven += 1
        if h == m:
            equal += 1
        elif abs(h-m) < 10:
            mostly_equal += 1
            if h-m > 0:
                h_slightly_longer += 1
            else:
                m_slightly_longer += 1
        elif h > m:
            human_longer += 1
        else:
            mouse_longer += 1

    print('num longs = {}'.format(num_longs))
    print('num unevens = {}'.format(num_uneven))
    print(
        'humans longer in {}% of cases'.format(
            100 * human_longer / (len(human_shared))
        )
    )
    print(
        'mouse longer in {}% of cases'.format(
            100 * mouse_longer / (len(human_shared))
        )
    )
    print(
        'neither longer in {}% of cases'.format(
            100 * equal / (len(human_shared))
        )
    )
    print(
        'mostly equal in {}% of cases'.format(
            100 * mostly_equal / (len(human_shared))
        )
    )
    print(
        'h slightly longer in {}% of cases'.format(
            100 * h_slightly_longer / (len(human_shared))
        )
    )
    print(
        'm slightly longer in {}% of cases'.format(
            100 * m_slightly_longer / (len(human_shared))
        )
    )
    scatter_with_regression(
        human_not_shared, mouse_not_shared, 'non_shared_jx', out_path, now
    )
    scatter_with_regression(
        human_shared + human_not_shared,
        mouse_shared + mouse_not_shared,
        'all_jx', out_path, now
    )