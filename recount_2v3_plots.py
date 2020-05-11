#!/usr/bin/env python3

"""
recount_2v3_jx_comp.py
Python 3 code for comparing RNA-seq junctions called across mutual samples by
the recount2 vs. recount3 Sequence Read Archive, TCGA, and GTEx analyses.

"""

import argparse
from datetime import datetime
import json
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import seaborn as sns

sns.set(color_codes=True)

_LOG_MODE = 'INFO'
_CANON_MOTIFS = ['GTAG', 'GCAG', 'ATAC']
_OTHER_MOT = 'other'
_V2_JX = 'v2_only'
_V3_JX = 'v3_only'
_MU_JX = 'mutual'
_COMPLEMENT = str.maketrans("ATCG", "TAGC")
_SRA = 'SRA'
_TCGA = 'TCGA'
_GTEX = 'GTEx'
_JXS = 'jxs'
_EVENTS = 'events'


def grouped_boxplots(data_dict, plot_dict, fig_file, fig_size=(7.0, 4.0),
                     logscale=False, y_label='recall', percent=True,
                     right_lim_shift=2, x_label='splice site motif'):
    """

    :param data_dict:
    :param plot_dict:
    :param fig_file:
    :param fig_size:
    :param logscale:
    :param y_label:
    :param percent:
    :param right_lim_shift:
    :param x_label:
    :return:
    """
    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = fig_size
    sns.set_context("paper")
    sns.set_style("whitegrid")
    sns.set(style="ticks")

    f, ax = plt.subplots(nrows=1, ncols=1)
    plt.sca(ax)
    light_cols = plot_dict['light colors']
    dark_cols = plot_dict['dark colors']
    num_boxes = len(light_cols)
    adjust_val = num_boxes + 1

    bp_pos = list(range(1, num_boxes + 1))

    label_locs = []
    labels = [str(key) for key in data_dict.keys()]
    curr_label_loc = sum(bp_pos) / len(bp_pos)

    for abbr, values in data_dict.items():
        labels.append(abbr)
        data = values['data']

        boxes = plt.boxplot(
            data, positions=bp_pos, widths=0.6, patch_artist=True
        )
        box_elements = zip(boxes['fliers'], boxes['medians'], boxes['boxes'])
        for i, (fly, med, box) in enumerate(box_elements):
            plt.setp(
                fly, color=light_cols[i], markerfacecolor=light_cols[i],
                marker='.', markersize=4.0, linestyle='none', linewidth=0.15,
                markeredgecolor=dark_cols[i]
            )
            plt.setp(
                med, color=dark_cols[i], linewidth=1.5
            )
            plt.setp(box, facecolor='white')

        label_locs.append(curr_label_loc)

        bp_pos = [x + adjust_val for x in bp_pos]
        curr_label_loc += adjust_val

    ax.set_xticklabels(labels)
    ax.set_xticks(label_locs)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)
    ax.set_xlim(left=0, right=curr_label_loc - right_lim_shift)
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')

    if logscale:
        plt.yscale('log')

    plt.ylabel(y_label, fontsize=7)
    plt.xlabel(x_label, fontsize=7)
    plt.setp(
        ax.xaxis.get_majorticklabels(), rotation=0, fontsize=7, color='black'
    )
    # plt.xticks(label_locs, labels)

    plt.setp(
        ax.yaxis.get_majorticklabels(), fontsize=7, color='black'
    )

    legend_list = []
    for p_col, p_lab in zip(plot_dict['dark colors'], plot_dict['row labels']):
        legend_list.append(mpatches.Patch(color=p_col, label=p_lab))

    plt.legend(handles=legend_list, prop={'size': 6})

    if percent:
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(
                lambda y, _: '{}%'.format('{:,g}'.format(100 * y))
            )
        )
    else:
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda y, _: '{:,g}'.format(y))
        )
    # plt.xticks([])
    fig = plt.gcf()
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


def prepare_boxplot_data(prepped_data_json, outpath, now, cohort):

    with open(prepped_data_json) as recover:
        sample_jx_dict = json.load(recover)

    plot_info_dict = {
        'light colors': ['#29A1E8', '#85C9F2', '#E8321F', '#F2867B'],
        'dark colors': ['#29A1E8', '#85C9F2', '#E8321F', '#F2867B'],
        'row labels': [
            'v3 junction recall (v2 ground truth)',
            'v3 detection event recall (v2 ground truth)',
            'v2 junction recall (v3 ground truth)',
            'v2 detection event recall (v3 ground truth)'
        ]
    }
    v3skips = 0
    v2skips = 0
    v2skip_dict = {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]}
    v3skip_dict = {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]}
    v3nonskips = 0
    v2nonskips = 0
    v3_gtag_blank = []
    # Begin processing of all jxs/all motifs
    plot_data_dict = {mot: {'data': [[], [], [], []]} for mot in _CANON_MOTIFS}
    plot_data_dict['all motifs'] = {'data': [[], [], [], []]}
    # Process each canonical motif
    for sample, data in sample_jx_dict.items():
        allmot_mut_jx = 0
        allmot_mut_ev = 0
        allmot_v2_jx = 0
        allmot_v2_ev = 0
        allmot_v3_jx = 0
        allmot_v3_ev = 0
        for motif in _CANON_MOTIFS + [_OTHER_MOT]:
            v2_jxs = data[_JXS][_V2_JX][motif]
            v3_jxs = data[_JXS][_V3_JX][motif]
            mutual_jxs = data[_JXS][_MU_JX][motif]
            v2_events = data[_EVENTS][_V2_JX][motif]
            v3_events = data[_EVENTS][_V3_JX][motif]
            mutual_events = data[_EVENTS][_MU_JX][motif]
            allmot_v2_jx += v2_jxs
            allmot_v3_jx += v3_jxs
            allmot_mut_jx += mutual_jxs
            allmot_v2_ev += v2_events
            allmot_v3_ev += v3_events
            allmot_mut_ev += mutual_events

            if motif != _OTHER_MOT:
                try:
                    jx_recall3 = mutual_jxs / (mutual_jxs + v3_jxs)
                    event_recall3 = mutual_events / (mutual_events + v3_events)
                    plot_data_dict[motif]['data'][2].append(jx_recall3)
                    plot_data_dict[motif]['data'][3].append(event_recall3)
                    v3nonskips += 1
                except ZeroDivisionError:
                    v3skips += 1
                    v3skip_dict[motif] += 1
                    if motif == 'GTAG':
                        v3_gtag_blank.append(sample)
                try:
                    jx_recall2 = mutual_jxs / (mutual_jxs + v2_jxs)
                    event_recall2 = mutual_events / (mutual_events + v2_events)
                    plot_data_dict[motif]['data'][0].append(jx_recall2)
                    plot_data_dict[motif]['data'][1].append(event_recall2)
                    v2nonskips += 1
                except ZeroDivisionError:
                    v2skips += 1
                    v2skip_dict[motif] += 1

        jx_recall3 = allmot_mut_jx / (allmot_mut_jx + allmot_v3_jx)
        event_recall3 = allmot_mut_ev / (allmot_mut_ev + allmot_v3_ev)
        jx_recall2 = allmot_mut_jx / (allmot_mut_jx + allmot_v2_jx)
        event_recall2 = allmot_mut_ev / (allmot_mut_ev + allmot_v2_ev)
        plot_data_dict['all motifs']['data'][2].append(jx_recall3)
        plot_data_dict['all motifs']['data'][3].append(event_recall3)
        plot_data_dict['all motifs']['data'][0].append(jx_recall2)
        plot_data_dict['all motifs']['data'][1].append(event_recall2)

    fig_file = os.path.join(outpath, 'allmotifs_{}_{}.pdf'.format(cohort, now))
    grouped_boxplots(plot_data_dict, plot_info_dict, fig_file)
    print('v3: {} OK, {} skipped'.format(v3nonskips, v3skips))
    print(v3skip_dict)
    print(v3_gtag_blank)
    print('v2: {} OK, {} skipped'.format(v2nonskips, v2skips))
    print(v2skip_dict)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plots recall of junctions and junction events between '
                    'recount v2 and v3 for GTEx, TCGA, and SRA samples.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give the path to store junction comparison output and plots.'
    )
    parser.add_argument(
        '--sra-prepped-data', '-s',
        help='Provide the json file with SRA jx and event counts.'
    )
    parser.add_argument(
        '--tcga-prepped-data', '-t',
        help='Provide the json file with TCGA jx and event counts.'
    )
    parser.add_argument(
        '--gtex-prepped-data', '-g',
        help='Provide the json file with GTEx jx and event counts.'
    )

    args = parser.parse_args()
    out_path = args.output_path
    sra_json = args.sra_prepped_data
    tcga_json = args.tcga_prepped_data
    gtex_json = args.gtex_prepped_data

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    if gtex_json:
        prepare_boxplot_data(gtex_json, out_path, now, cohort=_GTEX)
    if tcga_json:
        prepare_boxplot_data(tcga_json, out_path, now, cohort=_TCGA)
    if sra_json:
        prepare_boxplot_data(sra_json, out_path, now, cohort=_SRA)
