#!/usr/bin/env python3

"""
recount_2v3_jx_comp.py
Python 3 code for comparing RNA-seq junctions called across mutual samples by
the recount2 vs. recount3 Sequence Read Archive, TCGA, and GTEx analyses.

"""

import argparse
import csv
from datetime import datetime
import gzip
import json
import logging
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.ticker as ticker


_LOG_MODE = 'INFO'
_ACCEPTABLE_MOTIFS = ['GTAG', 'GCAG', 'ATAC']
_MIN_COVS = [0, 5, 10]
_MULTIPLIER = 100000000
_OTHER_MOTIF = 'other'


class JXIndexError(Exception):
    pass


def grouped_boxplots(data_dict, plot_dict, fig_file, fig_size=(3.0, 5.0),
                     logscale=False, y_label='precision and recall',
                     percent=True, right_lim_shift=2,
                     x_label='minimum scaled coverage threshold'):
    """

    Input:
        data_dict: (dictionary) should have the form
        plot_dict: (dictionary) should have the form
    :param data_dict:

    :param out_path:
    :param now:
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
    labels = [str(key) for key in data_dict]
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
    ax.set_xlim(left=0, right=curr_label_loc-right_lim_shift)
    ax.spines['bottom'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['left'].set_color('black')

    if logscale:
        plt.yscale('log')

    plt.ylabel(y_label, fontsize=7)
    plt.xlabel(x_label, fontsize=7)
    plt.setp(
        ax.xaxis.get_majorticklabels(), rotation=90, fontsize=5, color='black'
    )
    plt.setp(
        ax.yaxis.get_majorticklabels(), fontsize=7, color='black'
    )
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
    plt.xticks([])
    fig = plt.gcf()
    fig.savefig(fig_file, bbox_inches='tight', pad_inches=.1)
    return


def filter_and_compare_jxs(sample_jx_dict, outpath, now, cohort=''):
    plot_info_dict = {
        'light colors': ['xkcd:tangerine', 'xkcd:cerulean'],
        'dark colors': ['xkcd:pumpkin', 'xkcd:ocean blue']
    }

    # Begin processing of all jxs/all motifs
    all_motif_data_dict = {cov: {'data': [[], []]} for cov in _MIN_COVS}
    for sample, data in sample_jx_dict:
        tot_cov = data['tot_cov']
        v2_covs = [
            _MULTIPLIER * cov / tot_cov for cov in data['v2_jxs'][_OTHER_MOTIF]
        ]
        v3_covs = [
            _MULTIPLIER * cov / tot_cov for cov in data['v3_jxs'][_OTHER_MOTIF]
        ]
        mut_covs = [
            _MULTIPLIER * cov / tot_cov for cov in data['mutual'][_OTHER_MOTIF]
        ]
        for threshold in _MIN_COVS:
            v2_thrsh = list(filter(lambda x: x >= threshold, v2_covs))
            v3_thrsh = list(filter(lambda x: x >= threshold, v3_covs))
            mut_thrsh = list(filter(lambda x: x >= threshold, mut_covs))
            precision = len(mut_thrsh) / (len(mut_thrsh) + len(v3_thrsh))
            recall = len(mut_thrsh) / (len(mut_thrsh) + len(v2_thrsh))
            all_motif_data_dict[threshold]['data'][0].append(precision)
            all_motif_data_dict[threshold]['data'][1].append(recall)

    # Process each canonical motif
    for motif in _ACCEPTABLE_MOTIFS:
        grouped_data_dict = {cov: {'data': [[], []]} for cov in _MIN_COVS}
        for sample, data in sample_jx_dict:
            tot_cov = data['tot_cov']
            v2_covs = [
                _MULTIPLIER * cov / tot_cov for cov in data['v2_jxs'][motif]
            ]
            v3_covs = [
                _MULTIPLIER * cov / tot_cov for cov in data['v3_jxs'][motif]
            ]
            mut_covs = [
                _MULTIPLIER * cov / tot_cov for cov in data['mutual'][motif]
            ]
            for threshold in _MIN_COVS:
                v2_thrsh = list(filter(lambda x: x >= threshold, v2_covs))
                v3_thrsh = list(filter(lambda x: x >= threshold, v3_covs))
                mut_thrsh = list(filter(lambda x: x >= threshold, mut_covs))
                precision = len(mut_thrsh) / (len(mut_thrsh) + len(v3_thrsh))
                recall = len(mut_thrsh) / (len(mut_thrsh) + len(v2_thrsh))
                grouped_data_dict[threshold]['data'][0].append(precision)
                grouped_data_dict[threshold]['data'][1].append(recall)
                all_motif_data_dict[threshold]['data'][0].append(precision)
                all_motif_data_dict[threshold]['data'][1].append(recall)

        fig_file = os.path.join(
            outpath, '{}_{}_{}.pdf'.format(cohort, motif, now)
        )
        data_file = os.path.join(
            outpath, '{}_{}_plotinfo_{}.json'.format(cohort, motif, now)
        )
        with open(data_file, 'w') as output:
            dump_list = [
                grouped_data_dict, plot_info_dict, fig_file
            ]
            json.dump(dump_list, output)

        grouped_boxplots(
            grouped_data_dict, plot_info_dict, fig_file, fig_size=(5.0, 5.0)
        )

    fig_file = os.path.join(outpath, '{}_allmotifs_{}.pdf'.format(cohort, now))
    data_file = os.path.join(
        outpath, '{}_allmotifs_plotinfo_{}.json'.format(cohort, now)
    )
    with open(data_file, 'w') as output:
        dump_list = [
            all_motif_data_dict, plot_info_dict, fig_file
        ]
        json.dump(dump_list, output)

    grouped_boxplots(
        all_motif_data_dict, plot_info_dict, fig_file, fig_size=(5.0, 5.0)
    )
    return


def map_recount_ids_2to3(recount2_sample_ids, recount3_id_map):
    sample_ids = pd.read_table(
        recount2_sample_ids, header=None, usecols=[0, 2],
        names=['recount_id', 'universal_id']
    )
    sample_ids['recount3'] = sample_ids.universal_id.apply(
        lambda x: recount3_id_map.get(x, 0)
    )
    return sample_ids.set_index('recount_id')['recount3'].to_dict()


def accession_to_recount3_ids(gtex_ids, tcga_ids, sra_ids):
    uppercase = {'gdc_file_id': lambda ident: ident.upper()}
    tcga_df = pd.read_table(
        tcga_ids, sep='\t', usecols=['rail_id', 'gdc_file_id'],
        converters=uppercase, dtype=str
    )
    tcga_df.rename(columns={'gdc_file_id': 'project_id'}, inplace=True)
    gtex_df = pd.read_table(
        gtex_ids, sep='\t', usecols=['rail_id', 'run_acc'], dtype=str
    )
    gtex_df.rename(columns={'run_acc': 'project_id'}, inplace=True)
    sra_df = pd.read_table(
        sra_ids, sep='\t', usecols=['rail_id', 'run_acc'], dtype=str
    )
    sra_df.rename(columns={'run_acc': 'project_id'}, inplace=True)
    full_df = pd.concat([tcga_df, gtex_df, sra_df])
    return full_df.set_index('project_id')['rail_id'].to_dict()


def intropolis_firstpass(jx_file, recount_2to3_map):
    sample_cov_dict = {}
    with gzip.open(jx_file) as cov_file:
        jx_cov = csv.reader(cov_file, delimiter='\t')
        for line in jx_cov:
            ids = line[6].split(',')
            covs = line[7].split(',')
            for id, cov in zip(ids, covs):
                id = recount_2to3_map[id]
                if not id:
                    continue
                try:
                    sample_cov_dict[id] += cov
                except KeyError:
                    sample_cov_dict[id] = cov

    return sample_cov_dict


def recount_firstpass(jx_file, recount_2to3_map=None):
    sample_cov_dict = {}
    with gzip.open(jx_file) as cov_file:
        print('starting tcga junctions')
        jx_cov = csv.reader(cov_file, delimiter='\t')
        for line in jx_cov:
            samp_cov_info = line[11].split(',')
            samp_cov_info.pop(0)
            for entry in samp_cov_info:
                id, cov = entry.split(':')
                if recount_2to3_map:
                    id = recount_2to3_map[id]
                if not id:
                    continue
                try:
                    sample_cov_dict[id] += cov
                except KeyError:
                    sample_cov_dict[id] = cov

    return sample_cov_dict


def intropolis_collect_sample_jxs(jx_file, mutual_samples, v2_id_cov):
    samp_jx_dict = {}
    with gzip.open(jx_file) as cov_file:
        print('starting tcga junctions')
        jx_cov = csv.reader(cov_file, delimiter='\t')
        for line in jx_cov:
            ids = line[6].split(',')
            if len(mutual_samples.intersection(set(ids))) == 0:
                continue
            covs = line[7].split(',')
            chr, left, right, strand = line[0], line[1], line[2], line[3]
            jx = ';'.join([chr, left, right, strand])
            motif = line[4] + line[5]
            if motif not in _ACCEPTABLE_MOTIFS:
                motif = _OTHER_MOTIF

            for id, cov in zip(ids, covs):
                if id not in mutual_samples:
                    continue
                scaled_cov = _MULTIPLIER * cov / v2_id_cov
                try:
                    samp_jx_dict[id][motif][jx] = scaled_cov
                except KeyError:
                    samp_jx_dict[id] = {
                        mot: {} for mot in _ACCEPTABLE_MOTIFS + [_OTHER_MOTIF]
                    }
                    samp_jx_dict[id][motif][jx] = scaled_cov
    return


def collect_mutual_jxs(v2_jx_file, v3_jx_file, v2_id_cov, v3_id_cov,
                       mutual_samples):
    samp_dict = {}
    with gzip.open(v2_jx_file) as v2_file, gzip.open(v3_jx_file) as v3_file:
        curr_line_2 = next(v2_file)
        chr2, l2, r2 = curr_line_2[1], curr_line_2[2], curr_line_2[3]
        curr_chr_2 = int(chr2.strip('chr'))
        curr_line_3= next(v3_file)
        chr3, l3, r3 = curr_line_3[1], curr_line_3[2], curr_line_3[3]
        curr_chr_3 = int(chr3.strip('chr'))

        # TODO: fix while loops
        while curr_chr_2 <= curr_chr_3:
            while l2 < l3:
                for line in v2_file:
                    print('a')
                    # TODO run through v2 and add junctions to non-v3 lists

            if l2 == l3:
                print('a')
                # TODO: run through right coordinates and add mutual jxs

            else:
                while l2 > l3:
                    for line in v3_file:
                        print('a')
                        # TODO run through v3 and add jxs to non-v2 list

    return samp_dict


def main():
    parser = argparse.ArgumentParser(
        description='Parses recount2 and recount3 junction files and compares'
                    'junctions between samples from both cohorts.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give the path to store junction comparison output and plots.'
    )
    parser.add_argument(
        '--sra-v2-jxs', '-s',
        help='Provide the file containing recount2 junctions from the SRA.'
    )
    parser.add_argument(
        '--sra-v3-jxs', '-S',
        help='Provide the file containing recount3 junctions from the SRA.'
    )
    parser.add_argument(
        '--tcga-v2-jxs', '-t',
        help='Provide the file containing recount2 junctions from TCGA.'
    )
    parser.add_argument(
        '--tcga-v3-jxs', '-T',
        help='Provide the file containing recount3 junctions from TCGA.'
    )
    parser.add_argument(
        '--gtex-v2-jxs', '-g',
        help='Provide the file containing recount2 junctions from GTEx.'
    )
    parser.add_argument(
        '--gtex-v3-jxs', '-G',
        help='Provide the file containing recount3 junctions from GTEx.'
    )
    parser.add_argument(
        '--recount2-sample-ids', '-r',
        help='Provide sample_ids.tsv file for mapping recount2 sample IDs to '
             'SRA, GTEx, and TCGA project IDs.'
    )
    parser.add_argument(
        '--tcga-phenotypes-v3', '-P',
        help='Provide the phenotype file for TCGA recount3 samples.'
    )
    parser.add_argument(
        '--gtex-phenotypes-v3', '-p',
        help='Provide the phenotype file for GTEx recount3 samples.'
    )
    parser.add_argument(
        '--sra-phenotypes-v3', '-R',
        help='Provide the phenotype file for SRA recount3 samples.'
    )

    args = parser.parse_args()
    out_path = args.output_path
    sra_v2 = args.sra_v2_jxs
    sra_v3 = args.sra_v3_jxs
    tcga_v2 = args.tcga_v2_jxs
    tcga_v3 = args.tcga_v3_jxs
    gtex_v2 = args.gtex_v2_jxs
    gtex_v3 = args.gtex_v3_jxs
    ids_v2 = args.recount2_sample_ids
    gtex_ids3 = args.gtex_phenotypes_v3
    tcga_ids3 = args.tcga_phenotyes_v3
    sra_ids3 = args.sra_phenotypes_v3

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'recount_2v3_jx_comparison_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=_LOG_MODE)

    # map recount2 sample IDs to recount3 sample IDs
    recount3_id_map = accession_to_recount3_ids(gtex_ids3, tcga_ids3, sra_ids3)
    recount2_id_map = map_recount_ids_2to3(ids_v2, recount3_id_map)

    # Sequence Read Archive
    # junction coordinates: intropolis and recount3 match
    v2_id_cov = intropolis_firstpass(sra_v2, recount2_id_map)
    v3_id_cov = recount_firstpass(sra_v3)
    mutual_samples = set(v2_id_cov.keys()).intersection(set(v3_id_cov.keys()))
    v2_avg_cov = sum(v2_id_cov.values()) / len(v2_id_cov)
    v3_avg_cov = sum(v3_id_cov.values()) / len(v3_id_cov)
    logging.info(
        'SRA: v2 samples: {}, v3 samples: {}, mutual samples: {}'
        ''.format(len(v2_id_cov), len(v3_id_cov), len(mutual_samples))
    )
    logging.info(
        'SRA avg. total coverage: v2={}, v3={}'.format(v2_avg_cov, v3_avg_cov)
    )
    intropolis_collect_sample_jxs(sra_v2, mutual_samples, v2_id_cov)
    # TODO: read through v2 and v3 files - together or consecutively?
    # TODO: sort: sample: junction, coverage
    sample_jx_dict = collect_mutual_jxs(
        sra_v2, sra_v3, v2_id_cov, v3_id_cov, mutual_samples
    )
    data_file = os.path.join(
        out_path, 'SRA_sample_jx_dict_{}.json'.format(now)
    )
    with open(data_file, 'w') as output:
        json.dump(sample_jx_dict, output)

    filter_and_compare_jxs(sample_jx_dict, out_path, now, cohort='SRA')

    # GTEx
    # junction coordinates: v2 left + 1 = v3; v2 right + 1 = v3
    v2_id_cov = recount_firstpass(gtex_v2, recount2_id_map)
    v3_id_cov = recount_firstpass(gtex_v3)
    mutual_samples = set(v2_id_cov.keys()).intersection(set(v3_id_cov.keys()))
    v2_avg_cov = sum(v2_id_cov.values()) / len(v2_id_cov)
    v3_avg_cov = sum(v3_id_cov.values()) / len(v3_id_cov)
    logging.info(
        'GTEx: v2 samples: {}, v3 samples: {}, mutual samples: {}'
        ''.format(len(v2_id_cov), len(v3_id_cov), len(mutual_samples))
    )
    logging.info(
        'GTEx avg. total coverage: v2={}, v3={}'.format(v2_avg_cov, v3_avg_cov)
    )
    sample_jx_dict = collect_mutual_jxs(
        gtex_v2, gtex_v3, v2_id_cov, v3_id_cov, mutual_samples
    )
    data_file = os.path.join(
        out_path, 'GTEx_sample_jx_dict_{}.json'.format(now)
    )
    with open(data_file, 'w') as output:
        json.dump(sample_jx_dict, output)

    filter_and_compare_jxs(sample_jx_dict, out_path, now, cohort='GTEx')

    # TCGA
    # junction coordinates: v2 left + 1 = v3; v2 right + 1 = v3
    v2_id_cov = recount_firstpass(tcga_v2, recount2_id_map)
    v3_id_cov = recount_firstpass(tcga_v3)
    mutual_samples = set(v2_id_cov.keys()).intersection(set(v3_id_cov.keys()))
    v2_avg_cov = sum(v2_id_cov.values()) / len(v2_id_cov)
    v3_avg_cov = sum(v3_id_cov.values()) / len(v3_id_cov)
    logging.info(
        'TCGA: v2 samples: {}, v3 samples: {}, mutual samples: {}'
        ''.format(len(v2_id_cov), len(v3_id_cov), len(mutual_samples))
    )
    logging.info(
        'TCGA avg. total coverage: v2={}, v3={}'.format(v2_avg_cov, v3_avg_cov)
    )
    sample_jx_dict = collect_mutual_jxs(
        tcga_v2, tcga_v3, v2_id_cov, v3_id_cov, mutual_samples
    )
    data_file = os.path.join(
        out_path, 'TCGA_sample_jx_dict_{}.json'.format(now)
    )
    with open(data_file, 'w') as output:
        json.dump(sample_jx_dict, output)

    filter_and_compare_jxs(sample_jx_dict, out_path, now, cohort='TCGA')


if __name__ == '__main__':
    main()
