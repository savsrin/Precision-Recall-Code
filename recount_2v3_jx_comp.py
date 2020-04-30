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
import matplotlib.ticker as ticker
import os
import pandas as pd
import seaborn as sns
import sys

csv.field_size_limit(sys.maxsize)
sns.set(color_codes=True)

_LOG_MODE = 'INFO'
_CANON_MOTIFS = ['GTAG', 'GCAG', 'ATAC']
_MIN_COVS = [0, 5, 10]
_MULTIPLIER = 100000000
_OTHER_MOT = 'other'
_V2_JX = 'v2_jxs'
_V3_JX = 'v3_jxs'
_MU_JX = 'mutual'
_TOT_COV2 = 'total_v2_coverage'
_TOT_COV3 = 'total_v3_coverage'
_COMPLEMENT = str.maketrans("ATCG", "TAGC")


class JXIndexError(Exception):
    pass


def grouped_boxplots(data_dict, plot_dict, fig_file, fig_size=(3.0, 5.0),
                     logscale=False, y_label='precision and recall',
                     percent=True, right_lim_shift=2,
                     x_label='minimum scaled coverage threshold'):
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
        'light colors': ['#ED5C4D', '#57B5ED'],
        'dark colors': ['#ED5C4D', '#57B5ED']
    }
    # Begin processing of all jxs/all motifs
    all_motif_data_dict = {cov: {'data': [[], []]} for cov in _MIN_COVS}
    for sample, data in sample_jx_dict:
        v2_tot = data[_TOT_COV2]
        v3_tot = data[_TOT_COV3]
        v2_covs = [
            _MULTIPLIER * cov / v2_tot for cov in data[_V2_JX][_OTHER_MOT]
        ]
        v3_covs = [
            _MULTIPLIER * cov / v3_tot for cov in data[_V3_JX][_OTHER_MOT]
        ]
        v2_mutual_scaled = [
            _MULTIPLIER * cov / v2_tot for cov in data[_MU_JX][_OTHER_MOT][0]
        ]
        v3_mutual_scaled = [
            _MULTIPLIER * cov / v3_tot for cov in data[_MU_JX][_OTHER_MOT][1]
        ]

        for threshold in _MIN_COVS:
            v2_count = len(list(filter(lambda x: x >= threshold, v2_covs)))
            v3_count = len(list(filter(lambda x: x >= threshold, v3_covs)))
            mutual_count = 0
            for v2_val, v3_val in zip(v2_mutual_scaled, v3_mutual_scaled):
                if v2_val >= threshold:
                    if v3_val >= threshold:
                        mutual_count += 1
                    else:
                        v2_count += 1
                else:
                    if v3_val >= threshold:
                        v3_count += 1

            precision = mutual_count / (mutual_count + v3_count)
            recall = mutual_count / (mutual_count + v2_count)
            all_motif_data_dict[threshold]['data'][0].append(precision)
            all_motif_data_dict[threshold]['data'][1].append(recall)

    # Process each canonical motif
    for motif in _CANON_MOTIFS:
        grouped_data_dict = {cov: {'data': [[], []]} for cov in _MIN_COVS}
        for sample, data in sample_jx_dict:
            v2_tot = data[_TOT_COV2]
            v3_tot = data[_TOT_COV3]
            v2_covs = [
                _MULTIPLIER * cov / v2_tot for cov in data[_V2_JX][motif]
            ]
            v3_covs = [
                _MULTIPLIER * cov / v3_tot for cov in data[_V3_JX][motif]
            ]
            v2_mutual_scaled = [
                _MULTIPLIER * cov / v2_tot for cov in data[_MU_JX][motif][0]
            ]
            v3_mutual_scaled = [
                _MULTIPLIER * cov / v3_tot for cov in data[_MU_JX][motif][1]
            ]

            for threshold in _MIN_COVS:
                v2_count = len(list(filter(lambda x: x >= threshold, v2_covs)))
                v3_count = len(list(filter(lambda x: x >= threshold, v3_covs)))
                mutual_count = 0
                for v2_val, v3_val in zip(v2_mutual_scaled, v3_mutual_scaled):
                    if v2_val >= threshold:
                        if v3_val >= threshold:
                            mutual_count += 1
                        else:
                            v2_count += 1
                    else:
                        if v3_val >= threshold:
                            v3_count += 1

                precision = mutual_count / (mutual_count + v3_count)
                recall = mutual_count / (mutual_count + v2_count)
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
    with gzip.open(jx_file, 'rt') as cov_file:
        jx_cov = csv.reader(cov_file, delimiter='\t')
        for line in jx_cov:
            ids = line[6].split(',')
            covs = line[7].split(',')
            for samp, cov in zip(ids, covs):
                samp = recount_2to3_map[samp]
                if not samp:
                    continue
                try:
                    sample_cov_dict[samp] += cov
                except KeyError:
                    sample_cov_dict[samp] = cov

    return sample_cov_dict


def recount_firstpass(jx_file, recount_2to3_map=None):
    sample_cov_dict = {}
    with gzip.open(jx_file, 'rt') as cov_file:
        print('starting tcga junctions')
        jx_cov = csv.reader(cov_file, delimiter='\t')
        for line in jx_cov:
            samp_cov_info = line[11].split(',')
            samp_cov_info.pop(0)
            for entry in samp_cov_info:
                samp, cov = entry.split(':')
                if recount_2to3_map:
                    samp = recount_2to3_map[samp]
                if not samp:
                    continue
                try:
                    sample_cov_dict[samp] += cov
                except KeyError:
                    sample_cov_dict[samp] = cov

    return sample_cov_dict


def coordinates_from_recount_line(recount_line):
    """

    :param recount_line:
    :return:
    """
    chrom = recount_line[1]
    chrom_int = int(chrom.strip('chr'))
    left = int(recount_line[2])
    right = int(recount_line[3])
    strand = recount_line[5]
    return (chrom_int, left, right, strand)


def add_mutual_jx_to_samp_dict(samp_dict, line2, line3, mutual_samples,
                               recount2_id_map):
    """

    :param samp_dict:
    :param line2:
    :param line3:
    :param mutual_samples:
    :param recount2_id_map:
    :return:
    """
    joint_samps = {}
    motif = line2[4] + line2[5]
    v2_samps = line2[6].split(',')
    v2_covs = line2[7].split(',')
    samp_cov_3_list = line3[11].split(',')

    # process v2 samples and coverage
    for samp, cov in zip(v2_samps, v2_covs):
        samp = recount2_id_map.get(samp, None)
        if samp not in mutual_samples:
            continue
        if samp in samp_cov_3_list:
            joint_samps[samp] = {'v2': cov}
        else:
            samp_dict[samp][_V2_JX][motif].append(cov)

    # process v3 samples and coverage
    for samp_info in samp_cov_3_list:
        samp, cov = samp_info.split(':')
        if samp not in mutual_samples:
            continue
        if samp in joint_samps.keys():
            joint_samps[samp]['v3'] = cov
        else:
            samp_dict[samp][_V3_JX][motif].append(cov)

    # process coverage in mutual samples/jxs
    for samp, values in joint_samps.items():
        samp_dict[samp][_MU_JX][motif][0].append(values['v2'])
        samp_dict[samp][_MU_JX][motif][1].append(values['v3'])

    return samp_dict


def add_v2_jx_to_samp_dict(samp_dict, jx_line, mutual_samples, recount2_ids):
    """

    :param samp_dict:
    :param jx_line:
    :param mutual_samples:
    :param recount2_ids:
    :return:
    """
    motif = jx_line[4] + jx_line[5]
    samps = jx_line[6].split(',')
    covs = jx_line[7].split(',')
    for samp, cov in zip(samps, covs):
        samp = recount2_ids.get(samp, None)
        if samp not in mutual_samples:
            continue
        if motif not in _CANON_MOTIFS:
            motif = _OTHER_MOT
        samp_dict[samp][_V2_JX][motif].append(cov)

    return samp_dict


def add_v3_jx_to_samp_dict(samp_dict, jx_line, strand, mutual_samps):
    """

    :param samp_dict:
    :param jx_line:
    :param strand:
    :param mutual_samps:
    :return:
    """
    motif = jx_line[7] + jx_line[8]
    if strand == '-':
        motif = motif.translate(_COMPLEMENT)[::-1]
    samp_cov_list = jx_line[11].split(',')
    for samp_info in samp_cov_list:
        samp, cov = samp_info.split(':')
        if samp not in mutual_samps:
            continue
        if motif not in _CANON_MOTIFS:
            motif = _OTHER_MOT
        samp_dict[samp][_V3_JX][motif].append(cov)

    return samp_dict


def collect_jx_covs(v2_jxs, v3_jxs, v2_id_cov, v3_id_cov, mutual_samples,
                    recount2_id_map):
    """

    :param v2_jxs:
    :param v3_jxs:
    :param v2_id_cov:
    :param v3_id_cov:
    :param mutual_samples:
    :param recount2_id_map:
    :return:
    """
    samp_dict = {}
    for samp in mutual_samples:
        samp_dict[samp] = {
            _TOT_COV2: v2_id_cov[samp],
            _TOT_COV3: v3_id_cov[samp],
            _V2_JX: {mot: [] for mot in _CANON_MOTIFS + [_OTHER_MOT]},
            _V3_JX: {mot: [] for mot in _CANON_MOTIFS + [_OTHER_MOT]},
            _MU_JX: {mot: [[], []] for mot in _CANON_MOTIFS + [_OTHER_MOT]}
        }

    with gzip.open(v2_jxs, 'rt') as file2, gzip.open(v3_jxs, 'rt') as file3:
        v2_line = next(file2)
        v3_line = next(file3)
        v2_tup = coordinates_from_recount_line(v2_line)
        v3_tup = coordinates_from_recount_line(v3_line)
        # TODO: add coordinates from recount2
        while v2_line or v3_line:
            # Finish off v3 or v2 file after the other ends:
            if not v2_line:
                st3 = v3_tup[3]
                samp_dict = add_v3_jx_to_samp_dict(
                    samp_dict, v3_line, st3, mutual_samples
                )
                v3_line = next(file3)
                v3_tup = coordinates_from_recount_line(v3_line)
                continue

            if not v3_line:
                samp_dict = add_v2_jx_to_samp_dict(
                    samp_dict, v2_line, mutual_samples, recount2_id_map
                )
                v2_line = next(file2)
                v2_tup = coordinates_from_recount_line(v2_line)
                continue

            # Main comparison: running through lines in both files
            if v2_tup < v3_tup:
                samp_dict = add_v2_jx_to_samp_dict(
                    samp_dict, v2_line, mutual_samples, recount2_id_map
                )
                v2_line = next(file2)
                v2_tup = coordinates_from_recount_line(v2_line)
                continue
            elif v2_tup > v3_tup:
                st3 = v3_tup[3]
                samp_dict = add_v3_jx_to_samp_dict(
                    samp_dict, v3_line, st3, mutual_samples
                )
                v3_line = next(file3)
                v3_tup = coordinates_from_recount_line(v3_line)
            else:
                # add jx to mutual
                # TODO: check for next-line matching jxs
                # TODO: improve addition of jxs to dictionary
                # TODO: add "junction
                v2_line = next(file2)
                v2_tup = coordinates_from_recount_line(v2_line)
                v3_line = next(file3)
                v3_tup = coordinates_from_recount_line(v3_line)
                continue

    return samp_dict


def execute_jx_comp(v2_file, v3_file, recount2_id_map, out_path, now, flag):


    if flag == 'SRA':
        v2_id_cov = intropolis_firstpass(v2_file, recount2_id_map)

    else:
        v2_id_cov = recount_firstpass(v2_file, recount2_id_map)

    v3_id_cov = recount_firstpass(v3_file)
    mutual_samples = set(v2_id_cov.keys()).intersection(set(v3_id_cov.keys()))
    unwanted_keys_2 = set(v2_id_cov.keys) - mutual_samples
    for key in unwanted_keys_2:
        del v2_id_cov[key]

    unwanted_keys_3 = set(v3_id_cov.keys()) - mutual_samples
    for key in unwanted_keys_3:
        del v3_id_cov[key]

    v2_avg_cov = sum(v2_id_cov.values()) / len(v2_id_cov)
    v3_avg_cov = sum(v3_id_cov.values()) / len(v3_id_cov)
    logging.info(
        '{}: v2 samples: {}, v3 samples: {}, mutual samples: {}'
        ''.format(flag, len(v2_id_cov), len(v3_id_cov), len(mutual_samples))
    )
    logging.info(
        '{} avg. total coverage: v2={}, v3={}'
        ''.format(flag, v2_avg_cov, v3_avg_cov)
    )
    sample_jx_dict = collect_jx_covs(
        v2_file, v3_file, v2_id_cov, v3_id_cov, mutual_samples, recount2_id_map
    )
    data_file = os.path.join(
        out_path, '{}_sample_jx_dict_{}.json'.format(flag, now)
    )
    with open(data_file, 'w') as output:
        json.dump(sample_jx_dict, output)

    filter_and_compare_jxs(sample_jx_dict, out_path, now, cohort='SRA')
    return


if __name__ == '__main__':
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

    # compare v2 to v3 jxs for three cohorts
    execute_jx_comp(gtex_v2, gtex_v3, recount2_id_map, out_path, now, 'GTEx')
    execute_jx_comp(tcga_v2, tcga_v3, recount2_id_map, out_path, now, 'TCGA')
    execute_jx_comp(sra_v2, sra_v3, recount2_id_map, out_path, now, 'SRA')
