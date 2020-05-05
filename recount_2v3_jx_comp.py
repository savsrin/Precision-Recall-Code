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
import matplotlib.patches as mpatches
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


class JXIndexError(Exception):
    pass


def grouped_boxplots(data_dict, plot_dict, fig_file, fig_size=(7.0, 4.0),
                     logscale=False, y_label= 'recall', percent=True,
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

    legend_list = []
    for p_col, p_lab in zip(plot_dict['dark colors'], plot_dict['row labels']):
        legend_list.append(mpatches.Patch(color=p_col, label=p_lab))

    plt.legend(handles=legend_list)

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


def prepare_boxplot_data(sample_jx_dict, outpath, now, cohort):
    plot_info_dict = {
        'light colors': ['#29A1E8', '#85C9F2', '#E8321F', '#F2867B'],
        'dark colors':  ['#29A1E8', '#85C9F2', '#E8321F', '#F2867B'],
        'row labels': [
            'v3 junction recall: v2 ground truth',
            'v3 detection event recall; v2 ground truth',
            'v2 junction recall; v3 ground truth',
            'v2 detection event recall; v2 ground truth'
        ]
    }

    # Begin processing of all jxs/all motifs
    plot_data_dict = {mot: {'data': [[], [], [], []]} for mot in _CANON_MOTIFS}
    plot_data_dict['all motifs'] = {'data': [[], [], [], []]}

    # Process each canonical motif
    for motif in _CANON_MOTIFS + [_OTHER_MOT]:
        for sample, data in sample_jx_dict.items():
            v2_jxs = data[_JXS][_V2_JX][motif]
            v3_jxs = data[_JXS][_V3_JX][motif]
            mutual_jxs = data[_JXS][_MU_JX][motif]
            v2_events = data[_EVENTS][_V2_JX][motif]
            v3_events = data[_EVENTS][_V3_JX][motif]
            mutual_events = data[_EVENTS][_MU_JX][motif]

            jx_recall3 = mutual_jxs / (mutual_jxs + v3_jxs)
            jx_recall2 = mutual_jxs / (mutual_jxs + v2_jxs)
            event_recall3 = mutual_events / (mutual_events + v3_events)
            event_recall2 = mutual_events / (mutual_events + v2_events)

            plot_data_dict['all motifs']['data'][0].append(jx_recall2)
            plot_data_dict['all motifs']['data'][1].append(event_recall2)
            plot_data_dict['all motifs']['data'][2].append(jx_recall3)
            plot_data_dict['all motifs']['data'][3].append(event_recall3)

            if motif != _OTHER_MOT:
                plot_data_dict[motif]['data'][0].append(jx_recall2)
                plot_data_dict[motif]['data'][1].append(event_recall2)
                plot_data_dict[motif]['data'][2].append(jx_recall3)
                plot_data_dict[motif]['data'][3].append(event_recall3)

    fig_file = os.path.join(outpath, '{}_allmotifs_{}.pdf'.format(cohort, now))
    data_file = os.path.join(
        outpath, '{}_allmotifs_plotinfo_{}.json'.format(cohort, now)
    )
    with open(data_file, 'w') as output:
        dump_list = [
            plot_data_dict, plot_info_dict, fig_file
        ]
        json.dump(dump_list, output)

    grouped_boxplots(
        plot_data_dict, plot_info_dict, fig_file
    )
    return


def map_recount_ids_2to3(recount2_sample_ids, recount3_id_map):
    sample_ids = pd.read_table(
        recount2_sample_ids, header=None, usecols=[0, 2],
        names=['recount_id', 'universal_id'], dtype={'recount_id': str}
    )
    sample_ids['recount3'] = sample_ids.universal_id.apply(
        lambda x: recount3_id_map.get(x, 0)
    )
    return sample_ids.set_index('recount_id')['recount3'].to_dict()


def accession_to_recount3_ids(gtex_ids, tcga_ids, sra_ids):
    uppercase = {'gdc_file_id': lambda ident: ident.upper()}
    tcga_df = pd.read_table(
        tcga_ids, sep='\t', usecols=['rail_id', 'gdc_file_id'],
        converters=uppercase, dtype={'rail_id': str}
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
    sample_set = set()
    with gzip.open(jx_file, 'rt') as cov_file:
        jx_cov = csv.reader(cov_file, delimiter='\t')
        for line in jx_cov:
            ids = line[6].split(',')
            for samp in ids:
                samp = recount_2to3_map[samp]
                if not samp:
                    continue
                sample_set.add(samp)

    return sample_set


def recount_firstpass(jx_file, recount_2to3_map=None):
    sample_set = set()
    with gzip.open(jx_file, 'rt') as cov_file:
        jx_cov = csv.reader(cov_file, delimiter='\t')
        for line in jx_cov:
            samp_cov_info = line[11].split(',')
            samp_cov_info.pop(0)
            for entry in samp_cov_info:
                samp, cov = entry.split(':')
                if recount_2to3_map:
                    samp = recount_2to3_map.get(samp, None)
                if not samp:
                    continue
                sample_set.add(samp)

    return sample_set


def coordinates_from_jx_line(jx_line, coord_positions):
    """

    :param jx_line:
    :return:
    """
    if jx_line == None:
        return None
    chrom = jx_line[coord_positions[0]]
    left = int(jx_line[coord_positions[1]])
    right = int(jx_line[coord_positions[2]])
    strand = jx_line[coord_positions[3]]
    return (chrom, left, right, strand)


def add_mutual_jx_to_samp_dict(samp_dict, line2, line3, mutual_samples,
                               recount2_id_map, flag):
    """

    :param samp_dict:
    :param line2:
    :param line3:
    :param mutual_samples:
    :param recount2_id_map:
    :return:
    """
    joint_samps = {}
    if type(line3[0]) == str:
        motif = line3[7] + line3[8]
        strand = line3[5]
        samp_cov_3_list = line3[11].split(',')
        samp_cov_3_list.pop(0)
    else:
        motif = line3[0][7] + line3[0][8]
        strand = line3[0][5]
        samp_cov_3_list = []
        for line in line3:
            temp_list = line[11].split(',')
            temp_list.pop(0)
            samp_cov_3_list.extend(temp_list)

    if strand == '-':
        motif = motif.translate(_COMPLEMENT)[::-1]

    if flag == _SRA:
        v2_samps = line2[6].split(',')
        v2_covs = line2[7].split(',')
        for samp, cov in zip(v2_samps, v2_covs):
            samp = recount2_id_map.get(samp, None)
            if samp not in mutual_samples:
                continue
            if samp in samp_cov_3_list:
                joint_samps[samp] = {'v2': int(cov)}
            else:
                samp_dict[samp][_JXS][_V2_JX][motif] += 1
                samp_dict[samp][_EVENTS][_V2_JX][motif] += int(cov)
    else:
        samp_cov_2_list = line2[11].split(',')
        samp_cov_2_list.pop(0)
        for samp_info in samp_cov_2_list:
            samp, cov = samp_info.split(':')
            samp = recount2_id_map.get(samp, None)
            if samp not in mutual_samples:
                continue
            if samp in samp_cov_3_list:
                joint_samps[samp] = {'v2': int(cov)}
            else:
                samp_dict[samp][_JXS][_V2_JX][motif] += 1
                samp_dict[samp][_EVENTS][_V2_JX][motif] += int(cov)

    # process v3 samples and coverage
    for samp_info in samp_cov_3_list:
        samp, cov = samp_info.split(':')
        if samp not in mutual_samples:
            continue
        if samp in joint_samps.keys():
            joint_samps[samp]['v3'] = int(cov)
        else:
            samp_dict[samp][_JXS][_V3_JX][motif] += 1
            samp_dict[samp][_EVENTS][_V3_JX][motif] += int(cov)

    # process coverage in mutual samples/jxs
    for samp, values in joint_samps.items():
        v2_count = values['v2']
        v3_count = values['v3']
        if v2_count < v3_count:
            mutual_events = v2_count
            samp_dict[samp][_EVENTS][_V3_JX][motif] += v3_count - v2_count
        else:
            mutual_events = v3_count
            samp_dict[samp][_EVENTS][_V2_JX][motif] += v2_count - v3_count

        samp_dict[samp][_JXS][_MU_JX][motif] += 1
        samp_dict[samp][_EVENTS][_MU_JX][motif] += mutual_events

    return samp_dict


def add_v2_jx_to_samp_dict(samp_dict, jx_line, mutual_samples, recount2_ids,
                           flag):
    """

    :param samp_dict:
    :param jx_line:
    :param mutual_samples:
    :param recount2_ids:
    :return:
    """
    if flag == _SRA:
        samp_dict = add_intropolis_jx_to_samp_dict(
            samp_dict, jx_line, mutual_samples, recount2_ids
        )
    else:
        motif = jx_line[7] + jx_line[8]
        samp_dict = add_recount_jx_to_samp_dict(
            samp_dict, jx_line, _V2_JX, motif, mutual_samples
        )
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

    samp_dict = add_recount_jx_to_samp_dict(
        samp_dict, jx_line, _V3_JX, motif, mutual_samps
    )
    return samp_dict


def add_intropolis_jx_to_samp_dict(samp_dict, jx_line, mutual_samples,
                                   recount2_ids):
    motif = jx_line[4] + jx_line[5]
    samps = jx_line[6].split(',')
    covs = jx_line[7].split(',')
    for samp, cov in zip(samps, covs):
        samp = recount2_ids.get(samp, None)
        if samp not in mutual_samples:
            continue
        if motif not in _CANON_MOTIFS:
            motif = _OTHER_MOT
        samp_dict[samp][_JXS][_V2_JX][motif] += 1
        samp_dict[samp][_EVENTS][_V2_JX][motif] += int(cov)
    return samp_dict


def add_recount_jx_to_samp_dict(samp_dict, jx_line, v2v3, motif, mutual_samps):
    samp_cov_list = jx_line[11].split(',')
    samp_cov_list.pop(0)
    for samp_info in samp_cov_list:
        samp, cov = samp_info.split(':')
        if samp not in mutual_samps:
            continue
        if motif not in _CANON_MOTIFS:
            motif = _OTHER_MOT
        samp_dict[samp][_JXS][v2v3][motif] += 1
        samp_dict[samp][_EVENTS][v2v3][motif] += int(cov)
    return samp_dict


def collect_jx_covs(v2_jxs, v3_jxs, mutual_samples, recount2_id_map,
                    cohort_flag):
    """

    :param v2_jxs:
    :param v3_jxs:
    :param v2_id_cov:
    :param v3_id_cov:
    :param mutual_samples:
    :param recount2_id_map:
    :return:
    """
    if cohort_flag == _SRA:
        coord_locs = (0, 1, 2, 3)
    else:
        coord_locs = (1, 2, 3, 5)

    samp_dict = {}
    for samp in mutual_samples:
        samp_dict[samp] = {
            _JXS: {
                _V2_JX: {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]},
                _V3_JX: {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]},
                _MU_JX: {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]}
            },
            _EVENTS: {
                _V2_JX: {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]},
                _V3_JX: {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]},
                _MU_JX: {mot: 0 for mot in _CANON_MOTIFS + [_OTHER_MOT]}
            }
        }

    with gzip.open(v2_jxs, 'rt') as v2file, gzip.open(v3_jxs, 'rt') as v3file:
        file2 = csv.reader(v2file, delimiter='\t')
        file3 = csv.reader(v3file, delimiter='\t')
        v2_line = next(file2)
        v3_line = next(file3)
        v2_coords = coordinates_from_jx_line(v2_line, coord_locs)
        v3_coords = coordinates_from_jx_line(v3_line, coord_locs)
        while v2_line or v3_line:
            # Finish off v3 or v2 file after the other ends:
            if not v2_line:
                st3 = v3_coords[3]
                samp_dict = add_v3_jx_to_samp_dict(
                    samp_dict, v3_line, st3, mutual_samples
                )
                v3_line = next(file3, None)
                v3_coords = coordinates_from_jx_line(v3_line, coord_locs)
                continue

            if not v3_line:
                samp_dict = add_v2_jx_to_samp_dict(
                    samp_dict, v2_line, mutual_samples, recount2_id_map
                )
                v2_line = next(file2, None)
                v2_coords = coordinates_from_jx_line(v2_line, coord_locs)
                continue 

            # Main comparison: running through lines in both files
            if v2_coords < v3_coords:
                samp_dict = add_v2_jx_to_samp_dict(
                    samp_dict, v2_line, mutual_samples, recount2_id_map
                )
                v2_line = next(file2, None)
                v2_coords = coordinates_from_jx_line(v2_line, coord_locs)
                continue
            elif v2_coords > v3_coords:
                st3 = v3_coords[3]
                samp_dict = add_v3_jx_to_samp_dict(
                    samp_dict, v3_line, st3, mutual_samples
                )
                v3_line = next(file3, None)
                v3_coords = coordinates_from_jx_line(v3_line, coord_locs)
                continue
            else:
                # add jx to mutual
                v3_init_coords = v3_coords
                v3_init_line = v3_line
                v3_line_list = []
                while v3_coords == v3_init_coords:
                    v3_line_list.append(v3_line)
                    v3_line = next(file3, None)
                    v3_coords = coordinates_from_jx_line(v3_line, coord_locs)

                if len(v3_line_list) == 1:
                    samp_dict = add_mutual_jx_to_samp_dict(
                        samp_dict, v2_line, v3_init_line, mutual_samples,
                        recount2_id_map
                    )
                else:
                    samp_dict = add_mutual_jx_to_samp_dict(
                        samp_dict, v2_line, v3_line_list, mutual_samples,
                        recount2_id_map
                    )

                v2_line = next(file2, None)
                v2_coords = coordinates_from_jx_line(v2_line, coord_locs)
                continue

    return samp_dict


def execute_jx_comp(v2_file, v3_file, recount2_id_map, out_path, now, flag,
                    mutual_sample_file=None):
    if not mutual_sample_file:
        if flag == _SRA:
            v2_ids = intropolis_firstpass(v2_file, recount2_id_map)

        else:
            v2_ids = recount_firstpass(v2_file, recount2_id_map)

        v3_ids = recount_firstpass(v3_file)
        mutual_samples = v2_ids.intersection(v3_ids)

        logging.info(
            '{}: v2 samples: {}, v3 samples: {}'
            ''.format(flag, len(v2_ids), len(v3_ids))
        )
        sample_file = os.path.join(
            out_path, '{}_mutual_samples_{}.json'.format(flag, now)
        )
        with open(sample_file, 'w') as output:
            json.dump(list(mutual_samples), output)
    else:
        with open(mutual_sample_file) as recover:
            mutual_samples = set(json.load(recover))

    logging.info(
        '{}: {} mutual samples'.format(flag, len(mutual_samples))
    )    
    sample_jx_dict = collect_jx_covs(
        v2_file, v3_file, mutual_samples, recount2_id_map, flag
    )
    data_file = os.path.join(
        out_path, '{}_sample_jx_dict_{}.json'.format(flag, now)
    )
    with open(data_file, 'w') as output:
        json.dump(sample_jx_dict, output)

    prepare_boxplot_data(sample_jx_dict, out_path, now, cohort=flag)
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
    parser.add_argument(
        '--gtex-mutual-samples',
        help='Provide a json file of GTEx samples in recount2 and recount3.'
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
    tcga_ids3 = args.tcga_phenotypes_v3
    sra_ids3 = args.sra_phenotypes_v3
    gtex_both = args.gtex_mutual_samples

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'recount_2v3_jx_comparison_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=_LOG_MODE)

    # map recount2 sample IDs to recount3 sample IDs
    recount3_id_map = accession_to_recount3_ids(gtex_ids3, tcga_ids3, sra_ids3)
    recount2_id_map = map_recount_ids_2to3(ids_v2, recount3_id_map)
    
    # compare v2 to v3 jxs for three cohorts
    execute_jx_comp(
        gtex_v2, gtex_v3, recount2_id_map, out_path, now, _GTEX, gtex_both
    )    
    execute_jx_comp(tcga_v2, tcga_v3, recount2_id_map, out_path, now, _TCGA)
    execute_jx_comp(sra_v2, sra_v3, recount2_id_map, out_path, now, _SRA)

