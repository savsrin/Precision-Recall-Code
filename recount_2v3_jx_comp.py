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
import os
import pandas as pd
import sys
import time

csv.field_size_limit(sys.maxsize)

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


def coords_from_jx_line(jx_line, coord_positions):
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

    v3_samp_set = set()
    for samp_info in samp_cov_3_list:
        samp, _ = samp_info.split(':')
        v3_samp_set.add(samp)

    if strand == '-':
        motif = motif.translate(_COMPLEMENT)[::-1]

    if flag == _SRA:
        v2_samps = line2[6].split(',')
        v2_covs = line2[7].split(',')
        for samp, cov in zip(v2_samps, v2_covs):
            samp = recount2_id_map.get(samp, None)
            if samp not in mutual_samples:
                continue
            if samp in v3_samp_set:
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
            if samp in v3_samp_set:
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
                    cohort_flag, out_path):
    """

    :param v2_jxs:
    :param v3_jxs:
    :param v2_id_cov:
    :param v3_id_cov:
    :param mutual_samples:
    :param recount2_id_map:
    :return:
    """
    v3_coord_locs = (1, 2, 3, 5)
    if cohort_flag == _SRA:
        v2_coord_locs = (0, 1, 2, 3)
    else:
        v2_coord_locs = (1, 2, 3, 5)
    ticker = 0
    t_int = time.time()
    start_time = t_int
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
        v2_coords = coords_from_jx_line(v2_line, v2_coord_locs)
        v3_coords = coords_from_jx_line(v3_line, v3_coord_locs)
        while v2_line or v3_line:
            ticker += 1
            if (ticker % 10000000) == 0:
                t_curr = time.time()
                logging.info(
                    '{}th jx, intermediate time is {}'
                    ''.format(ticker, t_curr - t_int)
                )
                logging.info('total time is {}'.format(t_curr - start_time))
                t_int = time.time() 
                int_file = os.path.join(
                    out_path, '{}_int_dict_{}.json'.format(cohort_flag, ticker)
                )
                with open(int_file, 'w') as output:
                    json.dump(samp_dict, output)

            # Finish off v3 or v2 file after the other ends:
            if not v2_line:
                st3 = v3_coords[3]
                samp_dict = add_v3_jx_to_samp_dict(
                    samp_dict, v3_line, st3, mutual_samples
                )
                v3_line = next(file3, None)
                v3_coords = coords_from_jx_line(v3_line, v3_coord_locs)
                continue

            if not v3_line:
                samp_dict = add_v2_jx_to_samp_dict(
                    samp_dict, v2_line, mutual_samples, recount2_id_map, 
                    cohort_flag
                )
                v2_line = next(file2, None)
                v2_coords = coords_from_jx_line(v2_line, v2_coord_locs)
                continue 

            # Main comparison: running through lines in both files
            if v2_coords < v3_coords:
                samp_dict = add_v2_jx_to_samp_dict(
                    samp_dict, v2_line, mutual_samples, recount2_id_map,
                    cohort_flag
                )
                v2_line = next(file2, None)
                v2_coords = coords_from_jx_line(v2_line, v2_coord_locs)
                continue
            elif v2_coords > v3_coords:
                st3 = v3_coords[3]
                samp_dict = add_v3_jx_to_samp_dict(
                    samp_dict, v3_line, st3, mutual_samples
                )
                v3_line = next(file3, None)
                v3_coords = coords_from_jx_line(v3_line, v3_coord_locs)
                continue
            else:
                # add jx to mutual
                v3_init_coords = v3_coords
                v3_init_line = v3_line
                v3_line_list = []
                while v3_coords == v3_init_coords:
                    v3_line_list.append(v3_line)
                    v3_line = next(file3, None)
                    v3_coords = coords_from_jx_line(v3_line, v3_coord_locs)

                if len(v3_line_list) == 1:
                    samp_dict = add_mutual_jx_to_samp_dict(
                        samp_dict, v2_line, v3_init_line, mutual_samples,
                        recount2_id_map, cohort_flag
                    )
                else:
                    samp_dict = add_mutual_jx_to_samp_dict(
                        samp_dict, v2_line, v3_line_list, mutual_samples,
                        recount2_id_map, cohort_flag
                    )

                v2_line = next(file2, None)
                v2_coords = coords_from_jx_line(v2_line, v2_coord_locs)
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
        v2_file, v3_file, mutual_samples, recount2_id_map, flag, out_path
    )
    data_file = os.path.join(
        out_path, '{}_sample_jx_dict_{}.json'.format(flag, now)
    )
    with open(data_file, 'w') as output:
        json.dump(sample_jx_dict, output)

    # prepare_boxplot_data(sample_jx_dict, out_path, now, cohort=flag)
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
    parser.add_argument(
        '--sra-mutual-samples',
        help='Provide a json file of SRA samples in recoun3 and recount3.'
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
    sra_both = args.sra_mutual_samples

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'recount_2v3_jx_comparison_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=_LOG_MODE)

    # map recount2 sample IDs to recount3 sample IDs
    recount3_id_map = accession_to_recount3_ids(gtex_ids3, tcga_ids3, sra_ids3)
    recount2_id_map = map_recount_ids_2to3(ids_v2, recount3_id_map)
    
    # compare v2 to v3 jxs for three cohorts
    # execute_jx_comp(
    #     gtex_v2, gtex_v3, recount2_id_map, out_path, now, _GTEX, gtex_both
    # )    
    # execute_jx_comp(tcga_v2, tcga_v3, recount2_id_map, out_path, now, _TCGA)
    execute_jx_comp(
        sra_v2, sra_v3, recount2_id_map, out_path, now, _SRA, sra_both
    )

