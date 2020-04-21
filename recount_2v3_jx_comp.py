#!/usr/bin/env python3

"""
recount_2v3_jx_comp.py
Python 3 code for comparing RNA-seq junctions called across mutual samples by
the recount2 vs. recount3 Sequence Read Archive analysis.

"""

import argparse
from datetime import datetime
import logging
import os
import pandas as pd
import os


_LOG_MODE = 'INFO'
_ACCEPTABLE_MOTIFS = ['GTAG', 'GCAG', 'ATAC']


#TODO: first run to collect mutual sample_id list?  Or - parse as we go?
#TODO: switch from "lite" col names to full
#TODO: differential loading of v2/v3?
#TODO: update find_coverages to work with .apply instead of modifying input df
#TODO: improve cycling through samples - leverage pandas or parse as we go?

#user specified inputs needed for precision & recall calculations (prc)
sample_ids = ['18538','28168', '45103', '46189', '7613'] 
#motif to filter junctions by
filter_motif = 'GTAG'
#scaled coverage threshold to filter junctions by 
#cov_thresh = 0
sample_num = '7613'
#names of columns in the file that is read in as a data frame 
dataframe_col_names = ['chromosome', 'start', 'end', 'length', 'strand', 
                       '5prime_motif', '3prime_motif', 'sampleIDs:coverages', 
                       'num_samples', 'total_coverage']
#names of columns that specify junction coordinates & sample ID; used for prc 
col_names_for_prc = ['chromosome', 'start', 'end', 'strand', sample_num]


def find_coverages(sample_ids, cov_col_name, dataframe):
    #processes input data frame, calculates scaled coverage values 
    #for each sample, adds column with sample id and scaled coverages
    #for each junction to the input data frame 
    cov_strings = list(dataframe[cov_col_name]) 
    cov_values = []
    samples = dict.fromkeys(sample_ids)
    cum_cov = 0   #sum of total coverage values across samples
    for sample_id in sample_ids: 
        cov_values = []
        sample_id_string = ',' + sample_id + ':'                    
        for cov_string in cov_strings: 
            #logging.info(sample_string)\n',
            if sample_id_string in cov_string: 
                #isolates sample of interest from other samples found in junction
                first_split = cov_string.split(sample_id_string) 
                #isolates cov value for sample of interest
                second_split = first_split[1].split(',')
                cov_value = int(second_split[0]) 
            else:
                cov_value = 0

            cov_values.append(cov_value)
        total_cov = sum(cov_values)
        #logging.info(total_cov)
        dataframe[sample_id] = [cov_value/total_cov for cov_value in cov_values]
        cum_cov = cum_cov + total_cov
        
        #samples[sample_id] = (cov_values)
        #dataframe[sample_id]
    #value used to scale coverage threshold  
    return cum_cov/(len(sample_ids))
    #adds scaled coverages to dataframe to a new column for the sample 
    #for sample_id in sample_ids: 
       # dataframe[sample_id] = [cov_value/avg_total_cov for cov_value in samples[sample_id]]


def load_data(jx_file):
    # reading (sra data set) from file into jx data frame
    jx_df = pd.read_csv(
        jx_file, sep='\t', header=None, names=dataframe_col_names
    )
    # adding coverages for samples to jx data frame
    avg_total_cov_srav2 = find_coverages(
        sample_ids, 'sampleIDs:coverages', jx_df
    )
    jx_df.loc[:, sample_num] *= avg_total_cov_srav2
    jx_df['motif'] = jx_df['5prime_motif']  + jx_df['3prime_motif']
    return jx_df


def filter_and_compare_jxs(srav2_data, srav3_data):
    coverage_thresholds = [(0, 0), (5, 5), (10, 10)]
    logging.info(coverage_thresholds)
    for cov_thresh in coverage_thresholds:
        # creating new data frames for ground truth & test that are filtered by
        # user-specified thresholds for the motif and scaled coverage values
        srav2_data_filtered = srav2_data[
            (srav2_data['motif'] == filter_motif)
            & (srav2_data[sample_num] != 0)
            & (srav2_data[sample_num] >= cov_thresh[0])
            ]

        srav3_data_filtered = srav3_data[
            (srav3_data['motif'] == filter_motif)
            & (srav3_data[sample_num] != 0)
            & (srav3_data[sample_num] >= cov_thresh[1])
            ]

        # calculating precision and recall

        # number of samples in ground truth
        srav2_num_juncs = srav2_data_filtered.shape[0]

        # number of samples in test
        srav3_num_juncs = srav3_data_filtered.shape[0]
        merged_data = pd.merge(srav2_data_filtered[col_names_for_prc],
                               srav3_data_filtered[col_names_for_prc],
                               on=['chromosome', 'start', 'end', 'strand'],
                               how='inner')

        # number of samples shared by ground truth and test
        num_shared_juncs = merged_data.shape[0]
        logging.info("Coverage threshold: " + str(cov_thresh))
        if num_shared_juncs != 0:
            precision = num_shared_juncs / srav3_num_juncs
            recall = num_shared_juncs / srav2_num_juncs
            logging.info("Precision:" + str(precision))
            logging.info("Recall: " + str(recall))
        else:
            logging.info("No shared junctions.")

    logging.info("srav2 filtered data")
    logging.info(srav2_data_filtered[col_names_for_prc])
    logging.info("srav3 filtered data")
    logging.info(srav3_data_filtered[col_names_for_prc])
    return


def main():
    parser = argparse.ArgumentParser(
        description='Parses recount2 and recount3 junction files and compares'
                    'junctions between samples in both cohorts.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give the path to store junction comparison output and plots.'
    )
    parser.add_argument(
        '--v2-jxs', '-j',
        help='Provide the file containing recount2 junctions.'
    )
    parser.add_argument(
        '--v3-jxs', '-j',
        help='Provide the file containing recount3 junctions.'
    )

    args = parser.parse_args()
    out_path = args.output_path
    recount2 = args.v2_jxs
    recount3 = args.v3_jxs

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    log_file = os.path.join(
        out_path, 'recount_2v3_jx_comparison_log_{}.txt'.format(now)
    )
    logging.basicConfig(filename=log_file, level=_LOG_MODE)

    srav2_data = load_data(recount2)
    srav3_data = load_data(recount3)

    # logging.info ("before filtering\n")
    # logging.info(srav2_data.shape[0])
    # logging.info(srav3_data.shape[0])

    filter_and_compare_jxs(srav2_data, srav3_data)


if __name__ == '__main__':
    main()
