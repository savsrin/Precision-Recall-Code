import pandas as pd 
import os 
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


def find_coverages(sample_ids,cov_col_name, dataframe):
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
            #print(sample_string)\n',
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
        #print(total_cov)
        dataframe[sample_id] = [cov_value/total_cov for cov_value in cov_values]
        cum_cov = cum_cov + total_cov
        
        #samples[sample_id] = (cov_values)
        #dataframe[sample_id]
    #value used to scale coverage threshold  
    return cum_cov/(len(sample_ids))
    #adds scaled coverages to dataframe to a new column for the sample 
    #for sample_id in sample_ids: 
       # dataframe[sample_id] = [cov_value/avg_total_cov for cov_value in samples[sample_id]]
def main(): 

    #reading ground truth data (srav2 data set) from file into data frame
    srav2_data = pd.read_csv('srav2.junctions.test_set.tsv',   
                      sep='\t',
                      header=None, 
                      names=dataframe_col_names)
    #adding coverages for samples to ground truth data frame
    avg_total_cov_srav2 = find_coverages(sample_ids, 'sampleIDs:coverages', srav2_data)
    srav2_data.loc[:, sample_num] *=  avg_total_cov_srav2
     

    #reading test data (srav2 data set) from file into data drame
    srav3_data = pd.read_csv('srav3.junctions.test_set.tsv',   #reading test data
                      sep='\t',
                      header=None, 
                      names=dataframe_col_names)

    #adding coverages for samples to test data frame
    avg_total_cov_srav3 = find_coverages(sample_ids, 'sampleIDs:coverages', srav3_data) 
    srav3_data.loc[:,sample_num] *=  avg_total_cov_srav3

    #adding column with the full motif for each junction to dataframes
    srav2_data['motif'] = srav2_data['5prime_motif']  + srav2_data['3prime_motif'] 
    srav3_data['motif'] = srav3_data['5prime_motif']  + srav3_data['3prime_motif']
    
   # print ("before filtering\n")
   # print(srav2_data.shape[0])
   # print(srav3_data.shape[0]) 
   
    
    
    coverage_thresholds = [(0,0), (5,5), 
                           (10,10)]
    print(coverage_thresholds)
    for cov_thresh in coverage_thresholds: 

        #creating new data frames for ground truth & test that are filtered by user
        #specified thresholds for the motif and scaled coverage values 
        srav2_data_filtered = srav2_data[(srav2_data['motif'] == filter_motif) 
                                        & (srav2_data[sample_num] != 0)
                                        & (srav2_data[sample_num] >= cov_thresh[0])]
       
        srav3_data_filtered = srav3_data[(srav3_data['motif'] == filter_motif) 
                                         & (srav3_data[sample_num] != 0)
                                         & (srav3_data[sample_num] >= cov_thresh[1])]

        #calculating precision and recall 

        #number of samples in ground truth
        srav2_num_juncs = srav2_data_filtered.shape[0]
        
        #number of samples in test 
        srav3_num_juncs = srav3_data_filtered.shape[0] 
        merged_data = pd.merge(srav2_data_filtered[col_names_for_prc], 
                       srav3_data_filtered[col_names_for_prc], 
                       on=['chromosome', 'start', 'end', 'strand'],
                       how = 'inner')

        #number of samples shared by ground truth and test 
        num_shared_juncs = merged_data.shape[0]
        print ("Coverage threshold: " + str(cov_thresh))
        if num_shared_juncs != 0: 
            precision = num_shared_juncs/srav3_num_juncs
            recall = num_shared_juncs/srav2_num_juncs
            print("Precision:" + str(precision))
            print("Recall: " + str(recall))
        else: 
            print ("No shared junctions.")
    print("srav2 filtered data")
    print(srav2_data_filtered[col_names_for_prc])
    print("srav3 filtered data")
    print(srav3_data_filtered[col_names_for_prc])

if __name__ == '__main__':
    main()
