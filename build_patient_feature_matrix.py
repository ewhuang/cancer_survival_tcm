#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from file_operations import read_spreadsheet
import numpy as np
import time

# This script creates the feature matrices to get ready for experiments.

def build_feature_matrix(feature_dct_list, master_feature_list, patient_list):
    '''
    Takes the feature list and a dictionary, and build a feature matrix.
    '''
    feature_matrix = []
    # Loop through the global patient list, which is the intersection of all
    # feature matrix keys.
    for inhospital_id in patient_list:
        # Get the row for the patient.
        row = [0 for i in master_feature_list]
        # Get the values from each of the feature dictionaries.
        for feature_dct in feature_dct_list:
            tuple_list = feature_dct[inhospital_id]
            for (feature, feature_freq) in tuple_list:
                row[master_feature_list.index(feature)] += feature_freq
        feature_matrix += [row]

    feature_matrix = np.array(feature_matrix)
    # Remove columns that have fewer than 10 non-zero values.
    good_features = np.apply_along_axis(np.count_nonzero, 0, feature_matrix
        ) >= 10
    feature_matrix = feature_matrix[:, good_features]
    master_feature_list = list(np.array([master_feature_list])[:,good_features]
        [0])
    return feature_matrix, master_feature_list

def write_feature_matrix(feature_matrix, master_feature_list, patient_list,
    survival_dct):
    '''
    Writes the feature matrix out to file, along with column labels. First
    column should be inhospital_id's, and the second/third column should be the
    death/time event.
    '''
    out = open('./data/feature_matrix.txt', 'w')
    out.write('\tdeath\ttime\t%s\n' % '\t'.join(master_feature_list))
    for i, row in enumerate(feature_matrix):
        inhospital_id = patient_list[i]
        death, time = survival_dct[inhospital_id][0]
        out.write('%s\t%d\t%f\t%s\n' % (inhospital_id, death, time, '\t'.join(
            map(str, row))))
    out.close()

def main():
    survival_dct, dummy_return = read_spreadsheet('./data/cancer_life_days.txt')
    # Initialize our list of inhospital_id's.
    patient_list = set(survival_dct.keys())

    feature_dct_list, master_feature_list = [], []
    for fname in ('cancer_other_info_herbmed', 'cancer_other_info_mr_symp',
        'cancer_syndrome_syndromes', 'incase_check', 'cancer_drug_2017_sheet2'):
        feature_dct, feature_list = read_spreadsheet('./data/%s.txt' % fname)

        # Update the list of feature dictionaries.
        feature_dct_list += [feature_dct]
        master_feature_list += feature_list
        # Update the list of inhospital_id's.
        patient_list = patient_list.intersection(feature_dct.keys())

    # Structure the patient list for indexing when writing out to file.
    patient_list = list(patient_list)
    
    # Create the numpy array, and remove bad columns.
    feature_matrix, master_feature_list = build_feature_matrix(feature_dct_list,
        master_feature_list, patient_list)

    # Write out to file.
    write_feature_matrix(feature_matrix, master_feature_list, patient_list,
        survival_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))