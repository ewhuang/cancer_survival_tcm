#!/usr/bin/python
# -*- coding: utf-8 -*-
### Author: Edward Huang

from file_operations import read_spreadsheet
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
import sys
import time

# This script creates the feature matrices to get ready for experiments.

def build_feature_matrix(feature_dct_list, master_feature_list, patient_list):
    '''
    Takes the feature list and a dictionary, and build a feature matrix.
    '''
    feature_matrix = []

    for inhospital_id in patient_list:
        # Get the row for the patient.
        row = [0 for i in master_feature_list]
        # Get the values from each of the feature dictionaries.
        for feature_dct in feature_dct_list:
            # An inhospital ID might not be in feature_dct with prosnet.
            if inhospital_id not in feature_dct:
                continue
            for (feature, feature_freq) in feature_dct[inhospital_id]:
                # TODO: Decide whether or not to keep the numerical features.
                # row[master_feature_list.index(feature)] += feature_freq
                row[master_feature_list.index(feature)] = 1
        feature_matrix += [row]

    return np.array(feature_matrix)

def impute_missing_data(feature_matrix, master_feature_list):
    '''
    Given the feature matrix and the column labels (master_feature_list), impute
    the missing feature data by getting the Prosnet vectors.
    '''
    def read_prosnet_output(master_feature_list):
        '''
        Reads the output low-dimensional vectors created by prosnet.
        '''
        vector_dct = {}
        f = open('./data/prosnet_data/prosnet_node_vectors_%s_dims.vec' %
            num_dim, 'r')
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.split()
            feature, vector = line[0], map(float, line[1:])
            assert len(vector) == int(num_dim) and feature not in vector_dct
            vector_dct[feature] = vector
        f.close()
        # Reorganize the matrix according to the order of master_feature_list.
        vector_matrix = []
        for feature in master_feature_list:
            vector_matrix += [vector_dct[feature]]
        return np.array(vector_matrix)

    vector_matrix = read_prosnet_output(master_feature_list)
    similarity_matrix = cosine_similarity(vector_matrix)
    # Multiply the two feature matrix and the similarity matrix.
    enriched_feature_matrix = np.dot(feature_matrix, similarity_matrix)

    for row_idx, row in enumerate(enriched_feature_matrix):
        # Set threshold of a matrix equal to one standard deviation above mean.
        # TODO: Tune this threshold.
        # threshold = np.mean(row) + np.std(row)
        threshold = np.mean(row)
        row[row < threshold] = 0
        row[row > threshold] = 1
        enriched_feature_matrix[row_idx] = row

    return enriched_feature_matrix

def write_feature_matrix(feature_matrix, master_feature_list, patient_list,
    survival_dct):
    '''
    Writes the feature matrix out to file, along with column labels. First
    column should be inhospital_id's, and the second/third column should be the
    death/time event.
    '''
    results_folder = './data/feature_matrices'
    if isImputation:
        out_fname = '%s/feature_matrix_%s.txt' % (results_folder, num_dim)
    else:
        out_fname = '%s/feature_matrix.txt' % results_folder
    out = open(out_fname, 'w')
    out.write('patient_id\tdeath\ttime\t%s\n' % '\t'.join(master_feature_list))
    for i, row in enumerate(feature_matrix):
        inhospital_id = patient_list[i]
        death, time = survival_dct[inhospital_id][0] # one event per patient.
        out.write('%s\t%d\t%f\t%s\n' % (inhospital_id, death, time, '\t'.join(
            map(str, row))))
    out.close()

def main():
    if len(sys.argv) not in [1, 2]:
        print 'Usage:python %s num_dim<optional>' % sys.argv[0]
        exit()
    global isImputation
    isImputation = False
    if len(sys.argv) == 2:
        # Prosnet has a command line argument.
        global num_dim
        num_dim, isImputation = sys.argv[1], True
        assert num_dim.isdigit()

    survival_dct = read_spreadsheet('./data/cancer_life_days.txt')[0]
    # Initialize our list of inhospital_id's.
    patient_list = list(survival_dct.keys())

    feature_dct_list, master_feature_list = [], []
    for fname in ('cancer_other_info_herbmed', 'cancer_other_info_mr_symp',
        'cancer_syndrome_syndromes', 'incase_check', 'cancer_drug_2017_sheet2'):
        feature_dct, feature_list = read_spreadsheet('./data/%s.txt' % fname)
        # Update the list of feature dictionaries.
        feature_dct_list += [feature_dct]
        master_feature_list += feature_list

    # Create the numpy array, and remove bad columns.
    feature_matrix = build_feature_matrix(feature_dct_list, master_feature_list,
        patient_list)

    if isImputation:
        feature_matrix = np.add(feature_matrix, impute_missing_data(
            feature_matrix, master_feature_list))

    # Write out to file.
    write_feature_matrix(feature_matrix, master_feature_list, patient_list,
        survival_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))