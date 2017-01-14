#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

import numpy as np
from sklearn.cluster import KMeans
import sys
import time

# This script creates the feature matrices to get ready for clustering.

def read_spreadsheet(fname):
    '''
    Depending on the spreadsheet, return a dictionary mapping the inhospital_id
    to the feature list. All feature lists should be of the same length; some
    are binary (like the survival labels) and some are frequency vectors.
    Key: inhospital_id -> str
    Value: varies from spreadsheet to spreadsheet.
    '''
    feature_dct, unique_feature_list = {}, set([])
    f = open(fname, 'r')
    for i, line in enumerate(f):
        # Skip all header lines.
        if i == 0:
            continue
        line = line.strip().split('\t')
        if len(line) == 1:
            continue
        # Survival events. Spit out binary labels, along with the event length.
        if 'life_days' in fname:
            inhospital_id, feature, feature_freq  = line[0], line[3], line[6]
            if '死亡' in feature:
                feature = 1
            else:
                feature = 0
        # Different feature and frequency cases for different spreadsheets.
        elif 'herbmed' in fname:
            assert len(line) == 4
            inhospital_id, feature_freq, feature = line[:3]
        elif 'mr_symp' in fname:
            assert len(line) == 3
            inhospital_id, feature, feature_freq = line
        elif 'syndrome_syndromes' in fname:
            assert len(line) == 3
            inhospital_id, feature, feature_freq = line[1], line[2], 1
        elif 'incase_check' in fname:
            assert len(line) == 9
            inhospital_id, feature, feature_freq = line[0], line[3], line[4]
            if feature_freq == '无':
                feature_freq = 0
            elif not feature_freq.isdigit():
                feature_freq = 1
        elif 'drug_2017' in fname:
            assert len(line) == 10
            inhospital_id, feature, feature_freq = line[0], line[1], line[4]

        if inhospital_id not in feature_dct:
            feature_dct[inhospital_id] = []
        feature_dct[inhospital_id] += [(feature, float(feature_freq))]
        unique_feature_list.add(feature)
    f.close()
    return feature_dct, list(unique_feature_list)

def write_feature_list(feature_list):
    '''
    Record the order of feature lists for our matrix.
    '''
    out = open('./data/feature_list.txt', 'w')
    out.write('\n'.join(feature_list))
    out.close()

def build_feature_matrix(feature_dct_list, master_feature_list):
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
    return np.array(feature_matrix)

def cluster_and_write(feature_matrix, survival_dct):
    '''
    Perform clustering and write out the dataframe to file.
    '''
    classifier = KMeans(n_clusters=num_clusters, random_state=9305).fit(
        feature_matrix)
    labels = list(classifier.labels_)

    # Find the degenerate clusters.
    bad_clusters, min_threshold = [], 0.2 * len(feature_matrix)
    for i in range(num_clusters):
        if labels.count(i) < min_threshold:
            bad_clusters += [i]

    out = open('./data/patient_dataframes/%s_%d_df.txt' % (cluster_method,
        num_clusters), 'w')
    out.write('death\ttime\tcluster\n')
    for patient_idx, inhospital_id in enumerate(patient_list):
        death, time = survival_dct[inhospital_id][0]
        cluster = labels[patient_idx]
        if cluster in bad_clusters:
            continue
        out.write('%d\t%g\t%d\n' % (death, time, labels[patient_idx]))
    out.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage: python %s kmeans num_clusters' % sys.argv[0]
        exit()
    global patient_list, cluster_method, num_clusters
    cluster_method, num_clusters = sys.argv[1], int(sys.argv[2])

    survival_dct, dummy_return = read_spreadsheet('./data/cancer_life_days.txt')
    patient_list = set(survival_dct.keys())
    feature_dct_list, master_feature_list = [], []
    for fname in ('cancer_other_info_herbmed', 'cancer_other_info_mr_symp',
        'cancer_syndrome_syndromes', 'incase_check', 'cancer_drug_2017_sheet2'):
        feature_dct, feature_list = read_spreadsheet('./data/%s.txt' % fname)

        feature_dct_list += [feature_dct]
        master_feature_list += feature_list
        patient_list = patient_list.intersection(feature_dct.keys())

    write_feature_list(master_feature_list)
    feature_matrix = build_feature_matrix(feature_dct_list, master_feature_list)
    
    # Removing columns that are mostly zero.
    print feature_matrix.shape
    feature_matrix = feature_matrix[:, np.apply_along_axis(np.count_nonzero, 0,
        feature_matrix) >= 5]
    print feature_matrix.shape

    np.savetxt('./data/feature_matrix.txt', feature_matrix)

    # Clustering
    cluster_and_write(feature_matrix, survival_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))