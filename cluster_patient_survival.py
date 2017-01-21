#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from file_operations import read_feature_matrix
import numpy as np
import operator
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
import subprocess
import sys
import time

# This script takes the feature matrix, and, depending on the argument, uses
# a set of features to cluster the patients, then calls the R script to plot
# the survival models.

def feature_analysis(labels, feature_matrix, master_feature_list, bad_clusters):
    '''
    Given two clusters of patients, perform t-test on their feature values
    and sort by p-value. Write out to file. Only write for two clusters.
    '''
    # Maps features to p-values of t-tests for the two clusters.
    feature_t_test_dct = {}
    # Loop through the two possible clusters.
    cluster_index_lst = []
    for cluster in list(set(labels).difference(bad_clusters)):
        # First element are the patient indices in cluster 1, same for #2.
        cluster_index_lst += [np.nonzero(labels == cluster)]
    for feature_idx, feature_row in enumerate(feature_matrix.T):
        feature_list_a = feature_row[cluster_index_lst[0]]
        feature_list_b = feature_row[cluster_index_lst[1]]
        # Compute the unpaired t-test between the two samples of features.
        t_stat, p_value = ttest_ind(feature_list_a, feature_list_b)
        feature = master_feature_list[feature_idx]
        # Update the dictionary.
        if np.isnan(p_value):
            p_value = 1
        feature_t_test_dct[feature] = p_value

    feature_t_test_dct = sorted(feature_t_test_dct.items(),
        key=operator.itemgetter(1))

    out = open('./results/feature_analyses/%s_%d_features.txt' % (
        cluster_method, num_clusters), 'w')
    for (feature, p_value) in feature_t_test_dct:
        out.write('%s\t%g\n' % (feature, p_value))
    out.close()

def cluster_and_write(feature_matrix, survival_dct, master_feature_list):
    '''
    Perform clustering and write out the dataframe to file.
    '''
    if feature_type == 'multiple':
        # herb_list = ['咳痰', '白术', '阴影部位', '食欲不振']
        herb_list = ['食欲不振', '阴影部位', '其它费用', '乳酸脱氢酶(血清)']
        herb_indices = map(master_feature_list.index, herb_list)

    for feature_idx, feature in enumerate(master_feature_list):
        # This block puts together the known herbs and the current index.
        if feature_type == 'multiple':
            if feature_idx in herb_indices:
                continue
            selected_indices = [feature_idx] + herb_indices
        elif feature_type == 'single':
            selected_indices = [feature_idx]

        new_feature_matrix = feature_matrix[:,selected_indices]
        if cluster_method == 'kmeans':
            classifier = KMeans(n_clusters=num_clusters, random_state=9305
                ).fit(new_feature_matrix)
        labels = list(classifier.labels_)

        # Find the degenerate clusters.
        bad_clusters, min_threshold = [], 0.2 * len(feature_matrix)
        for i in range(num_clusters):
            if labels.count(i) < min_threshold:
                bad_clusters += [i]

        if len(set(labels)) - len(bad_clusters) <= 1:
            continue

        out = open('./data/patient_dataframes_%s_features/%s_%d_%s.txt' % (
            feature_type, cluster_method, num_clusters, feature), 'w')
        out.write('death\ttime\tcluster\n')
        for patient_idx, inhospital_id in enumerate(survival_dct.keys()):
            death, time = survival_dct[inhospital_id]
            cluster = labels[patient_idx]
            if cluster in bad_clusters:
                continue
            out.write('%d\t%g\t%d\n' % (death, time, labels[patient_idx]))
        out.close()

    # If we have two clusters, then perform t-test between the cluster features.
    # if num_clusters - len(bad_clusters) == 2:
    #     feature_analysis(labels, feature_matrix, master_feature_list,
    #         bad_clusters)

def main():
    if len(sys.argv) != 4:
        print ('Usage: python %s kmeans num_clusters single/multiple/num_dim'
            ) % sys.argv[0]
        exit()
    global cluster_method, num_clusters, feature_type
    cluster_method, num_clusters = sys.argv[1], int(sys.argv[2])
    feature_type = sys.argv[3]
    assert feature_type in ['single', 'multiple'] or feature_type.isdigit()

    # If feature_type is a digit, then we are running on Prosnet-enriched data.
    if feature_type.isdigit():
        f = './data/feature_matrices/feature_matrix_%s.txt' % feature_type
    else:
        f = './data/feature_matrices/feature_matrix.txt'
    feature_matrix, master_feature_list, survival_dct = read_feature_matrix(f)
    # Cluster the patients.
    cluster_and_write(feature_matrix, survival_dct, master_feature_list)
    # Call the R script.
    subprocess.call('Rscript survival_model.R %s' % feature_type, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))