#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from file_operations import read_feature_matrix
import os
import numpy as np
import operator
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
from sklearn.linear_model import LogisticRegression
import subprocess
import sys
import time

def generate_directories():
    global df_folder, feat_folder
    df_folder = './data/patient_dataframes_%s' % clus_meth
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)
    plot_folder = './results/survival_plots_%s' % clus_meth
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)
    feat_folder = './results/feature_p_values_%s' % clus_meth
    if not os.path.exists(feat_folder):
        os.makedirs(feat_folder)

def feature_analysis(labels, feature_matrix, master_feature_list, bad_clusters,
    suffix):
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
        # TODO: Add in whether the features are greater.
        a_mean, b_mean = np.mean(feature_list_a), np.mean(feature_list_b)
        a_size, b_size = len(feature_list_a), len(feature_list_b)
        if a_mean == b_mean:
            tag == 'same'
        elif a_size > b_size:
            if a_mean > b_mean:
                tag = 'larger'
            elif a_mean < b_mean:
                tag = 'smaller'
        elif a_mean > b_mean:
            tag = 'smaller'
        else:
            tag = 'larger'
        feature_t_test_dct[(feature, tag)] = p_value / 2.0

    feature_t_test_dct = sorted(feature_t_test_dct.items(),
        key=operator.itemgetter(1))

    # Determine the output filename.
    if suffix == '':
        out_name = '%s/without_prosnet.txt' % feat_folder
    else:
        out_name = '%s/prosnet_%s.txt' % (feat_folder, num_dim)
    out = open(out_name, 'w')
    for ((feature, tag), p_value) in feature_t_test_dct:
        out.write('%s\t%g\t%s\n' % (feature, p_value, tag))
    out.close()

def regression_analysis(survival_mat, feature_matrix, master_feature_list,
    labels, bad_clusters):
    '''
    Trains a classifier on examples with cluster labels as the training labels.
    Outputs the feature weights.
    '''
    bad_idx_list = [i for i in range(len(labels)) if labels[i] in bad_clusters]
    feature_matrix = [patient for i, patient in enumerate(feature_matrix
        ) if i not in bad_idx_list]
    labels = [label for i, label in enumerate(labels) if i not in bad_idx_list]
    logistic = LogisticRegression()
    logistic.fit(feature_matrix, labels)
    coeffs = logistic.coef_
    print coeffs.shape, len(master_feature_list)

def write_clusters(suffix):
    # Get the feature matrix.
    f = './data/feature_matrices/feature_matrix%s.txt' % suffix
    feature_matrix, master_feature_list, survival_mat = read_feature_matrix(f)

    # TODO: Currently clustering on a distance matrix.
    distance_matrix = squareform(pdist(feature_matrix, metric='cityblock'))
    # TODO: Set the number of clusters.
    num_clusters = 3

    # Perform the clustering.
    if clus_meth == 'kmeans':
        # TODO: Set clustering type and number of clusters.
        est = KMeans(n_clusters=num_clusters, n_init=1000, random_state=930519)
    elif clus_meth == 'spectral':
        est = SpectralClustering(n_clusters=num_clusters, affinity='cosine')
    # TODO: Clustering on distance matrix.
    est.fit(distance_matrix)
    labels = list(est.labels_)

    # Find the degenerate clusters. TODO: Set the minimum cluster size.
    bad_clusters, min_threshold = [], 10
    for i in range(num_clusters):
        if labels.count(i) < min_threshold:
            bad_clusters += [i]

    # TODO: Only write out runs that only have two clusters.
    if num_clusters - len(bad_clusters) == 2:
        # Determine the output filename.
        if suffix == '':
            out_name = '%s/without_prosnet.txt' % df_folder
        else:
            out_name = '%s/prosnet_%s.txt' % (df_folder, num_dim)

        out = open(out_name, 'w')
        out.write('death\ttime\tcluster\n')
        for patient_idx, patient_tup in enumerate(survival_mat):
            label = labels[patient_idx]
            # Don't write clusters with fewer than 10 patients.
            if labels.count(label) < 10:
                continue
            # TODO: Smaller group is 0.
            if labels.count(label) < 50:
                group_tag = '0'
            else:
                group_tag = '1'
            inhospital_id, death, time = patient_tup
            out.write('%d\t%g\t%s\n' % (death, time, group_tag))
        out.close()

        # TODO: Regression analysis of features.
        # regression_analysis(survival_mat, feature_matrix, master_feature_list,
        #     labels, bad_clusters)

        # TODO: Currently using un-imputed feature matrix to analyze features.
        f = './data/feature_matrices/feature_matrix.txt'
        feature_matrix, master_feature_list, survival_mat = read_feature_matrix(f)

        # Perform t-test between the cluster features.
        feature_analysis(labels, feature_matrix, master_feature_list,
            bad_clusters, suffix)
    else:
        print 'Error: results have %d clusters' % (num_clusters - len(
            bad_clusters))

def main():
    if len(sys.argv) != 3:
        print 'Usage: python %s cluster_method num_dim' % sys.argv[0]
        exit()
    global clus_meth, num_dim
    clus_meth, num_dim = sys.argv[1:]
    assert clus_meth in ['kmeans', 'spectral'] and num_dim.isdigit()

    generate_directories()

    # Cluster the patients.
    write_clusters('_%s' % num_dim) # Prosnet feature matrix suffix.
    write_clusters('') # Prosnet-less suffix.
    # # Call the R script.
    command = 'Rscript survival_model.R %s' % df_folder
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    # start_time = time.time()
    main()
    # print("--- %s seconds ---" % (time.time() - start_time))