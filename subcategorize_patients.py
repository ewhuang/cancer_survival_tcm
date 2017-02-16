#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from collections import Counter
from file_operations import read_feature_matrix, read_spreadsheet
import os
import numpy as np
import operator
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
import subprocess
import sys

def generate_directories():
    global df_folder, feat_folder
    df_folder = './data/patient_dataframes_%s' % matrix_type
    feat_folder = './results/feature_p_values_%s' % matrix_type
    plot_folder = './results/survival_plots_%s' % matrix_type
    for folder in [df_folder, feat_folder, plot_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

def get_col_idx_lst(feature_list, feature_type):
    '''
    Given a feature matrix and a feature type, only keep the columns of the
    matrix that are of that feature type.
    '''
    assert feature_type in ['symptoms', 'treatments']

    if feature_type == 'symptoms':
        symp_list = read_spreadsheet('./data/cancer_other_info_mr_symp.txt')[1]
        test_list = read_spreadsheet('./data/incase_check.txt')[1]
        feat_set = set(symp_list).union(test_list)
    elif feature_type == 'treatments':
        herb_list = read_spreadsheet('./data/cancer_other_info_herbmed.txt')[1]
        drug_list = read_spreadsheet('./data/cancer_drug_2017_sheet2.txt')[1]
        feat_set = set(herb_list).union(drug_list)
    # Get column indices of the current feature type.
    col_idx_lst = [i for i, e in enumerate(feature_list) if e in feat_set]
    return col_idx_lst

def get_cluster_labels(feature_matrix, num_clusters):
    '''
    Clusters using K-Means with the given number of clusters on the cityblock
    distance matrix of the given feature matrix.
    '''
    distance_matrix = squareform(pdist(feature_matrix, metric='cityblock'))
    est = KMeans(n_clusters=num_clusters, n_init=1000, random_state=930519)
    est.fit(distance_matrix)
    return list(est.labels_)

def write_clusters(labels, num_clusters, survival_mat, out_name):
    '''
    Given the labels, write the clusters out to file. For situations in which
    we only have two clusters, merge the smaller clusters and label it 0.
    '''
    assert len(labels) == len(survival_mat)
    # First, determine the index of the largest cluster.
    max_clus = Counter(labels).most_common(1)[0][0]

    if num_clusters == 3:
        tag_list = [1 if label == max_clus else 0 for label in labels]
    else:
        print 'Warning: Currently not merging clusters...'
        tag_list = labels[:]

    out = open(out_name, 'w')
    out.write('death\ttime\tcluster\n')
    for patient_idx, patient_tup in enumerate(survival_mat):
        inhospital_id, death, time = patient_tup
        tag = tag_list[patient_idx]
        out.write('%d\t%g\t%d\n' % (death, time, tag))
    out.close()

def get_cluster_symptoms(clus_idx_lst, feature_list, feature_matrix,
    symp_idx_lst):
    '''
    Given a list of patient indices, and a symptom feature matrix, find the
    symptoms that best characterize the patients.
    '''
    symptom_cands = []
    # Get the rows in the feature matrix belonging to the cluster.
    cluster_matrix = feature_matrix[clus_idx_lst]

    # Get the matrix of non-cluster patients.
    non_clus_idx_lst = [i for i in range(len(feature_matrix)
        ) if i not in clus_idx_lst]
    non_clus_matrix = feature_matrix[non_clus_idx_lst]

    assert len(clus_idx_lst) + len(non_clus_idx_lst) == len(feature_matrix)

    for symp_idx in symp_idx_lst:
        clus_col = cluster_matrix[:,symp_idx]
        non_clus_col = non_clus_matrix[:,symp_idx]

        with np.errstate(invalid='ignore'):
            t_stat, p_value = ttest_ind(clus_col, non_clus_col)

        if (p_value / 2.0) < 0.01:
            symptom = feature_list[symp_idx]
            clus_mean, non_clus_mean = np.mean(clus_col), np.mean(non_clus_col)
            if clus_mean == non_clus_mean:
                symptom += '='
            elif clus_mean > non_clus_mean:
                symptom += '>'
            else:
                symptom += '<'
            symptom_cands += [symptom]
    return ', '.join(symptom_cands) + '\n'

def feature_analysis(labels, feature_matrix, feature_list, out_name,
    symp_line=''):
    '''
    Given two clusters of patients, perform t-test on their feature values
    and sort by p-value. Write out to file. Only write for two clusters.
    Optional argument symptom_line is only used in sequential clustering.
    '''
    # Number of columns is equal to number of features.
    assert feature_matrix.shape == (len(labels), len(feature_list))
    # Get the label of the most common cluster.
    max_clus = Counter(labels).most_common(1)[0][0]
    # 'a' denotes the largest cluster.
    max_clus_idx_lst = [i for i, e in enumerate(labels) if e == max_clus]
    # Merge the smaller clusters.
    merged_clus_idx_lst = [i for i, e in enumerate(labels) if e != max_clus]

    # Maps features to p-values of t-tests for the two clusters.
    p_val_dct = {}
    for feature_idx, feature in enumerate(feature_list):
        column = feature_matrix[:,feature_idx]
        # This is the feature values for the larger cluster.
        max_feat_list = column[max_clus_idx_lst]
        assert len(max_feat_list) == len(max_clus_idx_lst)
        # Feature values for merged clusters.
        merged_feat_list = column[merged_clus_idx_lst]
        assert len(max_feat_list) + len(merged_feat_list) == len(labels)
        # Compute the unpaired t-test between the two samples of features.
        with np.errstate(invalid='ignore'):
            t_stat, p_value = ttest_ind(max_feat_list, merged_feat_list)
        if np.isnan(p_value):
            p_value = 1

        # > marker means that largest cluster has a larger mean than the
        # combined clusters.
        max_mean, merge_mean = np.mean(max_feat_list), np.mean(merged_feat_list)
        assert len(max_feat_list) > len(merged_feat_list)
        if max_mean == merge_mean:
            tag = '='
        elif max_mean > merge_mean:
            # TODO. Currently showing the size of the '>' sign.
            tag = '%d>' % len(max_feat_list)
        else:
            tag = '<'
        p_val_dct[(feature, tag)] = p_value / 2.0

    p_val_dct = sorted(p_val_dct.items(), key=operator.itemgetter(1))

    # Write out to file.
    out = open(out_name, 'w')
    out.write(symp_line) # Write out the optional symptom line. Sequential only.
    for ((feature, tag), p_value) in p_val_dct:
        out.write('%s\t%g\t%s\n' % (feature, p_value, tag))
    out.close()

def cluster_full_feature_matrix():
    '''
    Read the feature matrix, cluster, and then perform feature analysis.
    '''
    if isProsnet:
        suffix = '_' + num_dim
        df_fname = '%s/prosnet_%s.txt' % (df_folder, num_dim)
        feat_fname = '%s/prosnet_%s.txt' % (feat_folder, num_dim)
    else:
        suffix = ''
        df_fname = '%s/without_prosnet.txt' % df_folder
        feat_fname = '%s/without_prosnet.txt' % feat_folder

    feature_matrix, feature_list, survival_mat = read_feature_matrix(suffix)
    # TODO: Currently clustering with 3 clusters, merging the smaller 2.
    num_clusters = 3

    labels = get_cluster_labels(feature_matrix, num_clusters)
    write_clusters(labels, num_clusters, survival_mat, df_fname)

    # Always use the original feature matrix for feature analysis.
    base_feature_matrix, base_feat_list, base_surv_mat = read_feature_matrix()
    assert base_surv_mat == survival_mat and feature_list == base_feat_list
    if suffix == '':
        assert np.array_equal(base_feature_matrix, feature_matrix)

    feature_analysis(labels, base_feature_matrix, feature_list, feat_fname)

def sequential_cluster():
    '''
    First cluster on just symptom features, and then sub-cluster on treatment
    features.
    '''
    if isProsnet:
        suffix = '_' + num_dim
    else:
        suffix = ''

    feature_matrix, feature_list, survival_mat = read_feature_matrix(suffix)
    base_feature_matrix, base_feat_lst, base_surv_mat = read_feature_matrix()
    assert base_feat_lst == feature_list and survival_mat == base_surv_mat
    if suffix == '':
        assert np.array_equal(feature_matrix, base_feature_matrix)

    # First, only cluster on symptoms and tests for sequential clustering.
    symp_idx_lst = get_col_idx_lst(feature_list, 'symptoms')
    symp_feature_matrix = feature_matrix[:,symp_idx_lst]
    # TODO: Initial number of clusters when clustering on symptoms.
    # num_clusters = symp_feature_matrix.shape[1]
    num_clusters = 10
    print 'initial num clusters:', num_clusters
    labels = get_cluster_labels(symp_feature_matrix, num_clusters)

    # Cluster a second time, this time on drugs and herbs.
    drug_idx_lst = get_col_idx_lst(feature_list, 'treatments')
    drug_feature_matrix = feature_matrix[:,drug_idx_lst]
    drug_list = [feature_list[i] for i in drug_idx_lst]

    for i in range(num_clusters):
        if labels.count(i) < 20:
            continue
        # These are the indices of the patients in the current cluster.
        clus_idx_lst = [j for j, label in enumerate(labels) if label == i]
        assert len(clus_idx_lst) == labels.count(i)
        # Find the symptoms that best characterize the cluster.
        symptom_line = get_cluster_symptoms(clus_idx_lst, feature_list,
            base_feature_matrix, symp_idx_lst)
        sub_labels = get_cluster_labels(drug_feature_matrix[clus_idx_lst], 3)
        sub_survival_mat = [survival_mat[j] for j in clus_idx_lst]
        # Handling different dataframe filenames.
        if isProsnet:
            sub_df_fname = '%s/prosnet_%d_%s.txt' % (df_folder, i, num_dim)
            sub_feat_fname = '%s/prosnet_%d_%s.txt' % (feat_folder, i, num_dim)
        else:
            sub_df_fname = '%s/without_prosnet_%d.txt' % (df_folder, i)
            sub_feat_fname = '%s/without_prosnet_%d.txt' % (feat_folder, i)

        write_clusters(sub_labels, 3, sub_survival_mat, sub_df_fname)
        # Perform feature analysis on the original feature matrix, but only
        # consider the drugs.
        # TODO: may not be aligning the right features.
        feature_analysis(sub_labels, base_feature_matrix[clus_idx_lst][:,
            drug_idx_lst], drug_list, sub_feat_fname, symptom_line)

def main():
    if len(sys.argv) not in [2, 3]:
        print 'Usage: python %s matrix_type num_dim<optional>' % sys.argv[0]
        exit()
    global matrix_type, isProsnet
    matrix_type, isProsnet = sys.argv[1], False
    assert matrix_type in ['seq', 'full']

    # Optional dimensions argument for ProSNet experiments.
    if len(sys.argv) == 3:
        global num_dim
        num_dim, isProsnet = sys.argv[2], True
        assert num_dim.isdigit()

    generate_directories()

    if matrix_type == 'full':
        cluster_full_feature_matrix()
    else:
        sequential_cluster()

    # Call the R script.
    command = 'Rscript survival_model.R %s' % df_folder
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    main()