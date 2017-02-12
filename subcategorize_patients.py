#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

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

def feature_analysis(labels, feature_matrix, feature_list, max_clus, out_name,
    symptom_line=''):
    '''
    Given two clusters of patients, perform t-test on their feature values
    and sort by p-value. Write out to file. Only write for two clusters.
    Optional argument symptom_line is only used in sequential clustering.
    '''
    # Maps features to p-values of t-tests for the two clusters.
    feature_t_test_dct = {}
    # 'a' denotes the larger cluster.
    cluster_index_lst_a = [i for i, e in enumerate(labels) if e == max_clus]
    # Again, merge the smaller clusters.
    cluster_index_lst_b = [i for i, e in enumerate(labels) if e != max_clus]
    for feature_idx, feature_row in enumerate(feature_matrix.T):
        feature_list_a = feature_row[cluster_index_lst_a]
        feature_list_b = feature_row[cluster_index_lst_b]
        # Compute the unpaired t-test between the two samples of features.
        t_stat, p_value = ttest_ind(feature_list_a, feature_list_b)
        feature = feature_list[feature_idx]
        # Update the dictionary.
        if np.isnan(p_value):
            p_value = 1
        # TODO: Add in whether the features are greater.
        a_mean, b_mean = np.mean(feature_list_a), np.mean(feature_list_b)
        # assert len(feature_list_a) > len(feature_list_b)
        if a_mean == b_mean:
            tag = '='
        elif a_mean > b_mean:
            tag = '>'
        else:
            tag = '<'
        feature_t_test_dct[(feature, tag)] = p_value / 2.0

    feature_t_test_dct = sorted(feature_t_test_dct.items(),
        key=operator.itemgetter(1))

    out = open(out_name, 'w')
    out.write(symptom_line)
    for ((feature, tag), p_value) in feature_t_test_dct:
        out.write('%s\t%g\t%s\n' % (feature, p_value, tag))
    out.close()

def slice_feature_matrix(feature_matrix, master_feature_list, feature_type):
    '''
    Given a feature matrix and a feature type, only keep the columns of the
    matrix that are of that feature type.
    '''
    assert feature_type in ['symptoms', 'treatments']
    if feature_type == 'symptoms':
        feat_list = read_spreadsheet('./data/cancer_other_info_mr_symp.txt')[1]
    else:
        feat_list = read_spreadsheet('./data/cancer_other_info_herbmed.txt')[1]
        feat_list += read_spreadsheet('./data/cancer_drug_2017_sheet2.txt')[1]
    # Get column indices of the current feature type.
    col_idx_lst = [i for i in range(len(master_feature_list)
        ) if master_feature_list[i] in feat_list]
    return feature_matrix[:,col_idx_lst]

def get_cluster_labels(feature_matrix, num_clusters):
    '''
    Clusters using K-Means with the given number of clusters on the cityblock
    distance matrix of the given feature matrix.
    '''
    distance_matrix = squareform(pdist(feature_matrix, metric='cityblock'))
    est = KMeans(n_clusters=num_clusters, n_init=100, random_state=930519)
    est.fit(distance_matrix)
    return list(est.labels_)

def write_clusters(labels, num_clusters, survival_mat, out_name):
    '''
    Given the labels, write the clusters out to file.
    '''
    # First, determine the index of the largest cluster.
    max_clus, max_size = -1, -1
    for i in range(num_clusters):
        curr_size = labels.count(i)
        if curr_size > max_size:
            max_clus, max_size = i, curr_size

    out = open(out_name, 'w')
    out.write('death\ttime\tcluster\n')
    for patient_idx, patient_tup in enumerate(survival_mat):
        label = labels[patient_idx]
        # Merge the smaller clusters; all have index 0.
        group_tag = '0'
        # Largest cluster has index 1.
        if label == max_clus:
            group_tag = '1'
        # TODO: Currently not merging clusters unless we have 3 clusters.
        if num_clusters != 3:
            group_tag = label
        inhospital_id, death, time = patient_tup
        out.write('%d\t%g\t%s\n' % (death, time, group_tag))
    out.close()
    return max_clus

def cluster_full_feature_matrix(isProsnet):
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
    max_clus_id = write_clusters(labels, num_clusters, survival_mat, df_fname)

    # Always use the original feature matrix for feature analysis.
    base_feature_matrix = read_feature_matrix()[0]
    feature_analysis(labels, base_feature_matrix, feature_list, max_clus_id,
        feat_fname)

def characterize_cluster_symptoms(clus_patients, symptom_feature_matrix):
    '''
    Given a list of patient indices, and a symptom feature matrix, find the
    symptoms that best characterize the patients.
    '''
    symptom_cands = []
    symptom_list = read_spreadsheet('./data/cancer_other_info_mr_symp.txt')[1]
    cluster_matrix = symptom_feature_matrix[clus_patients]

    # Get the matrix of non-cluster patients.
    non_clus_patients = [i for i in range(len(symptom_feature_matrix)
        ) if i not in clus_patients]
    non_clus_matrix = symptom_feature_matrix[non_clus_patients]

    # Count the number of non-zero values in each column.
    for symp_idx, symptom in enumerate(symptom_list):
        clus_col = cluster_matrix[:,symp_idx]
        non_clus_col = non_clus_matrix[:,symp_idx]

        t_stat, p_value = ttest_ind(clus_col, non_clus_col)

        if p_value < 0.01:
            if np.mean(clus_col) == np.mean(non_clus_col):
                symptom += '='
            elif np.mean(clus_col) > np.mean(non_clus_col):
                symptom += '>'
            else:
                symptom += '<'
            symptom_cands += [symptom]
    return ', '.join(symptom_cands) + '\n'

def sequential_cluster(isProsnet):
    '''
    First cluster on just symptom features, and then sub-cluster on treatment
    features.
    '''
    if isProsnet:
        suffix = '_' + num_dim
    else:
        suffix = ''

    feature_matrix, feature_list, survival_mat = read_feature_matrix(suffix)

    # First, only cluster on symptoms for sequential clustering.
    symptom_feature_matrix = slice_feature_matrix(feature_matrix, feature_list,
        'symptoms')
    # TODO: Initial number of clusters when clustering on symptoms.
    num_clusters = symptom_feature_matrix.shape[1]

    labels = get_cluster_labels(symptom_feature_matrix, num_clusters)

    # Cluster a second time, this time on treatments.
    # First, only keep herbs/drugs as features.
    feature_matrix = slice_feature_matrix(feature_matrix, feature_list,
        'treatments')
    # Cluster on the trimmed feature matrices.
    base_feature_matrix = read_feature_matrix()[0]
    for i in range(num_clusters):
        # Only cluster again if there are least 20 patients.
        if labels.count(i) < 20:
            continue
        # These are the indices of the patients in the current cluster.
        clus_patients = [j for j, label in enumerate(labels) if label == i]
        # Find the symptoms that best characterize the cluster.
        symptom_line = characterize_cluster_symptoms(clus_patients,
            symptom_feature_matrix)
        sub_labels = get_cluster_labels(feature_matrix[clus_patients], 3)
        sub_survival_mat = [survival_mat[j] for j in clus_patients]
        # Handling different dataframe filenames.
        if isProsnet:
            sub_df_fname = '%s/prosnet_%d_%s.txt' % (df_folder, i, num_dim)
            sub_feat_fname = '%s/prosnet_%d_%s.txt' % (feat_folder, i, num_dim)
        else:
            sub_df_fname = '%s/without_prosnet_%d.txt' % (df_folder, i)
            sub_feat_fname = '%s/without_prosnet_%d.txt' % (feat_folder, i)

        max_sub_clus_id = write_clusters(sub_labels, 3, sub_survival_mat,
            sub_df_fname)
        feature_analysis(sub_labels, base_feature_matrix, feature_list,
            max_sub_clus_id, sub_feat_fname, symptom_line)

def main():
    if len(sys.argv) != 3:
        print 'Usage: python %s matrix_type num_dim' % sys.argv[0]
        exit()
    global matrix_type, num_dim
    matrix_type, num_dim = sys.argv[1:]
    assert matrix_type in ['seq', 'full'] and num_dim.isdigit()

    generate_directories()

    if matrix_type == 'full':
        cluster_full_feature_matrix(True)
        cluster_full_feature_matrix(False)
    else:
        sequential_cluster(True)
        sequential_cluster(False)

    # Call the R script.
    command = 'Rscript survival_model.R %s' % df_folder
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    main()