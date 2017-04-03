#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from collections import Counter
from file_operations import read_feature_matrix, read_spreadsheet
from file_operations import read_smoking_history
import itertools
import os
import numpy as np
import operator
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ttest_ind
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize
import subprocess
import sys

def get_subtype_labels(survival_mat):
    # Get the set of squamous and non-squamous patients.
    # TODO: Haven't gotten the small-cell patients.
    non_squamous, squamous = [], []
    f = open('./data/smoking_history.txt', 'r')
    for line in f:
        inhospital_id = line.rstrip('\n\r').split('\t')[0]
        if '腺癌' in line or '大细胞癌' in line or '肺泡' in line:
            non_squamous += [inhospital_id]
        elif '鳞癌' in line:
            squamous += [inhospital_id]
    f.close()
    # Assign labels to patients.
    labels = []
    for patient_idx, patient_tup in enumerate(survival_mat):
        inhospital_id, death, time = patient_tup
        if inhospital_id in non_squamous:
            labels += [1]
        elif inhospital_id in squamous:
            labels += [2]
        else:
            labels += [0]
    return labels

def generate_directories():
    global df_folder, feat_folder
    df_folder = './data/patient_dataframes_%s' % matrix_type
    feat_folder = './results/feature_p_values_%s' % matrix_type
    plot_folder = './results/survival_plots_%s' % matrix_type
    for folder in [df_folder, feat_folder, plot_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

def get_col_idx_lst(feature_list, cluster_features):
    '''
    Given a feature matrix and a feature type, only keep the columns of the
    matrix that are of that feature type.
    '''
    # TODO: First cluster on medical history, then cluster on current status.
    feature_dct = {}
    feature_dct['symptoms'] = read_spreadsheet('./data/cancer_other_info_mr_symp.txt')[1]
    feature_dct['history'] = read_smoking_history()[1]
    feature_dct['herbs'] = read_spreadsheet('./data/cancer_other_info_herbmed.txt')[1]
    feature_dct['drugs'] = read_spreadsheet('./data/cancer_drug_2017_sheet2.txt')[1]
    feature_dct['tests'] = read_spreadsheet('./data/cancer_check_20170324.txt')[1]

    feat_set = set([])
    for feat_type in cluster_features:
        feat_set = feat_set.union(feature_dct[feat_type])
    # Get column indices of the current feature type.
    col_idx_lst = [i for i, e in enumerate(feature_list) if e in feat_set]
    return col_idx_lst

def get_cluster_labels(feature_matrix, num_clusters):
    '''
    Clusters using K-Means with the given number of clusters on the cityblock
    distance matrix of the given feature matrix.
    '''
    # TODO: PCA.
    num_comp = int(feature_matrix.shape[1] * 0.2)
    # num_comp = 50
    pca = PCA(n_components=num_comp)
    distance_matrix = normalize(feature_matrix, norm='l1')
    distance_matrix = pca.fit_transform(distance_matrix)
    distance_matrix = squareform(pdist(distance_matrix, metric))

    est = KMeans(n_clusters=num_clusters, n_init=1000, random_state=930519)
    est.fit(distance_matrix)
    return list(est.labels_)

    # distance_matrix = squareform(pdist(feature_matrix, metric='cityblock'))
    # est = KMeans(n_clusters=num_clusters, n_init=1000, random_state=930519)
    # est.fit(distance_matrix)
    # return list(est.labels_)

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

        if (p_value / 2.0) < 0.1:
            symptom = feature_list[symp_idx]
            clus_mean, non_clus_mean = np.mean(clus_col), np.mean(non_clus_col)
            if clus_mean == non_clus_mean:
                symptom += '='
            elif clus_mean > non_clus_mean:
                symptom += '>'
            elif clus_mean < non_clus_mean:
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
        # assert len(max_feat_list) > len(merged_feat_list)
        if max_mean == merge_mean:
            tag = '='
        elif max_mean > merge_mean:
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

def sequential_cluster(cluster_features):
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
    labels = get_subtype_labels(survival_mat)

    # Cluster a second time, this time on drugs and herbs.
    # TODO: Currently initially clustering on medical tests.
    drug_idx_lst = get_col_idx_lst(feature_list, cluster_features)
    drug_feature_matrix = feature_matrix[:,drug_idx_lst]
    drug_list = [feature_list[i] for i in drug_idx_lst]
    for i in range(3):
        # These are the indices of the patients in the current cluster.
        clus_idx_lst = [j for j, label in enumerate(labels) if label == i]
        assert len(clus_idx_lst) == labels.count(i)

        sub_labels = get_cluster_labels(drug_feature_matrix[clus_idx_lst], 2)
        sub_survival_mat = [survival_mat[j] for j in clus_idx_lst]
        # Handling different dataframe filenames.
        if isProsnet:
            sub_df_fname = '%s/prosnet_%d_%s.txt' % (df_folder, i, num_dim)
            sub_feat_fname = '%s/prosnet_%d_%s.txt' % (feat_folder, i, num_dim)
        else:
            sub_df_fname = '%s/without_prosnet_%d.txt' % (df_folder, i)
            sub_feat_fname = '%s/without_prosnet_%d.txt' % (feat_folder, i)

        write_clusters(sub_labels, 2, sub_survival_mat, sub_df_fname)
        # Perform feature analysis on the original feature matrix.
        feature_analysis(sub_labels, base_feature_matrix[clus_idx_lst][:,
            drug_idx_lst], drug_list, sub_feat_fname)

def main():
    if len(sys.argv) not in [3, 4]:
        print ('Usage: python %s matrix_type metric num_dim<optional>' %
            sys.argv[0])
        exit()
    global matrix_type, metric, isProsnet
    matrix_type, metric, isProsnet = sys.argv[1], sys.argv[2], False
    assert matrix_type in ['seq', 'full']

    # Optional dimensions argument for ProSNet experiments.
    if len(sys.argv) == 4:
        global num_dim
        num_dim, isProsnet = sys.argv[3], True
        assert num_dim.isdigit()

    generate_directories()

    # TODO: Change iterations of features.
    # feature_list = ['tests', 'symptoms', 'herbs', 'drugs', 'history']
    # x = itertools.chain.from_iterable(itertools.combinations(feature_list,
        # r) for r in range(len(feature_list) + 1))
    x = [['tests', 'symptoms', 'history']]

    for i in x:
        if i == ():
            continue
        sequential_cluster(list(i))
        # Call the R script.
        command = 'Rscript survival_model.R %s' % df_folder
        subprocess.call(command, shell=True)

if __name__ == '__main__':
    main()