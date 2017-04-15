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

def generate_directories():
    global df_folder, feat_folder
    df_folder = './data/patient_dataframes_seq'
    feat_folder = './results/feature_p_values_seq'
    plot_folder = './results/survival_plots_seq'
    for folder in [df_folder, feat_folder, plot_folder]:
        if not os.path.exists(folder):
            os.makedirs(folder)

def get_subtype_labels(survival_mat):
    '''
    Split the patients into squamous-cell lung carcinoma and non-SQ NSCLC
    patients.
    '''
    # Get the set of squamous and non-squamous patients.
    squamous, non_squamous = [], []
    f = open('./data/smoking_history.txt', 'r')
    for line in f:
        inhospital_id = line.rstrip('\n\r').split('\t')[0]
        cancer = line.rstrip('\n\r').split('\t')[9]
        if ('腺癌' in cancer or '大细胞癌' in cancer or '肺泡' in cancer or
            '乳头状癌' in cancer):
            non_squamous += [inhospital_id]
        elif '鳞癌' in cancer:
            squamous += [inhospital_id]
    f.close()
    # Assign labels to patients.
    labels = []
    for (inhospital_id, death, time) in survival_mat:
        if inhospital_id in squamous:
            labels += [1]
        elif inhospital_id in non_squamous:
            labels += [2]
        else:
            labels += [0]
    return labels

def get_col_idx_lst(feature_list, feat_comb):
    '''
    Given a feature list and list of feature types, return the list of indices
    of the feature list that match the given feature types (feat_comb).
    '''
    # Map each feature type to its filename.
    fname_dct = ({'symptoms':'cancer_other_info_mr_symp', 'herbs':
        'cancer_other_info_herbmed', 'drugs':'cancer_drug_2017_sheet2',
        'tests':'cancer_check_20170324'})
    feat_set = set([])
    for feat_type in feat_comb:
        # Get the current feature set.
        if feat_type == 'history':
            curr_features = read_smoking_history()[1]
        else:
            curr_features = read_spreadsheet('./data/%s.txt' %
                fname_dct[feat_type])[1]
        # Update the overall feature set.
        feat_set = feat_set.union(curr_features)
    # Get the indices of all input feature types.
    col_idx_lst = [i for i, e in enumerate(feature_list) if e in feat_set]
    return col_idx_lst

def get_cluster_labels(feature_matrix):
    '''
    Clusters using K-Means with 2 clusters on the PCA, then distance matrix
    with the input metric.
    '''
    # l1-normalize.
    norm_matrix = normalize(feature_matrix, norm='l1')
    # Perform PCA.
    num_comp = int(feature_matrix.shape[1] * 0.2)
    pca = PCA(n_components=num_comp)
    pca_matrix = pca.fit_transform(norm_matrix)
    # Compute distance matrix.
    distance_matrix = squareform(pdist(pca_matrix, metric))
    # Always cluster with 2 clusters.
    est = KMeans(n_clusters=2, n_init=1000, random_state=930519)
    est.fit(distance_matrix)
    return list(est.labels_)

def get_vkps_labels(base_feature_matrix):
    '''
    For VKPS, one cluster of patients with score greater than 60, rest into
    the other cluster. Return a set of labels. Uses the base feature matrix.
    '''
    labels = [1 if kps > 60 else 0 for kps in base_feature_matrix]
    return labels

def write_clusters(labels, survival_mat, out_name):
    '''
    Given the labels, write the clusters out to file. For situations in which
    we only have two clusters, merge the smaller clusters and label it 0.
    '''
    assert len(labels) == len(survival_mat)
    # First, determine the index of the larger cluster.
    max_clus = Counter(labels).most_common(1)[0][0]
    # Larger cluster gets label 1.
    tag_list = [1 if label == max_clus else 0 for label in labels]

    out = open(out_name, 'w')
    out.write('death\ttime\tcluster\n')
    for patient_idx, patient_tup in enumerate(survival_mat):
        inhospital_id, death, time = patient_tup
        tag = tag_list[patient_idx]
        out.write('%d\t%g\t%d\n' % (death, time, tag))
    out.close()

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
            # Larger cluster gets label 1.
            tag = '%d>' % 1
        else:
            tag = '<'
        p_val_dct[(feature, tag, len(max_feat_list), max_mean,
            len(merged_feat_list), merge_mean)] = p_value / 2.0

    p_val_dct = sorted(p_val_dct.items(), key=operator.itemgetter(1))

    # Write out to file.
    out = open(out_name, 'w')
    out.write(symp_line) # Write out the optional symptom line. Sequential only.
    for ((feature, tag, max_len, max_mean, merge_len, merge_mean),
        p_value) in p_val_dct:
        out.write('%s\t%g\t%d\t%g\t%d\t%g\t%s\n' % (feature, p_value, max_len,
        max_mean, merge_len, merge_mean, tag))
    out.close()

def sequential_cluster(feat_comb):
    '''
    First split data by cancer subtype, and then sub-cluster on treatment
    features.
    '''
    if isImpute in ['prosnet', 'mean']:
        suffix = '_' + opt_arg
    else:
        suffix = ''

    feature_matrix, feature_list, survival_mat = read_feature_matrix(suffix)
    base_feature_matrix, base_feat_lst, base_surv_mat = read_feature_matrix(
        'unnormalized')
    assert base_feat_lst == feature_list and survival_mat == base_surv_mat

    # First, only cluster on symptoms and tests for sequential clustering.
    subtype_labels = get_subtype_labels(survival_mat)

    # Get the sliced feature matrix, depending on the feat_comb.
    if isImpute == 'vkps': # Use VKPS as the only feature.
        # Use the base feature matrix for VKPS, since we're thresholding at 60.
        sub_feature_matrix = base_feature_matrix[:,feature_list.index('VKPS')]
    else:
        feat_idx_lst = get_col_idx_lst(feature_list, feat_comb)
        sub_feature_matrix = feature_matrix[:,feat_idx_lst]
    # Skip the 0th subtype, since it's in the 'other' category.
    for i in [1, 2]:
        # These are the indices of the patients with the current subtype.
        clus_idx_lst = [j for j, label in enumerate(subtype_labels
            ) if label == i]
        assert len(clus_idx_lst) == subtype_labels.count(i)

        clus_feat_matrix = sub_feature_matrix[clus_idx_lst]
        if isImpute == 'vkps':
            sub_labels = get_vkps_labels(clus_feat_matrix)
        else:
            sub_labels = get_cluster_labels(clus_feat_matrix)
        sub_survival_mat = [survival_mat[j] for j in clus_idx_lst]
        # Handling different dataframe filenames.
        if isImpute == 'prosnet':
            sub_df_fname = '%s/prosnet_%d_%s.txt' % (df_folder, i, opt_arg)
            sub_feat_fname = '%s/prosnet_%d_%s.txt' % (feat_folder, i, opt_arg)
        elif isImpute in ['mean', 'vkps']:
            sub_df_fname = '%s/%s_%d.txt' % (df_folder, isImpute, i)
            sub_feat_fname = '%s/%s_%d.txt' % (feat_folder, isImpute, i)
        else:
            sub_df_fname = '%s/without_prosnet_%d.txt' % (df_folder, i)
            sub_feat_fname = '%s/without_prosnet_%d.txt' % (feat_folder, i)

        write_clusters(sub_labels, sub_survival_mat, sub_df_fname)
        # Perform feature analysis on all features.
        feature_analysis(sub_labels, base_feature_matrix[clus_idx_lst],
            feature_list, sub_feat_fname)

        # Call the R script.
        command = 'Rscript plot_kaplan_meiers.R %s' % sub_df_fname
        subprocess.call(command, shell=True)

def main():
    if len(sys.argv) not in [2, 3, 4]:
        print ('Usage: python %s metric num_dim/mean/vkps<optional> '
            'partial<optional>' % sys.argv[0])
        exit()
    global metric, isImpute
    metric, isImpute, isPartial = sys.argv[1], False, False

    # No imputation with partial features.
    if 'partial' in sys.argv and len(sys.argv) == 3:
        isPartial = True
    # Optional dimensions argument for ProSNet experiments or mean imputation.
    elif len(sys.argv) >= 3:
        global opt_arg # This argument is either num_dim, 'mean', or 'vkps'.
        opt_arg = sys.argv[2]
        assert opt_arg.isdigit() or opt_arg in ['mean', 'vkps']
        if opt_arg.isdigit():
            isImpute = 'prosnet'
        else:
            isImpute = opt_arg

    if 'partial' in sys.argv:
        isPartial = True

    generate_directories()

    if isImpute == 'vkps':
        feat_comb_list = ['VKPS']
    else:
        # If not just VKPS, iterate through all combinations of features.
        feature_list = ['tests', 'symptoms', 'herbs', 'drugs', 'history']
        feat_comb_list = itertools.chain.from_iterable(itertools.combinations(
            feature_list, r) for r in range(len(feature_list) + 1))
    if isPartial:
        feat_comb_list = [['symptoms', 'history']]

    for feat_comb in feat_comb_list:
        if feat_comb == ():
            continue
        sequential_cluster(list(feat_comb))

if __name__ == '__main__':
    main()