#!/usr/bin/python
# -*- coding: utf-8 -*-
# Author: Edward Huang

from file_operations import read_spreadsheet, read_feature_matrix
import numpy as np
import operator
import os
import subprocess

### This script goes through each of the markers, separating patients of each
### cancer subtype into one cluster that contains the marker, and another
### cluster that does not contain the marker. Writes out to file the markers
### that best separate patients into clusters with different survival functions.

def generate_directories():
    folder_lst = ['./results/top_markers', './results/top_marker_plots']
    for results_folder in folder_lst:
        if not os.path.exists(results_folder):
            os.makedirs(results_folder)

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

def write_marker_dataframes():
    '''
    Writes out the dataframes for each of the cancer subtypes.
    '''
    feature_matrix, mast_feature_list, survival_mat = read_feature_matrix('_raw')
    subtype_labels = get_subtype_labels(survival_mat)

    marker_feature_set, test_feature_set = set([]), set([])
    for fname in ('cancer_other_info_mr_symp', 'cancer_syndrome_syndromes',
        'cancer_check_20170324'):
        for feat in read_spreadsheet('./data/%s.txt' % fname)[1]:
            marker_feature_set.add(feat)
            if 'check' in fname: # Create set of medical test features.
                test_feature_set.add(feat)

    marker_patient_dct = [0, {}, {}]
    for cancer_idx in [1, 2]:
    # for cancer_idx in [2]:
        cancer_idx_lst = [i for i, e in enumerate(subtype_labels) if e == cancer_idx]
        # Slice the feature and survival matrices.
        cancer_mat = feature_matrix[cancer_idx_lst]
        surv_mat = [survival_mat[j] for j in cancer_idx_lst]

        for marker in marker_feature_set:
        # for marker in ['阻塞性肺炎']:
            if marker not in mast_feature_list:
                continue
            feat_idx = mast_feature_list.index(marker)
            # Slice the matrix.
            sliced_matrix = cancer_mat[:,feat_idx]
            # Skip features that no patients have (array of all zeros).
            if not np.any(sliced_matrix):
                continue
            # Get the threshold for which to separate patients into two clusters.
            threshold = 0.5 # Default: binary threshold.
            if marker in test_feature_set:
                threshold = np.median(sliced_matrix)
            marker_patient_tup = cluster_and_write(sliced_matrix, surv_mat, cancer_idx, marker, threshold)
            marker_patient_dct[cancer_idx][marker] = marker_patient_tup
    return marker_patient_dct

def cluster_and_write(feature_matrix, survival_mat, cancer_idx, marker, threshold):
    # First element is list of patients with the feature.
    marker_patient_tup = [[], []]
    # Parentheses giving R some trouble.
    marker = marker.replace('(', '')
    marker = marker.replace(')', '')
    # Write out the dataframe.
    fname = './results/top_markers/%s_%s.txt' % (marker, cancer_idx)
    out = open(fname, 'w')
    out.write('death\ttime\tcluster\n')
    for patient_idx, patient_tup in enumerate(survival_mat):
        inhospital_id, death, time = patient_tup
        # Separate the patients into two groups based on whether they have the marker.
        tag = feature_matrix[patient_idx]
        if tag <= threshold:
            tag = 0
        else:
            tag = 1
        marker_patient_tup[tag] += [inhospital_id]
        out.write('%d\t%g\t%d\n' % (death, time, tag))
    out.close()
    # Call the R script to compute survival function differences.
    command = 'Rscript compute_marker_p.R %s_%s' % (marker, cancer_idx)
    subprocess.call(command, shell=True)
    return marker_patient_tup

def main():
    generate_directories()
    marker_patient_dct = write_marker_dataframes()

    sq_dct, non_sq_dct = {}, {}
    folder = './results/top_markers'
    folder_list = os.listdir(folder)[:]
    for fname in folder_list:
        if 'squamous' in fname or '_p.txt' not in fname:
            continue
        f = open('%s/%s' % (folder, fname), 'r')
        try:
            p_val = float(f.readline())
            clus_size = map(int, f.readline().split())
            f.readline() # Skip an empty line because of coxph issues.
            cox_ratio = float(f.readline())
        except Exception:
            continue
        # fname[:-6] is the name of the marker.
        key = (fname[:-8], clus_size[0], clus_size[1], cox_ratio)
        if '_1_' in fname:
            sq_dct[key] = p_val
        else:
            non_sq_dct[key] = p_val
        f.close()
    # Sort the dictionaries.
    sq_dct = sorted(sq_dct.items(), key=operator.itemgetter(1))
    non_sq_dct = sorted(non_sq_dct.items(), key=operator.itemgetter(1))

    for (dct, patient_dct, fname) in [(sq_dct, marker_patient_dct[1],
        'squamous_markers'), (non_sq_dct, marker_patient_dct[2], 'non_squamous_markers')]:
        out = open('./results/top_markers/%s.txt' % fname, 'w')
        out.write('Marker\tSurvival difference p-value\tNon-marker cluster size\tNon-marker patients\tMarker cluster size\tMarker patients\tCox exp(coef)\n')
        for (marker, a_size, b_size, cox_ratio), p_val in dct:
            # Skip lopsided clusterings.
            if a_size < 0.05 * b_size or b_size < 0.05 * a_size:
                continue
            out.write('%s\t%g\t%d\t%s\t%d\t%s\t%g\n' % (marker, p_val, a_size,
                ';'.join(patient_dct[marker][0]), b_size,
                ';'.join(patient_dct[marker][1]), cox_ratio))
        out.close()

if __name__ == '__main__':
    main()