#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

import numpy as np
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

def main():
    survival_dct, dummy_return = read_spreadsheet('./data/cancer_life_days.txt')
    global patient_list
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
    np.savetxt('./data/feature_matrix.txt', feature_matrix)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))