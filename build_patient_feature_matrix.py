#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

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

def main():
    survival_dct, dummy_return = read_spreadsheet('./data/cancer_life_days.txt')
    herb_dct, unique_herb_list = read_spreadsheet('./data/cancer_other_info_'
        'herbmed.txt')
    symptom_dct, unique_symptom_list = read_spreadsheet('./data/cancer_other_'
        'info_mr_symp.txt')
    syndrome_dct, unique_syndrome_list = read_spreadsheet('./data/cancer_'
        'syndrome_syndromes.txt')
    incase_dct, unique_incase_list = read_spreadsheet('./data/incase_check.txt')
    drug_dct, unique_drug_list = read_spreadsheet('./data/cancer_drug_2017_'
        'sheet2.txt')
    # herb_dct has more patients than survival_dct. Only use survival_dct keys.
    # print set(herb_dct.keys()).difference(survival_dct.keys())

    write_feature_list(unique_herb_list + unique_syndrome_list + 
        unique_syndrome_list + unique_incase_list + unique_drug_list)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))