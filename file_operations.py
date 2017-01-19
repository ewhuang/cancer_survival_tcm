#!/usr/bin/python
# -*- coding: utf-8 -*-

### These functions are shared in a variety of scripts. Most are typically
### functions that read files.

from collections import OrderedDict
import numpy as np

# build_patient_feature_matrix.py
# drug_herb_synergistic_survival.py
# create_prosnet_input.py
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
            # Convert to months.
            feature_freq = float(feature_freq) / 30.0
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
            try:
                feature_freq = float(feature_freq)
            except:
                feature_freq = 1
        elif 'drug_2017' in fname:
            assert len(line) == 10
            inhospital_id, feature, feature_freq = line[0], line[1], line[4]

        if inhospital_id not in feature_dct:
            feature_dct[inhospital_id] = []
        # Deal with folder path / issues.
        if type(feature) == str and '/' in feature:
            feature = feature.replace('/', '_')
        feature_dct[inhospital_id] += [(feature, float(feature_freq))]
        unique_feature_list.add(feature)
    f.close()
    return feature_dct, list(unique_feature_list)

# cluster_patient_survival.py
def read_feature_matrix():
    '''
    Reads the feature matrix of the patient data.
    '''
    feature_matrix, master_feature_list, survival_dct = [], [], OrderedDict({})
    f = open('./data/feature_matrix.txt', 'r')
    for i, line in enumerate(f):
        line = line.strip().split('\t')
        if i == 0:
            master_feature_list = line[2:]
            continue
        feature_matrix += [map(float, line[3:])]
        survival_dct[line[0]] = (int(line[1]), float(line[2]))
    f.close()
    return np.array(feature_matrix), master_feature_list, survival_dct

# create_prosnet_input.py
def get_dictionary_symptom_herb_list():
    '''
    Returns a list of (symptom, herb) tuples.
    '''
    symptom_herb_list = []

    f = open('./data/herb_symptom_dictionary.txt', 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.strip().split('\t')

        line_length = len(line)
        # Some symptoms don't have good English translations.
        assert line_length == 2 or line_length == 5
        if line_length == 2:
            herb, symptom = line
        elif line_length == 5:
            herb, symptom, english_symptom, db_src, db_src_id = line
        symptom_herb_list += [(symptom, herb)]
    f.close()
    return list(set(symptom_herb_list))