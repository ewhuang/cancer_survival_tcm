#!/usr/bin/python
# -*- coding: utf-8 -*-

### These functions are shared in a variety of scripts. Most are typically
### functions that read files.

from collections import OrderedDict
import numpy as np

def read_case_info():
    '''
    Get only the inhospital_id values that correspond to a first-time visit.
    '''
    first_time_id_list = []
    f = open('./data/cancer_caseinfo.txt', 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        if line[1] != '1':
            continue
        first_time_id_list += [line[0]]
    f.close()
    return first_time_id_list

# build_patient_feature_matrix.py
# drug_herb_synergistic_survival.py
# run_prosnet.py
def read_spreadsheet(fname):
    '''
    Depending on the spreadsheet, return a dictionary mapping the inhospital_id
    to the feature list. All feature lists should be of the same length; some
    are binary (like the survival labels) and some are frequency vectors.
    Key: inhospital_id -> str
    Value: feature, feature frequency tuple -> (str, float)
    '''
    first_time_id_list = read_case_info()
    feature_dct, unique_feature_list = OrderedDict({}), []
    f = open(fname, 'r')
    for i, line in enumerate(f):
        # Skip all header lines.
        if i == 0:
            continue
        line = line.strip().split('\t')
        if len(line) == 1:
            continue
        # Survival events. Output binary labels, along with the event length.
        if 'life_days' in fname:
            assert len(line) == 7
            inhospital_id, feature, feature_freq  = line[0], line[3], line[6]
            if '死亡' in feature: # 1 is death. 0 is unknown.
                feature = 1
            else:
                feature = 0
            feature_freq = float(feature_freq) / 30.0 # Convert to months.
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
        # elif 'incase_check' in fname:
        #     assert len(line) == 9
        #     inhospital_id, feature, feature_freq = line[0], line[3], line[4]
        #     # Skip negative tests.
        #     # TODO: good results with without cost condition.
        #     if feature_freq in ['无', '否', '没有'] or '费' in feature or (feature
        #         in ['住院天数', '阴影部位', '影像学诊断类型']):
        #         inhospital_id = 'bad_id'
        #         continue
        #     elif '级' in feature_freq:
        #         feature_freq = feature_freq[:feature_freq.index('级')]
        #     try:
        #         feature_freq = float(feature_freq)
        #         # Skip 0 features.
        #         if feature_freq == 0.0:
        #             inhospital_id = 'bad_id'
        #             continue
        #     except:
        #         feature_freq = 1
        elif 'cancer_check_20170324' in fname:
            if len(line) != 9:
                continue
            inhospital_id, feature, feature_freq = line[0], line[6], line[8]
            if feature_freq == '有':
                feature_freq = 1
            elif feature_freq == '无':
                continue
            elif '级' in feature_freq:
                feature_freq = feature_freq[:feature_freq.index('级')]
            try:
                feature_freq = float(feature_freq)
            except:
                # TODO: What to do with unfloatable things.
                # feature_freq = 1
                continue
        elif 'drug_2017' in fname:
            assert len(line) == 10
            inhospital_id, feature, feature_freq = line[0], line[1], line[4]
        # TODO: Currently adding in senior citizen status.
        elif 'cancer_caseinfo' in fname:
            assert len(line) == 7
            inhospital_id, feature, feature_freq = line[0], 'senior', line[6]
            if int(feature_freq) >= 65:
                feature_freq = 1
            else:
                continue
        else:
            print 'file_operations.py line 97: No such file!'
            exit()
        # Don't use second or later visits. TODO.
        if inhospital_id not in first_time_id_list:
            continue
        if inhospital_id not in feature_dct:
            feature_dct[inhospital_id] = []
        # Deal with folder path / issues.
        if type(feature) == str:
            feature = feature.replace('/', '_')
            feature = feature.replace(' ', '_')
        feature_dct[inhospital_id] += [(feature, float(feature_freq))]
        if feature not in unique_feature_list:
            unique_feature_list += [feature]
    f.close()
    return feature_dct, unique_feature_list

# subcategorize_patients.py
# def read_feature_matrix(fname):
def read_feature_matrix(suffix=''):
    '''
    Reads the feature matrix of the patient data. Takes an optional argument in
    the form of '_50'.
    '''
    feature_matrix, master_feature_list, survival_matrix = [], [], []
    # Process the filename.
    fname = './data/feature_matrices/feature_matrix%s.txt' % suffix

    f = open(fname, 'r')
    for i, line in enumerate(f):
        line = line.strip().split('\t')
        feature_list = line[3:]
        if i == 0:
            master_feature_list = feature_list[:]
            continue
        # Survival matrix consists of (patient, death, time) tuples.
        survival_matrix += [(line[0], int(line[1]), float(line[2]))]
        assert len(feature_list) == len(master_feature_list)
        feature_matrix += [map(float, feature_list)]
    f.close()
    return np.array(feature_matrix), master_feature_list, survival_matrix

# run_prosnet.py
def get_dictionary_symptom_herb_set():
    '''
    Returns a list of (symptom, herb) tuples.
    '''
    symptom_herb_set = set([])

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
        # Special herb cases.
        if '(' in herb:
            herb = herb[:herb.index('(')]
        elif '银翘片' in herb:
            herb = '银翘片'
        # Reformatting potential bad strings.
        symptom = symptom.replace('/', '_').replace(' ', '_')
        herb = herb.replace('/', '_').replace(' ', '_')
        symptom_herb_set.add((symptom, herb))
    f.close()
    return symptom_herb_set

# run_prosnet.py
def get_entrez_to_hgnc_dct():
    '''
    Gets mappings from HGNC ID's to Entrez ID's.
    '''
    entrez_to_hgnc_dct = {}
    f = open('./data/hgnc_to_entrez.txt', 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.strip().split('\t')
        if len(line) != 2:
            continue
        hgnc_id, entrez_id = line
        assert entrez_id not in entrez_to_hgnc_dct
        entrez_to_hgnc_dct[entrez_id] = hgnc_id
    f.close()    
    return entrez_to_hgnc_dct

# run_prosnet.py
def read_smoking_history():
    '''
    Read the patient history, including the smoking history.
    '''
    first_time_id_list = read_case_info()
    feature_dct = OrderedDict({})
    binary_feature_list = (['V肝炎病史', 'V高血压病史', 'V冠心病史', 'V家族遗传性疾病',
                # 'V输血史', 'V吸烟史', 'V药物过敏史', 'V有无慢性肺部疾病史', 'V中毒史',
                # TODO: eliminating some features.
                'V输血史', 'V吸烟史', 'V有无慢性肺部疾病史',
                'V肿瘤家族史'])
    cancer_list = []
    f = open('./data/smoking_history.txt', 'r')
    for i, line in enumerate(f):
        line = line.rstrip('\n\r').split('\t')
        assert len(line) == 47
        if i == 0:
            # Map the important features to their indices in the dataframe.
            header_line, feature_index_dct = line[:], {}
            for feature in binary_feature_list:
                feature_index_dct[feature] = header_line.index(feature)
            continue
        inhospital_id = line[0]
        if inhospital_id not in first_time_id_list:
            continue

        # Don't use second or later visits.
        assert inhospital_id not in feature_dct
        feature_dct[inhospital_id] = []

        # Go through each of the binary history features.
        for feature in binary_feature_list:
            feature_freq = line[feature_index_dct[feature]]
            if feature_freq == '有':
                feature_dct[inhospital_id] += [(feature, 1.0)]

        # Add the cancer type. Binarize the V病理类型 column.
        cancer_type = line[header_line.index('V病理类型')]
        if '癌' in cancer_type:
            feature_dct[inhospital_id] += [(cancer_type, 1.0)]
            if cancer_type not in cancer_list:
                cancer_list += [cancer_type]

        # Add metastasis site. Binarize these columns.
        for metastasis_col in ['V侵及或转移部位', 'V肿瘤并发症']:
            metastasis_list = line[header_line.index(metastasis_col)].split('、')
            for metastasis in metastasis_list:
                # The last list element is something like a Chinese space.
                if metastasis in ['董海涛', '卢雯平', '不详', '', '　']:
                    continue
                feature_dct[inhospital_id] += [(metastasis, 1.0)]
                if metastasis not in cancer_list:
                    cancer_list += [metastasis]
        
        # Tumor stage column.
        stage = line[header_line.index('VTNM分期')]
        # Skip unknown or bogus values.
        if '期' in stage and '不详' not in stage:
            stage_dct = {'IA':1, 'IB':2, 'II':3, 'IIA':4, 'IIB':5, 'III':6, 'IIIA':7, 'IIIB':8, 'IV':9}
            stage = stage[:stage.index('期')]
            stage = stage_dct[stage]
            feature_dct[inhospital_id] += [('VTNM分期', stage)]
            if 'VTNM分期' not in cancer_list:
                cancer_list += ['VTNM分期']
        # Add numerical features.
        # numerical_features = ['VKPS', 'VM', 'VN', 'VT', 'VW']
        numerical_features = ['VKPS', 'VM', 'VN', 'VT']
        for numerical_col in numerical_features:
            numerical_val = line[header_line.index(numerical_col)]
            if numerical_val == '0':
                continue
            try:
                feature_dct[inhospital_id] += [(numerical_col,
                    float(numerical_val))]
            except Exception:
                continue
    f.close()
    return feature_dct, binary_feature_list + cancer_list + numerical_features