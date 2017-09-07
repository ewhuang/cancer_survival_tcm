#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

### This script preprocesses the lung cancer data for use with the topic model.
### Extracts tumor types, symptoms, drugs, and herbs.

import file_operations

def get_patient_dct(fname):
    if 'smoking_history' in fname:
        patient_dct = file_operations.read_smoking_history()[0]
    else:
        fname = './data/%s.txt' % fname
        patient_dct = file_operations.read_spreadsheet(fname)[0]
    return patient_dct

def main():
    symp_dct = get_patient_dct('cancer_other_info_mr_symp')
    herb_dct = get_patient_dct('cancer_other_info_herbmed')
    # history_dct = get_patient_dct('smoking_history')
    synd_dct = get_patient_dct('cancer_syndrome_syndromes')

    # disease_dct = {}
    # for inhospital_id in history_dct:
    #     history = history_dct[inhospital_id]
    #     for ele, freq in history:
    #         assert freq != 0
    #         if '癌' in ele:
    #             if inhospital_id in disease_dct:
    #                 disease_dct[inhospital_id] += [ele]
    #             else:
    #                 disease_dct[inhospital_id] = [ele]

    patient_set = set(symp_dct.keys()).intersection(herb_dct.keys())
    # patient_set = patient_set.intersection(disease_dct)
    patient_set = patient_set.intersection(synd_dct)
    
    out = open('./results/bmc_data.txt', 'w')
    for inhospital_id in patient_set:
        # disease_list = set(disease_dct[inhospital_id])
        # out.write(','.join(disease_list) + '\t')
        synd_list = set([synd for synd, freq in synd_dct[inhospital_id]])
        out.write(','.join(synd_list) + '\t')
        symp_list = set([symp for symp, freq in symp_dct[inhospital_id]])
        out.write(','.join(symp_list) + '\t')
        herb_list = set([herb for herb, freq in herb_dct[inhospital_id]])
        out.write(','.join(herb_list) + '\n')
    out.close()

if __name__ == '__main__':
    main()