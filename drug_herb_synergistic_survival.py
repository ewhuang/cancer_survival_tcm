#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from file_operations import read_spreadsheet
import numpy as np
import subprocess
import time

# This script separates patients into two clusters: patients that took a drug
# A and not an herb B, and patients that took both. Performs survival modeling
# on these two clusters.

def remove_frequency(dct):
    for key in dct:
        dct[key] = [pair[0] for pair in dct[key]]

def write_clusters(drug_only_list, drug_herb_list, survival_dct, drug, herb):
    '''
    Writes out the clusters for the R script.
    '''
    out = open('./data/patient_dataframes_synergy_features/%s_%s.txt' % (
        drug, herb), 'w')
    out.write('death\ttime\tcluster\n')
    # Cluster 0 for patients that have taken only the drug.
    for patient_id in drug_only_list:
        death, time = survival_dct[patient_id][0]
        out.write('%d\t%g\tdrug_only\n' % (death, time))

    # Cluster 1 for patients that have taken both the herb and the drug.
    for patient_id in drug_herb_list:
        death, time = survival_dct[patient_id][0]
        out.write('%d\t%g\tboth\n' % (death, time))
    out.close()

def cluster_patients():
    '''
    Separate the two types of patients for every possible drug-herb pair, and
    write out to file.
    '''
    survival_dct, dummy_return = read_spreadsheet('./data/cancer_life_days.txt')
    # Get the list of herbs and list of western drugs.
    drug_patient_dct, drug_list = read_spreadsheet('./data/cancer_drug_2017_'
        'sheet2.txt')
    herb_patient_dct, herb_list = read_spreadsheet('./data/cancer_other_info_'
        'herbmed.txt')

    # Disregard frequencies of attributes.
    remove_frequency(drug_patient_dct)
    remove_frequency(herb_patient_dct)

    # Loop through every possible drug-herb pair.
    for drug in drug_list:
        for herb in herb_list:
            # List of patient ID's.
            drug_only_list, drug_herb_list = [], []
            for patient in survival_dct:
                # Skip patient if they haven't taken drugs.
                if patient not in drug_patient_dct:
                    continue
                # Determine if a patient has taken the drug.
                patient_drugs = drug_patient_dct[patient]
                # Skip patient if they haven't taken this drug.
                if drug not in patient_drugs:
                    continue
                # Add patient to drug_only_list if patient hasn't taken herbs.
                if patient not in herb_patient_dct:
                    drug_only_list += [patient]
                    continue
                patient_herbs = herb_patient_dct[patient]
                # Check if the patient has taken this herb.
                if herb in patient_herbs:
                    drug_herb_list += [patient]
                else:
                    drug_only_list += [patient]
            # Only keep if both clusters have at least 10 patients each.
            if len(drug_only_list) < 10 or len(drug_herb_list) < 10:
                continue
            write_clusters(drug_only_list, drug_herb_list, survival_dct, drug,
                herb)

def main():
    cluster_patients()

    subprocess.call('Rscript survival_model.R synergy', shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))