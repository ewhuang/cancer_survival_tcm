#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from file_operations import read_feature_matrix, read_spreadsheet
import numpy as np
import os
import shutil
import subprocess
import sys
import time

# This script takes the feature matrix, and, depending on the argument, uses
# a set of features to cluster the patients, then calls the R script to plot
# the survival models.

def generate_directories():
    global df_folder
    df_folder = './data/patient_dataframes'
    plot_folder = './results/survival_plots'
    if isProsnet:
        df_folder += '_%s' % num_dim
        plot_folder += '_%s' % num_dim
    else:
        df_folder += '_none'
        plot_folder += '_none'
    if os.path.exists(df_folder):
        shutil.rmtree(df_folder)
    if os.path.exists(plot_folder):
        shutil.rmtree(plot_folder)
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

def cluster_and_write():
    '''
    Perform clustering and write out the dataframe to file. This clustering
    only creates two clusters: one in which 
    '''
    f = './data/feature_matrices/feature_matrix'
    if isProsnet:
        f += '_%s' % num_dim
    f += '.txt'

    feature_matrix, master_feature_list, survival_mat = read_feature_matrix(f)
    min_patients = 0.2 * len(feature_matrix)

    # Get symptom/syndrome list.
    condition_list = []
    for fname in ('cancer_other_info_mr_symp', 'cancer_syndrome_syndromes'):
        condition_list += read_spreadsheet('./data/%s.txt' % fname)[1]
    condition_list = set(condition_list).intersection(master_feature_list)
    # Get the herb/drug list.
    treatment_list = []
    for fname in ('cancer_other_info_herbmed', 'cancer_drug_2017_sheet2'):
        treatment_list += read_spreadsheet('./data/%s.txt' % fname) [1]
    treatment_list = set(treatment_list).intersection(master_feature_list)

    for condition in condition_list:
        condition_idx = master_feature_list.index(condition)
        for treatment in treatment_list:
            idx_list = [condition_idx, master_feature_list.index(treatment)]
            labels = []
            new_feature_matrix = feature_matrix[:,idx_list]
            for row in new_feature_matrix:
                # Invalid patient if they don't have the condition in question.
                if row[0] != 0:
                    labels += [-1]
                # Patient had condition, but not treated with the treatment.
                # TODO: Vary the definition of treated/not treated.
                elif row[1] == 0:
                    labels += ['not_treated']
                else:
                    labels += ['treated']
            # Skip lopsided clusterings.
            if labels.count('not_treated') < min_patients or labels.count(
                'treated') < min_patients:
                continue

            out = open('%s/%s_%s.txt' % (df_folder, condition, treatment), 'w')
            out.write('death\ttime\tcluster\n')
            for patient_idx, patient_tup in enumerate(survival_mat):
                label = labels[patient_idx]
                # Skip patients that don't have the condition.
                if label == -1:
                    continue
                inhospital_id, death, time = patient_tup
                out.write('%d\t%g\t%s\n' % (death, time, label))
            out.close()

def main():
    if len(sys.argv) not in [1, 2]:
        print 'Usage: python %s num_dim<optional>' % sys.argv[0]
        exit()
    global isProsnet
    isProsnet = False
    if len(sys.argv) == 2:
        global num_dim
        num_dim, isProsnet = sys.argv[1], True
    generate_directories()

    # Cluster the patients.
    cluster_and_write()
    
    # Call the R script.
    command = 'Rscript survival_model.R '
    if isProsnet:
        command += num_dim
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))