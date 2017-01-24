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
    df_folder = './data/patient_dataframes_%s' % plot_type
    plot_folder = './results/survival_plots_%s' % plot_type
    if isProsnet:
        df_folder += '_%s' % num_dim
        plot_folder += '_%s' % num_dim
    else:
        df_folder += '_none'
        plot_folder += '_none'
    # Delete the folder if this is a fresh run.
    if os.path.exists(df_folder):
        shutil.rmtree(df_folder)
    if os.path.exists(plot_folder):
        shutil.rmtree(plot_folder)
    # Remake the folder.
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

def write_clusters():
    '''
    Perform clustering and write out the dataframe to file.
    '''
    # Get the feature matrix.
    f = './data/feature_matrices/feature_matrix'
    if isProsnet:
        f += '_%s' % num_dim
    f += '.txt'
    feature_matrix, master_feature_list, survival_mat = read_feature_matrix(f)

    # For treatment plot type, independency list is set of both symptoms and 
    # syndromes, since both clusters must have it. For synergy, the independency
    # list is the set of drugs.
    indep_list = []
    if plot_type == 'treatment':
        for fname in ('cancer_other_info_mr_symp', 'cancer_syndrome_syndromes'):
            indep_list += read_spreadsheet('./data/%s.txt' % fname)[1]
        indep_list = set(indep_list).intersection(master_feature_list)
    elif plot_type == 'synergy':
        indep_list = read_spreadsheet('./data/cancer_drug_2017_sheet2.txt')[1]
    # Treatment plot type, dependency list is the set of treatments. For
    # synergy plot type, dependency list is the set of herbs.
    dep_list = []
    if plot_type == 'treatment':
        for fname in ('cancer_other_info_herbmed', 'cancer_drug_2017_sheet2'):
            dep_list += read_spreadsheet('./data/%s.txt' % fname) [1]
        dep_list = set(dep_list).intersection(master_feature_list)
    elif plot_type == 'synergy':
        dep_list = read_spreadsheet('./data/cancer_other_info_herbmed.txt')[1]

    for indep_node in indep_list:
        indep_idx = master_feature_list.index(indep_node)
        for dep_node in dep_list:
            idx_list = [indep_idx, master_feature_list.index(dep_node)]
            labels = []
            new_feature_matrix = feature_matrix[:,idx_list]
            for row in new_feature_matrix:
                # Invalid patient if they don't have the drug/condition.
                if row[0] == 0:
                    labels += [-1]
                # Patient had condition/drug, and treated with the dep_node.
                # TODO: Vary the definition of treated/not treated.
                elif row[1] == 1:
                    if plot_type == 'treatment':
                        labels += ['treated']
                    elif plot_type == 'synergy':
                        labels += ['both']
                else:
                    if plot_type == 'treatment':
                        labels += ['not_treated']
                    elif plot_type == 'synergy':
                        labels += ['drug_only']

            # Skip lopsided clusterings.
            num_indep = labels.count('treated') + labels.count('both')
            num_dep = labels.count('not_treated') + labels.count('drug_only')
            # TODO: Change the minimum number of clusters.
            min_patients = max(0.2 * (num_indep + num_dep), 10)
            if num_indep < min_patients or num_dep < min_patients:
                continue

            out = open('%s/%s_%s.txt' % (df_folder, indep_node, dep_node), 'w')
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
    if len(sys.argv) not in [2, 3]:
        print 'Usage: python %s treatment/synergy num_dim<optional>' % (
            sys.argv[0])
        exit()
    global plot_type, isProsnet
    plot_type, isProsnet = sys.argv[1], False
    assert plot_type in ['treatment', 'synergy']
    if len(sys.argv) == 3:
        global num_dim
        num_dim, isProsnet = sys.argv[2], True

    generate_directories()

    # Cluster the patients.
    write_clusters()
    # Call the R script.
    command = 'Rscript survival_model.R %s ' % plot_type
    if isProsnet:
        command += num_dim
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))