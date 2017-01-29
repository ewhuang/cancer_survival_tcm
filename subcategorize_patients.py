#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from file_operations import read_feature_matrix
import os
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans, SpectralClustering, AgglomerativeClustering
import subprocess
import sys
import time

def generate_directories():
    global df_folder
    df_folder = './data/patient_dataframes_%s' % clus_meth
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)
    plot_folder = './results/survival_plots_%s' % clus_meth
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

def write_clusters(suffix):
    # Get the feature matrix.
    f = './data/feature_matrices/feature_matrix%s.txt' % suffix
    feature_matrix, master_feature_list, survival_mat = read_feature_matrix(f)

    # TODO: Currently clustering on the cosine similarity matrix.
    feature_matrix = squareform(pdist(feature_matrix, metric='cityblock'))

    # Perform the clustering.
    if clus_meth == 'kmeans':
        # TODO: Set clustering type and number of clusters.
        est = KMeans(n_clusters=3, n_init=1000, random_state=930519)
    elif clus_meth == 'spectral':
        est = SpectralClustering(n_clusters=3, affinity='cosine')
    est.fit(feature_matrix)
    labels = list(est.labels_)

    # Determine the output filename.
    if suffix == '':
        out_name = '%s/without_prosnet.txt' % df_folder
    else:
        out_name = '%s/prosnet_%s.txt' % (df_folder, num_dim)

    out = open(out_name, 'w')
    out.write('death\ttime\tcluster\n')
    for patient_idx, patient_tup in enumerate(survival_mat):
        label = labels[patient_idx]
        # Don't write clusters with fewer than 10 patients.
        if labels.count(label) < 10:
            continue
        inhospital_id, death, time = patient_tup
        out.write('%d\t%g\t%s\n' % (death, time, label))
    out.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage: python %s cluster_method num_dim' % sys.argv[0]
        exit()
    global clus_meth, num_dim
    clus_meth, num_dim = sys.argv[1:]
    assert clus_meth in ['kmeans', 'spectral'] and num_dim.isdigit()

    generate_directories()

    # Cluster the patients.
    write_clusters('_%s' % num_dim) # Prosnet feature matrix suffix.
    write_clusters('') # Prosnet-less suffix.
    # # Call the R script.
    command = 'Rscript survival_model.R %s' % df_folder
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))