#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

import file_operations
import os
import subprocess
import string
import sys
import time

### This script prepares the input files for Prosnet. Nodes need to be in a
### separate file.

# Records the set of nodes already written to file. Also unique number of edges.
global_node_list, global_edge_list, num_edge_types = set([]), set([]), 0

def get_protein_herb_edge_list():
    '''
    Get a list of (protein, herb) edges.
    '''
    protein_herb_edge_list = []
    f = open('./data/herb_protein_relations.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        if len(line) < 3:
            continue
        herb, protein_list = line[1], line[2:]
        for protein in protein_list:
            protein_herb_edge_list += [(protein, herb)]
    f.close()
    return list(set(protein_herb_edge_list))

def get_ppi_edge_list():
    '''
    Get the protein-protein edge list from a PPI network.
    '''
    entrez_to_hgnc_dct = file_operations.get_entrez_to_hgnc_dct()

    ppi_edge_list = []
    # Gene ID's in this PPI network are Entrez ID's.
    f = open('./data/HumanNet.v1.benchmark.txt', 'r')
    for line in f:
        line = line.split()
        assert len(line) == 2
        node_a, node_b = line
        # Skip if no HGNC analogues.
        if node_a not in entrez_to_hgnc_dct or node_b not in entrez_to_hgnc_dct:
            continue
        ppi_edge_list += [(entrez_to_hgnc_dct[node_a],
            entrez_to_hgnc_dct[node_b])]
    f.close()
    return list(set(ppi_edge_list))

def get_coocc_edge_list(fname_a, fname_b):
    '''
    Get the symptom-herb relations from the co-occurrence counts of the patient
    records.
    '''
    coocc_edge_list = []
    # Read the respective patient dictionaries.
    patient_dct_a = file_operations.read_spreadsheet(fname_a)[0]
    patient_dct_b = file_operations.read_spreadsheet(fname_b)[0]
    # Get the intersecting list of patients in both dictionaries.
    patient_list = set(patient_dct_a.keys()).intersection(patient_dct_b.keys())
    for inhospital_id in patient_list:
        for node_a, node_a_freq in patient_dct_a[inhospital_id]:
            for node_b, node_b_freq in patient_dct_b[inhospital_id]:
                coocc_edge_list += [(node_a, node_b)]
    return list(set(coocc_edge_list))

def write_files(node_out, edge_out, edge_list, node_type_a, node_type_b):
    '''
    Write the edges out to file. edge_label is just a letter to differentiate
    amongst the different types of edges. We write nodes out as integers, so
    map_out contains a word on each line corresponding to these integers.
    '''
    global global_node_list, global_edge_list, num_edge_types
    # The order of node type a, b must match the edge order in edge_list.
    node_type_tup = (node_type_a, node_type_b)
    for edge in edge_list:
        node_a, node_b = [node.replace(' ', '_').replace('/', '_'
            ) for node in edge]
        # Check if edge has already been written.
        if edge in global_edge_list or (node_b, node_a) in global_edge_list:
            continue
        global_edge_list.add(edge)

        # Write out the nodes if they haven't been written out already.
        for i, node in enumerate(edge):
            if node not in global_node_list:
                global_node_list.add(node)
                node_out.write('%s\t%s\n' % (node, node_type_tup[i]))
        # Edge weights are all = 1. Map the edge type to a letter.
        edge_out.write('%s\t%s\t1\t%s\n' % (node_a, node_b,
            string.ascii_lowercase[num_edge_types]))
    num_edge_types += 1

def run_prosnet():
    results_folder = './data/prosnet_vectors'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    command = ('../simons_mouse/Sheng/prosnet/model/embed -node "%s/prosnet_'
        'node_list.txt" -link "%s/prosnet_edge_list.txt" -output "%s/prosnet_'
        'node_vectors_%s_dims.vec" -binary 0 -size %s -negative 5 -samples 1 '
        '-iters 100 -threads 24 -model 2 -depth 10 -restart 0.8 -edge_type_num '
        '%d -rwr_ppi 1 -rwr_seq 1 -train_mode 1' % (input_folder, input_folder,
        results_folder, num_dim, num_dim, num_edge_types))
    subprocess.call(command, shell=True)

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s num_dim' % sys.argv[0]
        exit()
    global num_dim, input_folder
    num_dim = sys.argv[1]
    assert num_dim.isdigit()

    # m: symptoms, n: syndromes, h: herbs, d: drugs, t: tests.
    f_tuples = [('m', './data/cancer_other_info_mr_symp.txt'), ('h', './data/'
        'cancer_other_info_herbmed.txt'), ('n', './data/cancer_syndrome_'
        'syndromes.txt'), ('d', './data/cancer_drug_2017_sheet2.txt'), ('t',
        './data/incase_check.txt')]

    input_folder = './data/prosnet_input'
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)
    node_out = open('%s/prosnet_node_list.txt' % input_folder, 'w')
    edge_out = open('%s/prosnet_edge_list.txt' % input_folder, 'w')

    # Start off by writing out the protein-herb list and PPI list.
    protein_herb_list = get_protein_herb_edge_list()
    write_files(node_out, edge_out, protein_herb_list, 'p', 'h')
    ppi_edge_list = get_ppi_edge_list()
    write_files(node_out, edge_out, ppi_edge_list, 'p', 'p')

    # Loop through every pair of node types.
    for i in range(len(f_tuples)):
        node_type_a, fname_a = f_tuples[i]
        for j in range(i + 1, len(f_tuples)):
            node_type_b, fname_b = f_tuples[j]
            # Get the co-occurrence edge list.
            edge_list = get_coocc_edge_list(fname_a, fname_b)

            # # symptom-herb edges should add the medical textbook's edges.
            # TODO: Currently not adding in textbook edges.
            # if node_type_a == 'm' and node_type_b == 'h':
            #     edge_list += file_operations.get_dictionary_symptom_herb_list()

            # Write the edges out to file.
            write_files(node_out, edge_out, edge_list, node_type_a, node_type_b)
            
    edge_out.close()
    node_out.close()

    # Run prosnet. Outputs the low-dimensional vectors into files.
    run_prosnet()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))