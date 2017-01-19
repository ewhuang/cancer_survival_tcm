#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

from file_operations import read_spreadsheet, get_dictionary_symptom_herb_list
import time

### This script prepares the input files for Prosnet. Nodes need to be in a
### separate file.

# Records the set of nodes already written to file.
written_nodes = set([])

def get_protein_herb_edge_list():
    '''
    Get a list of herb-protein edges.
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
    return protein_herb_edge_list

def get_coocc_edge_list(fname_a, fname_b):
    '''
    Get the symptom-herb relations from the co-occurrence counts of the patient
    records.
    '''
    coocc_edge_list = []
    # Read the respective patient dictionaries.
    patient_dct_a, feature_a_list = read_spreadsheet(fname_a)
    patient_dct_b, feature_b_list = read_spreadsheet(fname_b)
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
    amongst the different types of edges.
    '''
    global written_nodes
    for (node_a, node_b) in edge_list:
        # Edge weights are always 1.
        edge_out.write('%s\t%s\t1\t%s%s\n' % (node_a, node_b, node_type_a,
            node_type_b))
        if node_a not in written_nodes:
            written_nodes.add(node_a)
            node_out.write('%s\t%s\n' % (node_a, node_type_a))
        if node_b not in written_nodes:
            written_nodes.add(node_b)
            node_out.write('%s\t%s\n' % (node_b, node_type_b))

def main():
    protein_herb_list = get_protein_herb_edge_list()

    # Mapping letters to files.
    f_tuples = [('sym', './data/cancer_other_info_mr_symp.txt'), ('h', './data/'
        'cancer_other_info_herbmed.txt'), ('syn', './data/cancer_syndrome_'
        'syndromes.txt'), ('d', './data/cancer_drug_2017_sheet2.txt')]
    
    node_out = open('./data/dca_node_list.txt', 'w')
    edge_out = open('./data/dca_edge_list.txt', 'w')
    # Start off by writing out the protein-herb list.
    write_files(node_out, edge_out, protein_herb_list, 'p', 'h')
    # Loop through every pair of node types.
    for i in range(len(f_tuples)):
        node_type_a, fname_a = f_tuples[i]
        for j in range(i + 1, len(f_tuples)):
            node_type_b, fname_b = f_tuples[j]
            # Get the co-occurrence edge list.
            edge_list = get_coocc_edge_list(fname_a, fname_b)

            # symptom-herb edges should add the medical textbook's edges.
            if node_type_a == 'sym' and node_type_b == 'h':
                edge_list += get_dictionary_symptom_herb_list()

            # Write the edges out to file.
            write_files(node_out, edge_out, edge_list, node_type_a, node_type_b)
    edge_out.close()
    node_out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))