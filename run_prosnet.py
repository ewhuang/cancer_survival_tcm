#!/usr/bin/python
# -*- coding: utf-8 -*-

### Author: Edward Huang

import argparse
import file_operations
import os
import subprocess
import string

### This script prepares the input files for Prosnet. Nodes need to be in a
### separate file.

# Records the set of nodes already written to file. Also unique number of edges.
global_node_set, global_edge_set, num_edge_types = set([]), set([]), 0

def get_protein_herb_edge_set():
    '''
    Get the set of (protein, herb) edges.
    '''
    protein_herb_edge_set = set([])
    f = open('./data/herb_protein_relations.txt', 'r')
    for line in f:
        line = line.strip().split('\t')
        if len(line) < 3:
            continue
        # Format: Scientific name, Chinese name, list of proteins.
        herb, protein_list = line[1], line[2:]
        # Remove non-Chinese characters from the string.
        herb = [c for c in herb.decode('utf-8') if u'\u4e00' <= c <= u'\u9fff']        
        herb = ''.join(herb).encode('utf-8')
        # Check that the string is all in Chinese.
        assert all(u'\u4e00' <= c <= u'\u9fff' for c in herb.decode('utf-8'))
        # Add the protein-herb edge.
        for protein in protein_list:
            protein_herb_edge_set.add((protein, herb))
    f.close()
    return protein_herb_edge_set

def get_ppi_edge_set():
    '''
    Get the protein-protein edge set from a PPI network.
    '''
    entrez_to_hgnc_dct = file_operations.get_entrez_to_hgnc_dct()

    ppi_edge_set = set([])
    # Gene ID's in this PPI network are Entrez ID's.
    f = open('./data/HumanNet.v1.benchmark.txt', 'r')
    for line in f:
        line = line.split()
        assert len(line) == 2
        node_a, node_b = line
        # Skip if no HGNC analogues.
        if node_a not in entrez_to_hgnc_dct or node_b not in entrez_to_hgnc_dct:
            continue
        # Translate the Entrez ID to HGNC protein.
        ppi_edge_set.add((entrez_to_hgnc_dct[node_a], entrez_to_hgnc_dct[node_b]))
    f.close()
    return ppi_edge_set

def get_patient_dct(fname):
    if 'smoking_history' in fname:
        patient_dct = file_operations.read_smoking_history()[0]
        age_fname = './data/cancer_caseinfo.txt'
        age_dct = file_operations.read_spreadsheet(age_fname)[0]
        # Tack on age to the medical history.
        for inhospital_id in patient_dct:
            if inhospital_id in age_dct:
                patient_dct[inhospital_id] += [(age_dct[inhospital_id][0])]
    else:
        fname = './data/%s.txt' % fname
        patient_dct = file_operations.read_spreadsheet(fname)[0]
    return patient_dct

def get_coocc_edge_set(patient_dct_a, patient_dct_b):
    '''
    Get the symptom-herb relations from the co-occurrence counts of the patient
    records.
    '''
    coocc_edge_set = set([])
    # Get the intersecting set of patients in both dictionaries.
    patient_set = set(patient_dct_a.keys()).intersection(patient_dct_b.keys())
    for inhospital_id in patient_set:
        for (node_a, node_a_freq) in patient_dct_a[inhospital_id]:
            for (node_b, node_b_freq) in patient_dct_b[inhospital_id]:
                # We skip self edges, but keep other same-type edges.
                if node_a != node_b:
                    coocc_edge_set.add((node_a, node_b))
    return coocc_edge_set

# def write_edges(node_out, edge_out, edge_set, node_type_a, node_type_b):
def write_edges(edge_out, edge_set, node_type_tup):
    '''
    Write the edges out to file. edge_label is just a letter to differentiate
    amongst the different types of edges. We write nodes out as integers, so
    map_out contains a word on each line corresponding to these integers.
    '''
    global global_node_set, global_edge_set, num_edge_types

    for edge in edge_set:
        # Check if edge has already been written.
        if edge in global_edge_set or edge[::-1] in global_edge_set:
            continue
        global_edge_set.add(edge)

        # Write out the edge.
        edge_label = string.ascii_lowercase[num_edge_types]
        for i in range(2):
            global_node_set.add((edge[i], 'b')) # TODO
            # global_node_set.add((edge[i], node_type_tup[i]))
            edge_out.write('%s\t%s\t1\t%s\n' % (edge[i], edge[1 - i], edge_label))
    # num_edge_types += 1 # TODO: uncomment.

def write_nodes():
    node_out = open('./data/prosnet_data/prosnet_node_list.txt', 'w')
    for node_pair in global_node_set:
        node_out.write('%s\t%s\n' % node_pair)
    node_out.close()

def run_prosnet():
    os.chdir('../prosnet')
    network_folder = '../cancer_survival_tcm/data/prosnet_data'

    command = ('./embed -node %s/prosnet_node_list.txt -link '
        '%s/prosnet_edge_list.txt -meta_path %s/meta.txt -output %s/prosnet_vectors_%s -binary 0 -size %s -negative 5 -samples 1 '
        '-iters 501 -threads 12 -model 2 -depth 10 -restart 0.8 '
        '-edge_type_num %d -train_mode 2' % (
            network_folder, network_folder, network_folder, network_folder,
            args.num_dim, args.num_dim, num_edge_types + 1))
    print command
    subprocess.call(command, shell=True)

    # Rename the resulting file, depending on whether we exclude treatments.
    # if args.excl_treat == None:
    #     out_fname = '%s/embed_%s_450_treatments.txt' % (network_folder, args.num_dim)
    # else:
    #     out_fname = '%s/embed_%s_450_no_treatments.txt' % (network_folder, args.num_dim)
    # os.rename('%s/embed_%s_450.txt' % (network_folder, args.num_dim), out_fname)

def parse_args():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--num_dim', help='number of ProSNet dimensions',
        required=True, type=int)
    parser.add_argument('-e', '--excl_treat', help='whether to exclude treatments (drugs and herbs)')
    args = parser.parse_args()

def main():
    parse_args()

    # Symptom file must always come before herb file here.
    f_tuples = [('m', 'cancer_other_info_mr_symp'), ('h', 'cancer_other_info_herbmed'),
        ('n', 'cancer_syndrome_syndromes'), ('d', 'cancer_drug_2017_sheet2'),
        ('t', 'cancer_check_20170324'), ('v', 'smoking_history')]
    # Exclude treatments if excl_treat is not None.
    if args.excl_treat != None:
        f_tuples = [tup for tup in f_tuples if tup[0] not in ['h', 'd']]

    # Create the folder and files for the ProSNet input network.
    input_folder = './data/prosnet_data'
    if not os.path.exists(input_folder):
        os.makedirs(input_folder)
    # node_out = open('%s/prosnet_node_list.txt' % input_folder, 'w')
    edge_out = open('%s/prosnet_edge_list.txt' % input_folder, 'w')

    # Start off by writing out the protein-herb list and PPI list.
    protein_herb_edge_set = get_protein_herb_edge_set()
    write_edges(edge_out, protein_herb_edge_set, ('p', 'h'))
    ppi_edge_set = get_ppi_edge_set()
    write_edges(edge_out, ppi_edge_set, ('p', 'p'))

    # Loop through every pair of node types.
    for i in range(len(f_tuples)):
        node_type_a, fname_a = f_tuples[i]
        patient_dct_a = get_patient_dct(fname_a)
        for j in range(i, len(f_tuples)):
            node_type_b, fname_b = f_tuples[j]
            patient_dct_b = get_patient_dct(fname_b)
            # Get the co-occurrence edge set.
            edge_set = get_coocc_edge_set(patient_dct_a, patient_dct_b)

            if node_type_a == 'm' and node_type_b == 'h':
                edge_set = edge_set.union(file_operations.get_dictionary_symptom_herb_set())

            # Write the edges out to file.
            write_edges(edge_out, edge_set, (node_type_a, node_type_b))
            
    edge_out.close()

    write_nodes()

    # Run prosnet. Outputs the low-dimensional vectors into files.
    run_prosnet()

if __name__ == '__main__':
    main()