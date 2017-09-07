### Author: Edward Huang

import argparse
from file_operations import read_spreadsheet, read_smoking_history
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# This script creates the feature matrix inputs for clustering.

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--num_dim', help='Optional. Number of ProSNet dimensions.')
    parser.add_argument('-s', '--sim_thresh', help='Optional. Threshold for cosine similarity between ProSNet vectors.')
    args = parser.parse_args()
    if args.num_dim == None:
        assert args.sim_thresh == None
    else:
        assert args.num_dim.isdigit() and args.sim_thresh != None
    return args

def create_dct_lst():
    '''
    Creates a list of feature dictionaries for each spreadsheet. Also returns
    a list maintaining the unique features across all spreadsheets.
    '''
    feature_dct_list, master_feature_lst = [], []
    for fname in ('cancer_other_info_herbmed', 'cancer_other_info_mr_symp',
        'cancer_syndrome_syndromes', 'cancer_check_20170324',
        'cancer_drug_2017_sheet2', 'smoking_history', 'cancer_caseinfo'):
        # Smoking history has a separate reading function.
        if fname == 'smoking_history':
            feature_dct, feature_list = read_smoking_history()
        else:
            feature_dct, feature_list = read_spreadsheet('./data/%s.txt' %
                fname)
        # Update the list of feature dictionaries.
        feature_dct_list += [feature_dct]
        # Update the master feature list. Some symptoms occur in the tests, so
        # we must take care of duplicates.
        for feature in feature_list:
            if feature not in master_feature_lst:
                master_feature_lst += [feature]
    return feature_dct_list, master_feature_lst

def build_feature_matrix(feature_dct_list, master_feature_lst, patient_list):
    '''
    Takes the feature list and a dictionary, and build a feature matrix.
    '''
    feature_matrix = []

    for inhospital_id in patient_list:
        # Initialize the row for the patient.
        row = [0 for i in master_feature_lst]
        # Get the values from each of the feature dictionaries.
        for feature_dct in feature_dct_list:
            # An inhospital ID might not be in feature_dct with prosnet.
            if inhospital_id not in feature_dct:
                continue
            for (feature, feature_freq) in feature_dct[inhospital_id]:
                row[master_feature_lst.index(feature)] += feature_freq
        feature_matrix += [row]
    feature_matrix = np.array(feature_matrix)

    # Remove the bad columns.
    zero_column_indices = np.where(~feature_matrix.any(axis=0))[0]
    good_indices = [i for i in range(len(master_feature_lst)
        ) if i not in zero_column_indices]
    master_feature_lst = [e for i, e in enumerate(master_feature_lst
        ) if i in good_indices]

    return feature_matrix[:,good_indices], master_feature_lst

def impute_missing_data(feature_matrix, master_feature_lst, num_dim, sim_thresh):
    '''
    Given the feature matrix and the column labels (master_feature_lst), impute
    the missing feature data by getting the Prosnet vectors.
    '''
    def read_prosnet_output(master_feature_lst, num_dim):
        '''
        Reads the output low-dimensional vectors created by prosnet.
        '''
        vector_dct = {}
        # TODO: Currently taking last iteration number, 500. Can do other files.
        f = open('./data/prosnet_data/prosnet_vectors_%s_500' % num_dim, 'r')
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.split()
            feature, vector = line[0], map(float, line[1:])
            assert len(vector) == int(num_dim) and feature not in vector_dct
            vector_dct[feature] = vector
        f.close()
        # Reorganize the matrix according to the order of master_feature_lst.
        vector_matrix = []
        for feature in master_feature_lst:
            vector_matrix += [vector_dct[feature]]
        return np.array(vector_matrix)

    vector_matrix = read_prosnet_output(master_feature_lst, num_dim)
    similarity_matrix = np.abs(cosine_similarity(vector_matrix))

    similarity_matrix[similarity_matrix < sim_thresh] = 0
    # Refill diagonals with 1s in cases of rounding errors.
    np.fill_diagonal(similarity_matrix, 1)
    print similarity_matrix # TODO
    print np.array_equal(similarity_matrix, np.identity(len(similarity_matrix)))
    # Multiply the feature matrix and the similarity matrix.
    enriched_feature_matrix = np.dot(feature_matrix, similarity_matrix)

    return enriched_feature_matrix

def write_feature_matrix(feature_matrix, master_feature_lst, patient_list,
    survival_dct, fname_suffix):
    '''
    Writes the feature matrix out to file, along with column labels. First
    column should be inhospital_id's, and the second/third column should be the
    death/time event.
    '''
    out_fname = './data/feature_matrices/feature_matrix%s.txt' % fname_suffix
    out = open(out_fname, 'w')
    out.write('patient_id\tdeath\ttime\t%s\n' % '\t'.join(master_feature_lst))
    for i, row in enumerate(feature_matrix):
        inhospital_id = patient_list[i]
        death, time = survival_dct[inhospital_id][0] # one event per patient.
        out.write('%s\t%d\t%f\t%s\n' % (inhospital_id, death, time, '\t'.join(
            map(str, row))))
    out.close()

def print_sparsity_stats(feature_matrix):
    '''
    Prints average number of zero values for each patient.
    '''
    num_zeros = np.count_nonzero(feature_matrix==0)
    print 'average number of zeros:', num_zeros / float(feature_matrix.shape[0])
    print 'feature matrix shape:', feature_matrix.shape

def main():
    survival_dct = read_spreadsheet('./data/cancer_life_days.txt')[0]
    # Set keys as a list so it's consistent across dictionary reads.
    patient_list = list(survival_dct.keys())

    feature_dct_list, master_feature_lst = create_dct_lst()

    # Create the raw matrix, and remove bad columns.
    feature_matrix, master_feature_lst = build_feature_matrix(feature_dct_list,
        master_feature_lst, patient_list)

    args = parse_args()
    if args.num_dim == None:
        fname_suffix = '_raw'
    else:
        # Get the number of ProSNet dimensions and the cosine similarity threshold.
        num_dim, sim_thresh = args.num_dim, float(args.sim_thresh)
        feature_matrix = impute_missing_data(feature_matrix, master_feature_lst,
            num_dim, sim_thresh)
        # Add num_dim and sim_thresh to the filename suffix.
        fname_suffix = '_%s_%g' % (num_dim, sim_thresh)

    # Write out matrix out to file.
    write_feature_matrix(feature_matrix, master_feature_lst, patient_list,
        survival_dct, fname_suffix)

    # print_sparsity_stats(feature_matrix)

if __name__ == '__main__':
    main()