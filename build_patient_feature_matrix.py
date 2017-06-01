### Author: Edward Huang

import argparse
from file_operations import read_spreadsheet, read_smoking_history
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
# from sklearn.preprocessing import normalize
# import sys

# This script creates the feature matrix inputs for clustering.

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
        # f = open('./data/prosnet_data/prosnet_node_vectors_%s_dims.vec' %
        #     num_dim, 'r')
        # TODO: iteration number.
        # f = open('../simons_mouse/Sheng/prosnet/model/embed_%s_450.txt' %
        #     num_dim, 'r')
        f = open('./data/prosnet_data/embed_%s_450.txt' % num_dim, 'r')
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
    # TODO: absolute value.
    similarity_matrix = np.abs(cosine_similarity(vector_matrix))
    # similarity_matrix = cosine_similarity(vector_matrix)

    similarity_matrix[similarity_matrix < sim_thresh] = 0
    # Remove non-diagonal 1s. TODO.
    # similarity_matrix[similarity_matrix == 1] = 0
    np.fill_diagonal(similarity_matrix, 1)

    # Multiply the feature matrix and the similarity matrix.
    enriched_feature_matrix = np.dot(feature_matrix, similarity_matrix)

    return enriched_feature_matrix

# def write_feature_matrix(feature_matrix, master_feature_lst, patient_list,
#     survival_dct, out_fname=''):
def write_feature_matrix(feature_matrix, master_feature_lst, patient_list,
    survival_dct, fname_suffix):
    '''
    Writes the feature matrix out to file, along with column labels. First
    column should be inhospital_id's, and the second/third column should be the
    death/time event.
    '''
    # if out_fname == '':
    #     results_folder = './data/feature_matrices'
    #     if isImputation:
    #         out_fname = '%s/feature_matrix_%s.txt' % (results_folder, num_dim)
    #     else:
    #         out_fname = '%s/feature_matrix.txt' % results_folder
    out_fname = './data/feature_matrices/feature_matrix%s.txt' % fname_suffix
    out = open(out_fname, 'w')
    out.write('patient_id\tdeath\ttime\t%s\n' % '\t'.join(master_feature_lst))
    for i, row in enumerate(feature_matrix):
        inhospital_id = patient_list[i]
        death, time = survival_dct[inhospital_id][0] # one event per patient.
        out.write('%s\t%d\t%f\t%s\n' % (inhospital_id, death, time, '\t'.join(
            map(str, row))))
    out.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--num_dim', help='Number of ProSNet dimensions.')
    parser.add_argument('-s', '--sim_thresh', help='Threshold for cosine similarity between ProSNet vectors.')

    args = parser.parse_args()
    if args.num_dim == None:
        assert args.sim_thresh == None
    else:
        assert args.sim_thresh != None

    return args

def main():
    # if len(sys.argv) not in [1, 3]:
    #     print 'Usage:python %s num_dim<optional> sim_thresh<optional>' % (
    #         sys.argv[0])
    #     exit()
    # global isImputation
    # isImputation = False

    # args = parse_args()
    # # if len(sys.argv) == 3:
    # if args.num_dim != None:
    #     # isImputation = True
    #     global num_dim, sim_thresh
    #     num_dim, sim_thresh = args.num_dim, float(args.sim_thresh)
    #     # num_dim, sim_thresh = sys.argv[1], float(sys.argv[2])
    #     assert num_dim.isdigit()

    survival_dct = read_spreadsheet('./data/cancer_life_days.txt')[0]
    # Set keys as a list so it's consistent across dictionary reads.
    patient_list = list(survival_dct.keys())

    feature_dct_list, master_feature_lst = create_dct_lst()

    # Create the numpy array, and remove bad columns.
    feature_matrix, master_feature_lst = build_feature_matrix(feature_dct_list,
        master_feature_lst, patient_list)

    # # Print the number of 0 values.
    # num_zeros = 0.0
    # for row in feature_matrix:
    #     for ele in row:
    #         if ele == 0:
    #             num_zeros += 1
    # print 'average number of zeros', num_zeros / float(feature_matrix.shape[0])
    # print 'feature matrix shape', feature_matrix.shape

    # Write out to file a unnormalized file.
    args = parse_args()
    if args.num_dim == None:
    # write_feature_matrix(feature_matrix, master_feature_lst, patient_list,
    #     survival_dct, './data/feature_matrices/unnormalized_feature_matrix.txt')
        write_feature_matrix(feature_matrix, master_feature_lst, patient_list,
            survival_dct, '_raw')

    # Normalize feature matrix. Currently not normalizing. Save it for clustering.
    # feature_matrix = normalize(feature_matrix, norm='l1')

    # Decide whether to perform ProSNet imputation.
    # if len(sys.argv) == 3:
    else:
        # isImputation = True
        # global num_dim, sim_thresh
        num_dim, sim_thresh = args.num_dim, float(args.sim_thresh)
        # num_dim, sim_thresh = sys.argv[1], float(sys.argv[2])
        assert num_dim.isdigit()
        feature_matrix = impute_missing_data(feature_matrix, master_feature_lst,
            num_dim, sim_thresh)
        # Add num_dim and sim_thresh to the filename suffix.
        fname_suffix = '_%s_%g' % (num_dim, sim_thresh)

        # Write out matrix out to file.
        write_feature_matrix(feature_matrix, master_feature_lst, patient_list,
            survival_dct, fname_suffix)

    # print feature_matrix.shape

if __name__ == '__main__':
    main()