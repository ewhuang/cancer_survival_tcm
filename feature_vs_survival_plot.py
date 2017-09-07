### Author: Edward Huang

from file_operations import read_feature_matrix
import matplotlib
import os
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

### This script plots the top 10 features of with/without prosnet as computed
### by the t-test, vs. the survival times of the patients.

def generate_folder():
    global out_folder
    out_folder = './results/feature_plots'
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

def read_feature_file(fname):
    '''
    Returns the top 10 features for the given filename's method.
    '''
    feature_list = []
    f = open(fname, 'r')
    for i, line in enumerate(f):
        if i == 10:
            break
        line = line.strip().split('\t')
        assert len(line) == 3
        feature_list += [line[0]]
    f.close()
    return feature_list

def plot_features_vs_survival(combined_feature_list):
    '''
    Given the list of features, for each feature, we plot their values vs.
    survival rate for each patient.
    '''
    # Only look at base feature matrix when plotting.
    f = './data/feature_matrices/feature_matrix.txt'
    feature_matrix, master_feature_list, survival_mat = read_feature_matrix(f)

    # Iterate through each feature.
    for feature in combined_feature_list:
        feat_idx = master_feature_list.index(feature)
        alive_list, dead_list = [], []
        # Iterate through each patient.
        for pat_idx, row in enumerate(feature_matrix):
            feature_value = row[feat_idx]
            if feature_value == 0: # Skip patients that don't have the feature
                continue
            inhospital_id, death, time = survival_mat[pat_idx]
            if death == 0: # TODO: Skipping survivors
                alive_list += [(time, feature_value)]
            else:
                dead_list += [(time, feature_value)]
        alive_list = sorted(alive_list, key=lambda x: x[0])
        dead_list = sorted(dead_list, key=lambda x: x[0])
        plt.title('Feature value vs. survival time')
        plt.xlabel('Survival time (months)')
        plt.ylabel('Feature value')
        plt.plot(*zip(*alive_list), color='r', label='alive')
        plt.plot(*zip(*dead_list), color='b', label='dead')
        plt.ylim(0, 1.2 * max([pt[1] for pt in alive_list + dead_list]))
        plt.xlim(0, 1.2 * max([pt[0] for pt in alive_list + dead_list]))
        plt.legend(loc='upper right')
        plt.show()
        pylab.savefig('%s/%s.png' % (out_folder, feature))
        plt.close()

def main():
    generate_folder()

    feat_folder = './results/feature_p_values_full'
    prosnet_features = read_feature_file('%s/prosnet_50.txt' % feat_folder)
    base_features = read_feature_file('%s/without_prosnet.txt' % feat_folder)

    combined_feature_list = list(set(prosnet_features).union(base_features))

    plot_features_vs_survival(combined_feature_list)

if __name__ == '__main__':
    main()