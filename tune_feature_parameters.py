### Author: Edward Huang

import numpy as np
import itertools
import subprocess

### This script reads the output of cluster_cancer_subtypes.py and outputs
### the feature types that produce significant results.

def script_call(sim_thresh):
    # for num_dim in [100, 200, 300, 400, 500]:
    # TODO: add in num_dim after cluster_cancer_subtypes.py seq 50.
    subprocess.call('python build_patient_feature_matrix.py -d 500 -s %s' %
        sim_thresh, shell=True)
    # subprocess.call('python cluster_cancer_subtypes.py %s %s partial > '
    #     './results/%s_%s_out' % (metric, num_dim, metric, num_dim), shell=True)
    subprocess.call('python cluster_cancer_subtypes.py -d 500 -s %s' % (sim_thresh),
        shell=True)


def main():

    # for sim_thresh in np.arange(0.1, 1.0, 0.1):
    #     script_call(sim_thresh)

    f = open('out', 'r')
    # line = f.readline()
    for sim_thresh in np.arange(0.1, 1.0, 0.1):
        feature_list = ['tests', 'symptoms', 'herbs', 'drugs', 'history']
        feat_comb_list = itertools.chain.from_iterable(itertools.combinations(
            feature_list, r) for r in range(len(feature_list) + 1))
        for feat_comb in feat_comb_list:
            if feat_comb == ():
                continue
            one = f.readline().split()[1]
            two = f.readline().split()[1]
            if float(one) < 0.01 and float(two) < 0.01:
                print sim_thresh, feat_comb, one, two
    # while line != '':
    #     feat_set = feat_comb_list.next()
    #     print line.split()
    #     one = line.split()[1]
    #     two = f.readline().split()[1]
    #     if float(one) < 0.2 and float(two) < 0.2:
    #         print feat_set, one, two
    f.close()

if __name__ == '__main__':
    main()