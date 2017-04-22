### Author: Edward Huang

import itertools
import subprocess

### This script reads the output of cluster_cancer_subtypes.py and outputs
### the feature types that produce significant results.

def script_call(metric):
    for num_dim in [100, 200, 300, 400, 500]:
    # TODO: add in num_dim after cluster_cancer_subtypes.py seq 50.
    # subprocess.call('python build_patient_feature_matrix.py %s %s' % (num_dim,
    #     s), shell=True)
        subprocess.call('python cluster_cancer_subtypes.py %s %s partial > '
            './results/%s_%s_out' % (metric, num_dim, metric, num_dim), shell=True)


def main():
    for metric in ('euclidean', 'minkowski', 'cityblock', 'sqeuclidean',
        'cosine', 'correlation', 'braycurtis'):
        print metric
        script_call(metric)

if __name__ == '__main__':
    main()