### Author: Edward Huang

import itertools
import subprocess

### This script reads the output of cluster_cancer_subtypes.py and outputs
### the feature types that produce significant results.

def script_call(metric):
    num_dim = '300'
    # TODO: add in num_dim after cluster_cancer_subtypes.py seq 50.
    # subprocess.call('python build_patient_feature_matrix.py %s %s' % (num_dim,
    #     s), shell=True)
    subprocess.call('python cluster_cancer_subtypes.py %s %s > '
        './results/%s_%s_out' % (metric, num_dim, metric, num_dim), shell=True)

    feature_list = ['tests', 'symptoms', 'herbs', 'drugs', 'history']
    x = itertools.chain.from_iterable(itertools.combinations(feature_list,
        r) for r in range(len(feature_list) + 1))
    x.next()
    with open('./results/%s_%s_out' % (metric, num_dim)) as f:
        while True:
            # Read 9 lines at a time.
            next_n_lines = list(itertools.islice(f, 6))
            if not next_n_lines:
                break
            feat_list = x.next()
            no_squam = next_n_lines[2].split()[1]
            squam = next_n_lines[5].split()[1]
            if float(no_squam) < 0.01 and float(squam) < 0.01:
                print metric, feat_list, no_squam, squam

def main():
    for metric in ('euclidean', 'minkowski', 'cityblock', 'sqeuclidean',
        'cosine', 'correlation', 'hamming', 'jaccard', 'chebyshev',
        'canberra', 'braycurtis'):
        script_call(metric)

if __name__ == '__main__':
    main()