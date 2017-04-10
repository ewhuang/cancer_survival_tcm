### Author: Edward Huang

import itertools
import subprocess

### This script reads the output of cluster_cancer_subtypes.py and outputs
### the feature types that produce significant results.

def script_call(metric, s):
    # TODO: add in 50 after cluster_cancer_subtypes.py seq 50.
    # subprocess.call('python build_patient_feature_matrix.py 100 %s' % s,
    #     shell=True)
    subprocess.call('python cluster_cancer_subtypes.py seq %s > ./results/%s_out'
        % (metric, metric), shell=True)

    feature_list = ['tests', 'symptoms', 'herbs', 'drugs', 'history']
    x = itertools.chain.from_iterable(itertools.combinations(feature_list,
        r) for r in range(len(feature_list) + 1))
    x.next()
    i = 0
    with open('./results/%s_out' % (metric)) as f:
        while True:
            # Read 6 lines at a time.
            next_n_lines = list(itertools.islice(f, 6))
            if not next_n_lines:
                break
            i += 1
            feat_list = x.next()
            no_squam = next_n_lines[3].split()[1]
            squam = next_n_lines[5].split()[1]
            if float(no_squam) < 0.01 and float(squam) < 0.01:
                print metric, s, feat_list, next_n_lines[1].split()[1], no_squam, squam

def main():
    # for s in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
    for s in ['nah']:
        s = str(s)
        for metric in ('euclidean', 'minkowski', 'cityblock', 'sqeuclidean',
            'cosine', 'correlation', 'hamming', 'jaccard', 'chebyshev',
            'canberra', 'braycurtis'):
            script_call(metric, s)

if __name__ == '__main__':
    main()