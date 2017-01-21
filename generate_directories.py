### Author: Edward Huang

import os

### Generates the data and results directories.

def main():
    for folder_name in ('./results/', './data/', './results/feature_analyses',
        './results/survival_plots_single_features', './data/feature_matrices',
        './results/survival_plots_multiple_features',
        './results/survival_plots_synergy_features',
        './data/patient_dataframes_single_features',
        './data/patient_dataframes_multiple_features',
        './data/patient_dataframes_synergy_features'):
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

if __name__ == '__main__':
    main()