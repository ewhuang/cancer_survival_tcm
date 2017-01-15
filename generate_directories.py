### Author: Edward Huang

import os

### Generates the data and results directories.

def main():
    results_folder = './results'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    png_folder = './results/survival_plots_single_features'
    if not os.path.exists(png_folder):
        os.makedirs(png_folder)
    png_folder = './results/survival_plots_multiple_features'
    if not os.path.exists(png_folder):
        os.makedirs(png_folder)
    png_folder = './results/survival_plots_synergy_features'
    if not os.path.exists(png_folder):
        os.makedirs(png_folder)
    features_folder = './results/feature_analyses'
    if not os.path.exists(features_folder):
        os.makedirs(features_folder)
    data_folder = './data'
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)
    df_folder = './data/patient_dataframes_single_features'
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)
    df_folder = './data/patient_dataframes_multiple_features'
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)
    df_folder = './data/patient_dataframes_synergy_features'
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)

if __name__ == '__main__':
    main()