### Author: Edward Huang

import os

### Generates the data and results directories.

def main():
    results_folder = './results'
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    data_folder = './data'
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)
    df_folder = './data/patient_dataframes'
    if not os.path.exists(df_folder):
        os.makedirs(df_folder)
    png_folder = './data/survival_plots'
    if not os.path.exists(png_folder):
        os.makedirs(png_folder)

if __name__ == '__main__':
    main()