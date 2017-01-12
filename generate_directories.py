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

if __name__ == '__main__':
    main()