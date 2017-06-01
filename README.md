# cancer_survival_tcm
Cancer survival models for TCM patients.


## Data preprocessing.

1.  Download cancer_liujie.zip and cancer_drug_2017.1.12.xlsx from e-mails.
    Get cancer_other_info_mr_symp.txt and cancer_other_info_herbmed from the
    corresponding shsets of cancer_other_info.xlsx.
    Get cancer_life_days.txt from cancer_life_days.xlsx.
    Get cancer_syndrome_syndromes.txt from the syndromes sheet of cancer_syndrome.xlsx.
    Get incase_check.txt from incase_check.xlsx. Merge lines 2228/2229 and 3184/3185, 3328/3329.
    Get cancer_drug_2017_sheet2.txt from the .xlsx of the same name.

2.  Generate directories.
    
    ```bash
    python generate_directories.py
    ```

3.  Compiling ProSNet.

    ```bash
    cd ./prosnet/model
    make
    ```

## HEMnet Data Imputation.
Convert the HGNC IDs of hit_herb_target relations_zhou.xlsx by going to
http://www.genenames.org/cgi-bin/download, and then checking only Approved Symbol and Entrez Gene ID, then clicking submit. Save results into ./data/hgnc_to_entrez.txt.
To get the actual PPI network, download http://www.functionalnet.org/humannet/HumanNet.v1.benchmark.txt. Move to ./data/.

1.  Create the input for prosnet. From hit_herb_target relations_zhou.xlsx, 
    copy the filtered herb-target relations spreadsheet into a new file,
    ./data/herb_protein_relations.txt.

    ```bash
    python run_prosnet.py num_dim
    ```

2.  Build the patient feature matrix. -d and -s arguments optional, only for
    ProSNet enrichment.

    ```bash
    python build_patient_feature_matrix.py [-h] [-d NUM_DIM] [-s SIM_THRESH]
    ```

    Best run:
    ```bash
    python build_patient_feature_matrix.py -d 500 -s 0.3
    ```
    No arguments creates ./data/feature_matrices/feature_matrix_raw.txt
    ```bash
    python build_patient_feature_matrix.py
    ```
    
3.  Cancer subtyping. First categorize into cancer subtypes, then cluster into
    two clusters. 'partial' argument means we only use tests and symptoms in the
    clustering phase. Runs plot_kaplan_meiers.R internally.

    ```bash
    python cluster_cancer_subtypes.py [-h] [-d NUM_DIM] [-s SIM_THRESH] [-o OTHER_FEAT] [-p PARTIAL]
    ```

    Best run:
    ```bash
    python cluster_cancer_subtypes.py -d 500 -s 0.3 -p 1
    ```
    In plot_kaplan_meiers.R, change lines 25/27 in order to change the legend labels.

## Paper scripts.

1.  Sub-typing pipeline. First classify patients into the cancer groups, and then
    perform the clustering. Iterates through different similarity score thresholds
    and distance metrics. Holds the dimension constant, must change in file.
    ```bash
    python tune_feature_parameters.py
    ```
    Runs cluster_cancer_subtypes.py internally.

2.  Plots the dimension tuning effects with pre-entered plotting points.

    ```bash
    Rscript plot_dimension_sensitivity.R
    ```

## Extra experimental scripts. Nothing to do with HEMnet.
1.  Results appear in ./results/survival_plots. Incorporates survival_model.R.
    Using the additional num_dim argument means we are using the prosnet
    feature matrices. 'treatment' argument means we're studying patients with
    and without a treatment given a condition. 'synergy' argument means we're
    studying drug effects in/outside of the presence of an herb.
    'translateCharUTF8' must be called on a CHARSXP error might happen for
    R error corruption issues. So far, only way is to increase the threshold
    for the variable min_patients (currently 20).

    ```bash
    python dependency_survival_analysis.py treatment/synergy num_dim<optional>
    ```

    Best results when using binary features (line 30) and mean threshold (line 70) of build_patient_feature_matrix.py, as well as all prior information in run_prosnet.py

2.  Cluster patients using the command line argument clustering method. Clusters
    for both types of feature matrices: with and without Prosnet. Also runs the
    R script for plotting survival curves.
    
    ```bash
    python subcategorize_patients.py seq/full num_dim<optional>
    ```

    In the resulting feature analysis files in ./results/feature_p_values,
    the first line in sequential clustering indicates the symptoms/tests that
    were significantly different from other sub-clusters. <> signs indicate
    whether they were less than or greater than in value, respectively. The
    list shows the features that are significantly different between the two
    clusters within the sub-cluster.

## Post-experiment analysis

1.  Tunes the similarity threshold from 0.0 to 0.5 (any similarity score < this
    threshold is set to 0 in the similarity matrix). Plots the p-values of the
    survival fit curves vs. this parameter.

    ```bash
    python plot_threshold_sensitivity.py
    ```

2.  Plots the the top 10 features of with/without prosnet as computed by the t-
    test, vs. the survival times of the patients.

    ```bash
    python feature_vs_survival_plot.py
    ```

## Script for BIBM special issue BMC topic modeling.

1.  Generates a file with only patient syndromes, symptoms, and herbs.
    ```bash
    python bmc_preprocess.py
    ```