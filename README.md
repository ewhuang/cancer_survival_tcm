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

2.  Build the patient feature matrix. Currently removes features that appear in
    fewer than 10 patients (line 87/88).

    ```bash
    python build_patient_feature_matrix.py
    ```

3.  Run with 'single' argument first, then find the herbs that appear in
    ./results/survival_plots_single_features, then add them to line 56. Then,
    run this script with the 'multiple' argument. Results appear in ./results/
    survival_plots_multiple_features. Incorporates survival_model.R.

    ```bash
    python cluster_patient_survival.py kmeans num_clusters single/multiple
    ```

4.  Survival model on two clusters: the first contains patients that take a drug
    A, but not an herb B. The second cluster contains patients tha take both.
    Also incorporates survival_model.R.

    ```bash
    python drug_herb_synergistic_survival.py
    ```

## Prosnet Data Imputation.
Download PPI edges by going to 

Convert the HGNC IDs of hit_herb_target relations_zhou.xlsx by going to
http://www.genenames.org/cgi-bin/download, and then checking only Approved Symbol and Entrez Gene ID, then clicking submit. Save results into ./data/hgnc_to_entrez.txt.
To get the actual PPI network, download http://www.functionalnet.org/humannet/HumanNet.v1.benchmark.txt. Move to ./data/.

1.  Create the input for prosnet. From hit_herb_target relations_zhou.xlsx, 
    copy the filtered herb-target relations spreadsheet into a new file,
    ./data/herb_protein_relations.txt.

    ```bash
    python create_prosnet_input.py num_dim
    ```