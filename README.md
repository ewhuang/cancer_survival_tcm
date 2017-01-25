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

## Prosnet Data Imputation.
Convert the HGNC IDs of hit_herb_target relations_zhou.xlsx by going to
http://www.genenames.org/cgi-bin/download, and then checking only Approved Symbol and Entrez Gene ID, then clicking submit. Save results into ./data/hgnc_to_entrez.txt.
To get the actual PPI network, download http://www.functionalnet.org/humannet/HumanNet.v1.benchmark.txt. Move to ./data/.

1.  Create the input for prosnet. From hit_herb_target relations_zhou.xlsx, 
    copy the filtered herb-target relations spreadsheet into a new file,
    ./data/herb_protein_relations.txt.

    ```bash
    python run_prosnet.py num_dim
    ```

## Clustering patients for survival model
1.  Build the patient feature matrix. Adding the optional num_dim argument, we
    build a feature matrix by using prosnet vector similarity scores.

    ```bash
    python build_patient_feature_matrix.py num_dim<optional>
    ```

2.  Results appear in ./results/survival_plots. Incorporates survival_model.R.
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