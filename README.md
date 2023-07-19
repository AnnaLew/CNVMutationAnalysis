# Mutation and CNV co-occurrence in cancer analysis

## Purpose of this code

This code was developed as part of the Bioinformatics Research Network project entitled "Mutation / CNV co-occurrence in cancer". The project aimed to determine whether mutations and copy number variations(CNVs) co-occur in cancer cells and in which cancer-gene pairs is this phenomenon most prevalent. 

## Data used

The initial data can be found in the `data/original` subfolder. The data used for this analysis are:
- Mutation data from DNA-Seq: mc3.v0.2.8.PUBLIC.maf.gene_vclass_HGVSp_sample_likelyDriverLoose.tsv
- CNV data: all_thresholded.by_genes_whitelisted.tsv.gz
- Clinical data: clinical_PANCAN_patient_with_followup.tsv

## Pipeline

The analysis is divided into 6 main steps. They make computing and understanding the script easier.

### 1. Preprocess

This script performs the initial preprocessing of the data. It makes the column names uniform and merges all the mutations and CNV data with corresponding patient data. Output files can be found in the `data/preprocessed` subfolder.

### 2. Process

This script processes the CNV data further by extracting the Hugo symbol, barcode and CNV information and merging it with the corresponding barcode data. The output file can be found in the `data/processed` subfolder.

### 3. MergeMutationsCNV

This script merges the mutation and CNV data into one data table. The output file can be found in the `data/merged` subfolder.

### 4. CountGeneCancerPairs

For each gene-cancer pair, this script counts the number of:
- mutations occurring alone,
- mutations occurring alongside positive CNVs,
- mutations occurring alongside negative CNVs,
- positive CNVs occurring alone,
- negative CNVs occurring alone,
- neither mutations nor CNV present.
 
The output file can be found in the `data/countedPairs` subfolder.

### 5. MakeGraphs

This script creates graphs showing the share of each abovementioned group in the total number of samples per gene-cancer pair. The graphs  can be found in the `data/graphs` subfolder.

### 6. Fishers test

This script performs Fisher's test to assess the statistical significance of the findings. The output file can be found in the `data/fishers` subfolder.
