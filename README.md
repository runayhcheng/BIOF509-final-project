# BIOF509-final-project
# IBD/Crohn's Disease and cancer

### Files Included
1. Data: 
  -combined (folder):with all patient micrarrays (see How to Run Script for link to download)
  -codeToGene: code of genes on microarray 
  -sampleInfo: patient to sample match
  -rankedDataNoIndexoTransp.csv: ranked the expression of each gene wth percentile ranking (see How to Run Script for link to download)
2. Results:
  -codeOOP.py.ipynb : OOP coding with final output of the table showing dysregulatd gene
  -final_withallDataDisplayed.py : code with exploratory data analysis and all inter-step graphs shown
  -genes and cancer association: dysregulated gene found and references to association with cancer, if applicable
  -IEEE format final project

### Table of Contents
- Introduction 
- Data Description
- How to run scripts

### Introduction
IBD/Crohn's Disease is associated with colorectal cancer. The goal of this project is to use DBSCAN density clustering after applying UMAP
dimensionality reduction to find dysregulated genes that could lead to cancer, particularly if they happen in a particular subgroup.

### Data Description
I have microarray data from blood samples of 338 individuals, 261  of whom has IBD/Crohn's Disease. The gene labels to actual genes can be found in codeToGene.txt.

### How to run scripts
Download 'combined' folder and unzip from this link:https://drive.google.com/file/d/1Cd4VzmU8Jwh0EK9h9Sg5WPEaB1qyMJ5a/view?usp=sharing
Download 'rankedataNoIndexoTransp.csv from this link:https://drive.google.com/file/d/1YhkL1phAcq1gDl0gI0FhJKxWojimN8KJ/view?usp=sharing
To run codeOOP on the microarray data, run: `python codeOOP.py`. The paths of the different data have already been set but the os path require reset depending on the location of the downloads (i.e. 'combined' folder).
To see indivudal steps and output, please look at final_withallDataDisplayed.py (HTML file).
