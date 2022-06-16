# Using Principal Component Analysis (PCA) for Clustering of Metagenomic Sequences


> The objective of this project is to develop a simple general computational method for characterizing the viral communities in a metagenomic sample inhabiting a certain environment. Our aim is to computationally be able to cluster metagenomic DNA sequences into their family groups accurately through finding the optimum size of K-mers that are ferquent enough to represent a certain taxonomy to calculate the Principal Components Analysis (PCA) to be able to visualize the similarities and differences between these communities. for mode details, please check the report project.

---

### Table of Contents

- [Description](#description)
- [How To Use](#how-to-use)
- [Example](#example)
---

## Description

![A workflow diagram.](https://github.com/Nourah-Salem/2022-CPBS-Preliminary-Exam-Day-3/blob/main/Images/2022-CPBS-Preliminary-Exam-Day-3%20workflow.png)

A workflow diagram, representing the steps for solving our problem, starting with data collection and preparation, followed by generating  our clustering feature, which is k-mer frequency, followed by using our clustering method, which is the principal component analysis (PCA).

#### Main Concepts Implemented 

- Normalised freqeuncy of k-mers
- Principal Components Analysis (PCA)


#### Model Inputs 
1. One FASTA file for bacterial metagenomic 500 samples
2. One FASTQ file for viral (SARS-CoV2) metagenomic 500 samples
#### Model Onputs
A graph represnting the clusters of both communities on the top two principal components

[Back To The Top](#read-me-template)

---

## How To Use

#### Installation
In order to run all the scripts (including the validation one), please make sure the following packages are installed:
1. Pandas
```html
    conda install -c anaconda pandas
```
2. Seaborn

```html
    conda install -c anaconda seaborn
```
3. sklearn (optional, in case of runing the validation method)
```html
    conda install -c anaconda scikit-learn
```

#### Steps to run the scripts:
First, please clone the repository directly from GitHub and make sure that you're in the project directory
Second, Convert the FASTA and FASTQ data files for the viral and bacterial communities (repectively) to CSVs:   
```html
    python ./Preprocessing/FASTA_CSV.py
```
```html
    python ./Preprocessing/FASTQ_CSV.py
```
this will generate the CSV files in the Processed_data folder
third, run the K-mer builder/ PCA model (default Kmer size is 3, can be changed on the functional call):
```html
    python PCA_of_Kmer_Frequency.py
```
This will generate the PCs for the 2 communities and save them in the Output folder.
Next, if you would like to validate the PCA manual implementation, please run the Sklearn PCA:
```html
    python PCA_Sklearn.py
```

[Back To The Top](#read-me-template)

---

## Unit Test
All our Unit tests are implemented in the following script:
```html
    python ./Tests/Method_Uint_tests.py
```

[Back To The Top](#read-me-template)

---

## Example:

We used our model on 2 main datasets, one representing the bacterial community and another representing the Viral SARS-Cov2, the size of K-mers selected was 3 to measure their frequencies in the genomes and apply the PCA. the following graph represent the first 2 PCs, discriminating the 2 communities:  
![output](https://github.com/Nourah-Salem/2022-CPBS-Preliminary-Exam-Day-3/blob/main/Images/pca_sklearn.png)
(we have done the same experiemnt but with different sizes of the kmers and they showed different behavior, results on them and discussion are presented in the report)

[Back To The Top](#read-me-template)
