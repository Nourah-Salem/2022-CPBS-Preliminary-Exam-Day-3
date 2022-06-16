# Using Principal Component Analysis (PCA) for Clustering of Metagenomic Sequences


> The objective of this project is to develop a general computational method for characterizing the viral communities in a metagenomic sample inhabiting a certain environment. Our aim is to computationally be able to classify metagenomic DNA sequences into their family groups accurately through finding the optimum size of K-mers to help us calculate the Principal Components Analysis (PCA) for the normalized frequencies of the K-mers. for mode details, please check the report project.

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

[Back To The Top](#read-me-template)

---

## How To Use

#### Installation
In order to run all the scripts (including the validation one), please make sure the following packages are installed:
1. Pandas
2. Seaborn
3. sklearn (optional, in case of reuning the validation method)


#### Steps to run the scripts:
First, Convert the FASTA and FASTQ data files for the viral and bacterial communities (repectively) to CSVs:   
```html
    python ./Preprocessing/FASTA_CSV.py
```
```html
    python ./Preprocessing/FASTQ_CSV.py
```
this will generate the CSV files in the Processed_data folder
Second, run the K-mer builder/ PCA model:
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
```html
    python ./Tests/Method_Uint_tests.py
```

[Back To The Top](#read-me-template)

---

## Author Info

- Twitter - [@jamesqquick](https://twitter.com/jamesqquick)
- Website - [James Q Quick](https://jamesqquick.com)

[Back To The Top](#read-me-template)
