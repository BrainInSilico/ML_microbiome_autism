# SVM_microbiome_autism
predicting Autism disorder spectrum with SVM and Random Forest classifiers

This project was a quick start project for a customer who wanted to see the capabilities of "classical" machine learning predictors on microbiome data

## Data
The original dataset comes [from Kaggle](https://www.kaggle.com/datasets/antaresnyc/human-gut-microbiome-with-asd?select=GSE113690_Autism_16S_rRNA_OTU_assignment_and_abundance.csv) : 

## Scripts
dataCleaning.R contains a preprocessing of the the dataset
* extracting abundances
* transforming them to relative abundances
* removing OTU with less than 5% non-zero data over all dataset
* selecting the most informative OTU

evaluatePrediction.R build 10 subsets for a 10-fold cross-validation. 
For each fold, it create a SVM model and a random Forest model, make a prediction on the test set and evaluate it with AUC

## Results
The average expected AUC is 0.693 for SVM and 0.976 for Random Forest
