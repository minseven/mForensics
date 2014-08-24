This repository is built for publicly providing the code used in the paper "Microbial Forensics: Predicting phenotypic characteristics and environmental conditions from large-scale gene expression profiles" by Kim, Zorraquino-Salvo, and Tagkopoulos under revision. The code is written in MATLAB and under GNU GPL V2 License. 

1. initializeEcoGEC.m : it preprocesses the compendium to get normalized, bias-adjusted and transformed expression levels into categorical values. The code requires to have the following files in the current working directory.

1.a ecogec_v1.conditions_v3.txt : condition file where row represents sample and column represents condition
1.b ecogec_v1.partition.txt : platform annotation for each sample. The order must be equivalent to 1.a
1.c ecogec_v1.dataset.adj.txt : the compendium data where row represents sample and column represents gene expression level

2. categorizeDataset.m : it categorizes expression levels by being called from initializeEcoGEC.m

3. crossValidateDataset.m : it reports training and testing of consensus-based prediction from the given dataset

useage: crossValidateCondition($dataset,$pheno,$nfold,$nclasses,$algorithm) where

$dataset for the preprocssed compendium produced from 1.
$pheno is correponding phenotype annotation.
$nfold is number of folds for cross-validation.
$nclasses is number classes for building a classifier.
$algorithm indicates the algorithm to be used for building a classifier (type 'All' for a consensus-based prediction, 'NB' for Naive-Bayes, 'SVM' for Support-Vector Machine, 'DT' for Decision Tree, 'KNN' for K-nearest neighbors)   

The command produces the following outputs:

$train_accuracy for training accuracies of $nfold folds
$test_accuracy for testing accuracies for $nfold folds
$parameter_spaces for parameters for each of gene sets selected by MI. 
$post_vector for posterior probability for each prediction if applicable
$r_vector for confidence score of each prediction
$totalcm for confusion matrix
$mi for mutual information for all genes

3. trainCondition.m : it trains with given dataset to find the gene set that maximizes the training performance. 

4. mutualInformation.m : it measures mutual information for each gene by being called from trainCondition.m

5. iterativeLearning.m : it imputes missing attributes by iteratively learning the attribute. 

