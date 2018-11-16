#!python3
"""
Data download: TCGA_expression_pipeline.py
R preprocessing code: preprocess_classifier.R
Use multiple classification methods to classfy stage i-iv
breast cancer, using their fpkm expression values.
"""
import os
import random
import numpy as np
from libs.data_process import *
from libs.visualize import *
from libs.logistic import *
from libs.neural_network import *
from libs.svm import *
from libs.fm import *

os.chdir("C:/Users/jzhou/Desktop/")
random.seed(317)
np.random.seed(317)


def run():
    # Load data
    annot = load_annot(project="COAD")
    features, feature_names, labels, label_names = load_FPKM(
        filepath=f'./{annot[0]}_tumor_FPKM.csv', annot=annot, scaling=True)
    features, feature_names = feature_selection(
        features, labels, feature_names, anova=True, tree_based=False)
    train, test, train_labels, test_labels = split_dataset(
        features, labels, resample=False)
    # Logistic model
    logModel = logistic_regression(train, train_labels)
    logFig, logAx, logROC = visualize_model(
        logModel, test, test_labels, roc=True)
    # Neural network
    nnModel = neural_network(train, train_labels)
    nnLoss, nnAcc = evaluate_nn(nnModel, test, test_labels)
    # SVM model
    svcModel = svm(train, train_labels)
    svmFig, svmAx = visualize_model(svcModel, test, test_labels)
    # Factorization machine
    fmModel = fmc(train, train_labels)
    fmFig, fmAx = visualize_model(fmModel, test, test_labels)


if __name__ == '__main__':
    run()
