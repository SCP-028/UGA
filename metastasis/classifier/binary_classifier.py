#!python3
"""
Data download: TCGA_expression_pipeline.py
R preprocessing code: preprocess_classifier.R
Use multiple classification methods to classfy stage i-iii
and stage iv breast cancer, using their fpkm expression values.
"""

import os
import random
import re
import math
import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score, classification_report

random.seed(1005)
os.chdir("C:/Users/jzhou/Desktop/")

# annotation table #
annot = pd.read_csv("./expression_FPKM/annotation/annot.tsv")
annot = annot.loc[(annot.project == 'BRCA') &
                  (annot.sample_type == 'tumor')]
annot_iv = annot.loc[annot.tumor_stage.str.contains(
    r'iv', na=False), 'barcode']
annot_i_iii = annot.loc[annot.tumor_stage.str.contains(
    r'\si{1,3}[abc]?', na=False), 'barcode']

# expression data #
df = pd.read_csv("./BRCA_DEG.csv", index_col='hgnc_symbol')

# separate stage i-iii and stage iv groups
df_iv = df.loc[:, df.columns.isin(annot_iv)]
df_i_iii = df.loc[:, df.columns.isin(annot_i_iii)]
df = pd.concat([df_i_iii, df_iv], axis=1)

# change form to what sklearn needs and split out test set #
label_names = np.array(['i_iii', 'iv'])
labels = np.array([0 for x in df_i_iii.columns] +
                  [1 for x in df_iv.columns])
feature_names = np.array(df_iv.index)
features = df.transpose().as_matrix()
features = scale(features)
train, test, train_labels, test_labels = train_test_split(features,
                                                          labels,
                                                          test_size=0.20,
                                                          random_state=1005)

# over-sampling on training set #
sm = SMOTE(random_state=1005, ratio='minority', kind='svm')
train_res, train_labels_res = sm.fit_sample(train, train_labels)


def pretty_print_linear(coefs, names=None, sort=False, dataFrame=False):
    """ Output function for lineaer models.

    Args:
        coefs [array]: coefficients extracted from models
        names [array]: feature names
        sort [bool]: whether to remove zeros in the output
        dataFrame [bool]: determine output as dataFrame or formula

    Return:
        either a formula (default) or a dataFrame of features and cofficients.
    """
    coefs = coefs.flatten()
    if names is None:
        names = [f"X{x}" for x in range(len(coefs))]
    lst = dict(zip(names, coefs))
    if sort:
        lst = {name: coef for name, coef in lst.items() if coef != 0}
    if dataFrame:
        formula = pd.Series(lst, name='coef')
        formula.index.name = 'symbol'
    else:
        formula = " + ".join(f"{coef:{1}.{4}} * {name}" for name,
                             coef in lst.items())
        formula = re.sub(r'\+\s-', '- ', formula)
    return(formula)


# logistic regression #
logistic = LogisticRegressionCV(penalty='l1', solver='saga', scoring='f1_weighted',
                                class_weight='balanced', max_iter=300,
                                n_jobs=3, verbose=1, random_state=1005)
logModel = logistic.fit(train_res, train_labels_res)
logPreds = logModel.predict(test)
print(f"AUC: {roc_auc_score(test_labels, logPreds)}")
print(classification_report(test_labels, logPreds))

# Elastic Net #
glmnet = ElasticNetCV(l1_ratio=0.5, max_iter=10000, verbose=1,
                      n_jobs=3, random_state=1005, selection='cyclic')
netModel = glmnet.fit(train_res, train_labels_res)
netPreds = netModel.predict(test)
print(f"AUC: {roc_auc_score(test_labels, netPreds)}")


def svc_param_selection(X, y, nfolds=3):
    """Perform grid search on SVM hyperparameters.

    Args:
        X [array]: training set
        y [array]: training labels
        nfolds [int]: cross validation folds

    Return:
        a [GridSearchCV] object containing the best parameters.
    """
    tuned_parameters = [
        {
            'kernel': ['rbf'],
            # 'C': [2 ** i for i in range(-5, 15)],
            'C': [2 ** i for i in np.linspace(-5, -2, endpoint=True)],
            # 'gamma': [2 ** i for i in range(-15, 3)]
            'gamma': [2 ** i for i in np.linspace(-15, -11, endpoint=True)]
        },
        {
            'kernel': ['linear'],
            'C': [1, 10, 100, 1000]
        },
        {
            'kernel': ['poly'],
            'degree': [2, 3],
            'gamma': [0.001, 0.01, 0.1, 1]
        }
    ]
    grid_search = GridSearchCV(SVC(), tuned_parameters, cv=nfolds,
                               scoring='roc_auc', n_jobs=3, verbose=1)
    grid_search.fit(X, y)
    print(f"{grid_search.best_params_}\nBest score: {grid_search.best_score_}")
    return grid_search


# SVM #
svc_params = svc_param_selection(train_res, train_labels_res)
svc = svc_params.best_estimator_
svcModel = svc.fit(train_res, train_labels_res)
svcPreds = svcModel.predict(test)
print(f"AUC: {roc_auc_score(test_labels, svcPreds)}")
print(classification_report(test_labels, svcPreds))

# Neural Network #
