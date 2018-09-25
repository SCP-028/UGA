#!python3

import math
import pandas as pd
import numpy as np

from imblearn.over_sampling import SMOTE
from collections import Counter
from sklearn.preprocessing import scale
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import (VarianceThreshold, SelectFromModel,
                                       SelectFwe, SelectKBest, f_classif)
from sklearn.ensemble import ExtraTreesClassifier


def load_annot(filepath="./expression_FPKM/annotation/annot.tsv",
               project="COAD"):
    """Load annotation file and separate different stages.

    Parameters
    ----------
        filepath: str, optional
            Filepath to annot.tsv from `TCGA_expression_pipeline.py`.
        project : str, optional
            The project name to extract.

    Return
    ------
        A [list] of length 5, containing the project extracted, and
        annotation of different cancer stages.
    """
    annot = pd.read_csv(filepath)
    annot = annot.loc[(annot.project == project) &
                      (annot.sample_type == 'tumor')]
    annot1 = annot.loc[annot.tumor_stage.str.contains(
        r'\si[abc]?$', na=False), 'barcode']
    annot2 = annot.loc[annot.tumor_stage.str.contains(
        r'\si{2}[abc]?$', na=False), 'barcode']
    annot3 = annot.loc[annot.tumor_stage.str.contains(
        r'\si{3}[abc]?$', na=False), 'barcode']
    annot4 = annot.loc[annot.tumor_stage.str.contains(
        r'iv', na=False), 'barcode']
    return [project, annot1, annot2, annot3, annot4]


def load_FPKM(filepath, annot, scaling=True):
    """Load RNA-Seq FPKM data and separate into stages.

    Parameters
    ----------
        filepath: str
            The filepath to the expression values.
        annot: list
            Returned list from function `load_annot`.
        scaling: bool, optional
            Should the feature values be scaled or not.

    Return
    ------
        Feature values, feature names, labels, and label names.
        All of the above are numpy arrays.
    """
    df = pd.read_csv(filepath, index_col='hgnc_symbol')
    df1 = df.loc[:, df.columns.isin(annot[1])]
    df2 = df.loc[:, df.columns.isin(annot[2])]
    df3 = df.loc[:, df.columns.isin(annot[3])]
    df4 = df.loc[:, df.columns.isin(annot[4])]
    df = pd.concat([df1, df2, df3, df4], axis=1)
    features = df.transpose().as_matrix()
    feature_names = np.array(df4.index)
    labels = np.array([0 for x in df1.columns] +
                      [1 for x in df2.columns] +
                      [2 for x in df3.columns] +
                      [3 for x in df4.columns])
    label_names = np.array(['i', 'ii', 'iii', 'iv'])
    if scaling:
        features = scale(features)
    print(f'Original features shape: {features.shape}')
    return features, feature_names, labels, label_names


def feature_selection(features, labels, feature_names,
                      anova=True, tree_based=False):
    """Tree-based / K-best feature selection.

    Parameters
    ----------
        features: array
            From function `load_FPKM`.
        labels: array
            From function `load_FPKM`.
        feature_names: array
            From function `load_FPKM`.
        anova: bool, optional
            Whether to apply f_classif method.
        tree_based:bool, optional
            Whether to apply tree-based method.

    Return
    ------
        Selected features and their corresponding names.
    """
    sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
    features = sel.fit_transform(features, labels)
    feature_names = feature_names[sel.get_support()]
    print(f'{features.shape[1]} features selected using VarianceThreshold.')
    if anova:
        if tree_based:
            sel = SelectKBest(score_func=f_classif,
                              k=5 * round(math.sqrt(features.shape[1])))
        else:
            sel = SelectFwe(score_func=f_classif,
                            alpha=0.05)
        features = sel.fit_transform(features, labels)
        feature_names = feature_names[sel.get_support()]
        print(f'{features.shape[1]} Features selected with ANOVA.')
    if tree_based:
        clf = ExtraTreesClassifier(max_features='sqrt',
                                   class_weight='balanced',
                                   verbose=1)
        clf = clf.fit(features, labels)
        sel = SelectFromModel(clf, prefit=True)
        features = sel.transform(features)
        feature_names = feature_names[sel.get_support()]
        print(f'{features.shape[1]} Features selected with Extra Tree.')
    return features, feature_names


def split_dataset(features, labels, resample=False):
    """Spliting training and test set, and resample training set if needed.

    Parameters
    ----------
        features: array
            From function `load_FPKM`.
        labels: array
            From function `load_FPKM`.
        resample: bool, optional
            Whether to perform over-sampling on training set.

    Return
    ------
        training set, test set, and their corresponding labels.
    """
    train, test, train_labels, test_labels = train_test_split(features,
                                                              labels,
                                                              test_size=0.20)
    if resample:
        sm = SMOTE(ratio='auto', kind='svm')
        train, train_labels = sm.fit_sample(train, train_labels)
    print(f"Training data shape:\n{Counter(sorted(train_labels)).items()}")
    return train, test, train_labels, test_labels
