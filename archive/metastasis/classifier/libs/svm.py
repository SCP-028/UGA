#!python3
import numpy as np
from math import log2
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC


def svc_param_selection(X, y, nfolds=3,
                        C_range=(-5, 15), gamma=(-15, 3)):
    """Perform grid search on SVM hyperparameters.

    Parameters
    ----------
        X: array
            Training set
        y: array
            Training labels
        nfolds: int, optional
            Cross validation folds
        C_range: tuple, optional
            A tuple of two integers indicating the range of C values to be
            searched.
        gamma: tuple, optional
            The range of gamma values to be searched for kernel rbf.

    Return
    ------
        A [GridSearchCV] object containing the best parameters.
    """
    tuned_parameters = [
        {
            'kernel': ['rbf'],
            'C': [2 ** i for i in range(C_range[0], C_range[1])],
            'gamma': [2 ** i for i in range(gamma[0], gamma[1])]
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
    grid_search = GridSearchCV(SVC(decision_function_shape='ovo'),
                               tuned_parameters, cv=nfolds,
                               scoring='accuracy', n_jobs=3, verbose=1)
    grid_search.fit(X, y)
    params = grid_search.best_params_
    if params['kernel'] == 'rbf':
        tuned_parameters = [
            {
                'kernel': ['rbf'],
                'C': [2 ** i for i in np.linspace(log2(params['C']) - 1,
                                                  log2(params['C']) + 1,
                                                  endpoint=True)],
                'gamma': [2 ** i for i in np.linspace(log2(params['gamma']) - 1,
                                                      log2(params['gamma']) + 1,
                                                      endpoint=True)]
            }
        ]
        grid_search = GridSearchCV(SVC(decision_function_shape='ovo'),
                                   tuned_parameters, cv=nfolds,
                                   scoring='accuracy', n_jobs=3, verbose=1)
        grid_search.fit(X, y)
        params = grid_search.best_params_
    print(f"{params}\nBest score: {grid_search.best_score_}")
    return grid_search


# SVM #
def svm(train, train_labels):
    """Tune the hyperparameters and classify using SVM.

    Parameters
    ----------
        train: array
            Input data from function `split_dataset`.
        train_labels: array
            Labels from function `split_dataset`.

    Return
    ------
        A tuned and trained SVM model.
    """
    svc_params = svc_param_selection(train, train_labels)
    svc = svc_params.best_estimator_
    svcModel = svc.fit(train, train_labels)
    return svcModel
