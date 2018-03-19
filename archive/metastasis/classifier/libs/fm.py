#!python3
from sklearn.multiclass import OneVsRestClassifier
from polylearn import FactorizationMachineClassifier
from fastFM import als
from scipy import sparse


def fmc(train, train_labels):
    """Use factorization machines for multi-classification.

    Parameters
    ----------
        train: array
            Input data from function `split_dataset`.
        train_labels: array
            Labels from function `split_dataset`.

    Return
    ------
        A fitted factorization machine model.
    """
    model = FactorizationMachineClassifier()
    fmModel = OneVsRestClassifier(model, n_jobs=3)
    fmModel.fit(train, train_labels)
    return fmModel


def fast_FM(train, train_labels,
            n_iter=1000, init_stdev=0.1, rank=2,
            l2_reg_w=0.1, l2_reg_V=0.5):
    """Binary classification with factorization machines.

    Parameters
    ----------
    train: array
        Input data from function `split_dataset`.
    train_labels: array
        Labels from function `split_dataset`.
    n_iter: int, optional
        The number of samples for the MCMC sampler, number or iterations over
        the training set for ALS and number of steps for SGD.
    init_stdev: float, optional
           Sets the stdev for the initialization of the parameter
   rank: int, optional
       The rank of the factorization used for the second order interactions.
   l2_reg_w : float, optional
       L2 penalty weight for linear coefficients.
   l2_reg_V : float, optional
       L2 penalty weight for pairwise coefficients.
    """
    X = train
    X = sparse.csc_matrix(X)
    y = train_labels
    y[y != 3] = -1
    y[y == 3] = 1
    # TODO: Grid search to tune hyperparameters
    model = als.FMClassification(n_iter=n_iter, init_stdev=init_stdev,
                                 rank=rank, l2_reg_w=l2_reg_w, l2_reg_V=l2_reg_V)
    model.fit(X, y)
    # TODO: See if can deal with multiclass problems
    return model
