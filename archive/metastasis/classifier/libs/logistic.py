#!python3
from sklearn.linear_model import LogisticRegressionCV


def logistic_regression(train, train_labels, n_jobs=3,
                        score_method='f1_weighted', max_iter=4000):
    """Train a logistic regression model for multi-class classification.

    Parameters
    ----------
        train: array
            Input data from function `split_dataset`.
        train_labels: array
            Labels from function `split_dataset`.
        n_jobs: int, optional
            Number of jobs when fitting model.
        score_method: str, optional
            Choose from "f1_weighted", "roc_auc_score", "log_loss".
            Default is 'f1_weighted'.
        max_iter: int, optional
            Maximum number of iterations of the optimization algorithm.
            Default is 4000 (slow! LogisticRegressionCV default is 100).

    Return
    ------
        A fitted logistic model.
    """
    logistic = LogisticRegressionCV(penalty='l1', multi_class='multinomial',
                                    solver='saga', scoring=score_method,
                                    class_weight='balanced', max_iter=4000,
                                    n_jobs=n_jobs, verbose=1)
    model = logistic.fit(train, train_labels)
    return model
