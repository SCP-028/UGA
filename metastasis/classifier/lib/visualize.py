#!python3
import re
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix
from mlxtend.plotting import plot_confusion_matrix
from yellowbrick.classifier import ROCAUC


def visualize_model(model, test, test_labels, roc=False):
    """Visualization of fitted models.
    Parameters
    ----------
        model: model
            Must have model.predict method builtin.
        test: array
            From function `split_dataset`.
        test_labels: array
            From function `split_dataset`.
        roc: bool, optional
            Plot ROC-AUC curve.

    Return
    ------
        The confusion matrix and the ROC curve figures.
    """
    prediction = model.predict(test)
    # precision, recall, and f1-score
    print(classification_report(test_labels, prediction))
    # confusion matrix
    conf_mat = confusion_matrix(test_labels, prediction)
    fig, ax = plot_confusion_matrix(
        conf_mat, figsize=(7.5, 7.5), alpha=0.4)
    plt.show()
    # ROC curve(s)
    if roc:
        visualizer = ROCAUC(model)
        visualizer.score(test, test_labels)
        visualizer.poof()
        return fig, ax, visualizer
    else:
        return fig, ax


def pretty_print_linear(coefs, names=None, sort=False, dataFrame=False):
    """ Output function for lineaer models.

    Parameters
    ----------
        coefs: array
            Coefficients extracted from models
        names: array
            Feature names
        sort: bool, optional
            Whether to remove zeros in the output
        dataFrame: bool, optional
            Determine output as dataFrame or formula

    Return
    ------
        either a formula (default) or a dataFrame of features & cofficients.
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
