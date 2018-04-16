import os

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as offline
import seaborn as sns
import statsmodels.api as sm
from sklearn.decomposition import PCA

FEVER_CUTOFF = 39
TIME_WINDOW = 2
FOURIER_COEF_NUM = 10
HODRICK_PRESCOTT_LAMBDA = 3000


def hodrick_prescott(df, lamb=1600):
    """Use the Hodrick-Prescott Filter to estimate the trend.

    Parameters
    ----------
        df: pandas DataFrame
            The temperature data at different timepoints
        lamb: float, optional
            The Hodrick-Prescott smoothing parameter. The larger it is the
            smoother the outcome gets.

    Returns
    -------
        df with a new column containing the trend, and the de-trended data.
    """
    df = df.sort_values("timepoint")
    _, df["trend"] = sm.tsa.filters.hpfilter(df.temperature, lamb=lamb)
    df["detrended"] = df["temperature"] - df["trend"]
    return df


def tsplot(df, uuid, save_image=False):
    """
    Parameters
    ----------
        df: pandas DataFrame
            Output of function `hodrick_prescott`.
        uuid: str
            Experiment and case identifier.
        save_image: bool, optional
            Whether to download the figure as a static file.

    Returns
    -------
        A long DataFrame containing the temperature, trend,
        and variance and mean of the deviation.
    """
    dfm = pd.melt(df, id_vars="timepoint")
    fig = {
        'data': [
            {
                'x': dfm[dfm['variable'] == "temperature"]['timepoint'],
                'y': dfm[dfm['variable'] == "temperature"]['value'],
                'name': "Temperature",
                'opacity': 0.7
            },
            {
                'x': dfm[dfm['variable'] == "trend"]['timepoint'],
                'y': dfm[dfm['variable'] == "trend"]['value'],
                'name': "Hodrick-Prescott filter trend",
                'opacity': 0.7
            },
            {
                'x': dfm[dfm['variable'] == "detrended"]['timepoint'],
                'y': dfm[dfm['variable'] == "detrended"]['value'],
                'name': "De-trended",
                "yaxis": "y2",
                'opacity': 0.7
            }
        ],
        'layout': {
            "title": uuid,
            'xaxis': {'title': 'Date'},
            'yaxis': {'title': "Temperature"},
            "yaxis2": {
                "title": "Temperature",
                "overlaying": "y",
                "side": "right"
            }
        }
    }
    if save_image:
        offline.iplot(fig, filename=f"{uuid}", show_link=False, image="png")
    else:
        offline.iplot(fig, filename=f"{uuid}", show_link=False)
    return dfm


def fourier_transform(dfs, t=2, n=10, detrend=True):
    """Perform a Fourier transform on the de-trended data.

    Parameters
    ----------
        df: pandas DataFrame
            trend_list: a series of output from the function `hodrick_prescott`.
        t: int, optional
            The time window to scan df.
        n: int, optional
            The number of coefficients in the output.
        detrend: bool, optional
            Whether to use de-trended data or the original temperature data.
    """
    M_matrix = pd.DataFrame(columns=[f"c{i+1}" for i in range(n)])
    for df, uuid in zip(dfs, ids):
        if detrend:
            df = df.loc[:, ["timepoint", "detrended"]]
        else:
            df = df.loc[:, ["timepoint", "temperature"]]
        df["measure_date"] = df.timepoint.apply(lambda x: x.date())
        time_windows = df["measure_date"].unique()
        time_windows = [time_windows[i:i + t] for i in range(len(time_windows) - t + 1)]
        for time_window in time_windows:
            data = df.loc[df["measure_date"].isin(time_window)]
            date_range = f"{str(time_window.min())[-5:]}_{str(time_window.max())[-5:]}"
            if detrend:
                M_matrix.loc[f"{uuid}_{date_range}"] = abs(np.fft.fft(a=data["detrended"], n=n))
            else:
                M_matrix.loc[f"{uuid}_{date_range}"] = abs(np.fft.fft(a=data["temperature"], n=n))
    return M_matrix


if __name__ == "__main__":
    os.chdir("D:/Documents/2018/UGA/courses/BINF 8950 - Systems Biology/capstone")
    df = pd.read_csv("malaria.csv")
    df["timepoint"] = pd.to_datetime(df["timepoint"])
    df["id"] = df[["experiment", "case"]].apply("_".join, axis=1)
    ids = df["id"].unique()  # all monkeys in different experiments
    df["condition"] = np.where(df.temperature >= FEVER_CUTOFF, "fever", "normal")  # mark fevers for later coloring in PCA

    trend_list = []
    for uuid in ids:
        x = df[df["id"] == uuid][["timepoint", "temperature"]]
        trend_list.append(hodrick_prescott(x, lamb=HODRICK_PRESCOTT_LAMBDA))

    # dfs = [tsplot(x, uuid) for x, uuid in zip(trend_list, ids)]

    M_matrix = fourier_transform(trend_list, t=TIME_WINDOW, n=FOURIER_COEF_NUM)

    # M_tilde = (M_matrix - M_matrix.mean()) / M_matrix.std()
    # V_matrix = M_tilde.T.dot(M_tilde)
    pca = PCA(n_components=2, random_state=317)
    c = pca.fit(M_matrix.T)  # take transpose of matrix to get the correct dimensions
    x, y = c.components_  # PC1 and PC2

    obs = M_matrix.index.tolist()
    df["newID"] = df["id"] + [f"_{str(x.date())[-5:]}_{str((x + pd.DateOffset(TIME_WINDOW - 1)).date())[-5:]}" for x in df["timepoint"]]
    fevers = df.loc[df["condition"] == "fever", "newID"]
    data = pd.DataFrame(dict(PC1=x, PC2=y), index=obs)
    data["group"] = np.where(data.index.isin(fevers), "fever", "normal")
    # data
    fig = sns.lmplot(data=data, x="PC1", y="PC2", hue="group", size=10, fit_reg=False)
