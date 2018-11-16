import os
import sys


import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as offline
import seaborn as sns
import statsmodels.api as sm
from sklearn.decomposition import PCA

TEST_METHOD: bool = False
TEST_FOUR_DAYS: bool = False  # cannot be True if TEST_METHOD is True

INTERPOLATE_DATA: bool = True  # cannot be True if TEST_METHOD is True
TIME_WINDOW: int = 4
USE_FOURIER_TRANSFORM: bool = False
FOURIER_COEF_NUM: int = 15
HODRICK_PRESCOTT_LAMBDA: int = 15000


try:
    assert sys.version_info.major == 3
    assert sys.version_info.minor > 5
except AssertionError:
    raise RuntimeError("This code requires Python 3.6+.")


def interpolate_data(df):
    """Interpolate to hourly data so that FFT works correctly.

    Parameters
    ----------
        df: pandas DataFrame
            Index([
                'timepoint', 'elapsed_time', 'temperature',
                'experiment', 'case', 'id', 'condition'])

    Returns
    -------
        A interpolated long DataFrame.
    """
    df_interp = pd.DataFrame(columns=df.columns)
    for uuid in ids:
        df_uuid = df[df["id"] == uuid][["timepoint", "temperature"]]
        timepoint = pd.date_range(
            start=df_uuid.timepoint.min(),
            end=df_uuid.timepoint.max(),
            freq="30min"
        )
        df_return = pd.DataFrame(dict(timepoint=timepoint, id=uuid))
        df_return = df_return.join(df_uuid.set_index("timepoint"), on="timepoint", how="left")
        df_return = df_return.interpolate(method="cubic")
        df_interp = pd.concat([df_interp, df_return], ignore_index=True)
    return df_interp


def hodrick_prescott(df: pd.DataFrame, lamb: int=1600) -> pd.DataFrame:
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


def fourier_transform(dfs, t: int=2, n: int=10, detrend: bool=True) -> pd.DataFrame:
    """Perform a Fourier transform on the de-trended data.

    Parameters
    ----------
        dfs: list of pandas DataFrames
            trend_list: a series of output from the function `hodrick_prescott`.
        t: int, optional
            The time window to scan df.
        n: int, optional
            The number of coefficients in the output.
        detrend: bool, optional
            Whether to use de-trended data or the original temperature data.

    Returns
    -------
        The Fast Fourier Transformed matrix of coefficients.
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


def fourier_series(dfs, t: int=2, n: int=10, detrend: bool=True) -> pd.DataFrame:
    """Perform a Fourier cosine series on the data.

    Parameters
    ----------
        dfs: list of pandas DataFrames
            trend_list: a series of output from the function `hodrick_prescott`.
        t: int, optional
            The time window (days) to scan df.
        n: int, optional
            The number of coefficients in the output.
        detrend: bool, optional
            Whether to use de-trended data or the original temperature data.

    Returns
    -------
        S_matrices: A long DataFrame w/ columns timepoint, id, fcs, and temperature / detrended.
        M_matrix: The Fourier Cosine Transformed matrix.
    """
    def integrator(y, x):
        I = 0
        for i in range(0, len(y) - 2, 2):
            a = x.iloc[i]
            b = x.iloc[i + 2]
            beginning = y.iloc[i]
            middle = y.iloc[i + 1]
            ending = y.iloc[i + 2]
            I += (b - a) / 6 * (beginning + 4 * middle + ending)
        return I

    def calc_series(n, x, y):
        L = x.iloc[-1] - x.iloc[0]
        # S = y.mean()
        S = 1 / L * integrator(y, x)
        c = []
        for i in range(1, n + 1):
            p = y * np.cos(i * np.pi * x / L)
            q = np.cos(i * np.pi * x / L) ** 2
            c.append(integrator(p, x) / integrator(q, x))
            S += c[i - 1] * np.cos(i * np.pi * x / L)  # S should be an array w/ the same len as x
        return dict(S=S, c=c)

    S_matrices = pd.DataFrame()
    M_matrix = pd.DataFrame(columns=[f"c{i+1}" for i in range(n)])
    for df, uuid in zip(dfs, ids):
        df = df.sort_values("timepoint")
        df["measure_date"] = df.timepoint.apply(lambda x: x.date())
        time_windows = df["measure_date"].unique()
        time_windows = [time_windows[i:i + t] for i in range(len(time_windows) - t + 1)]
        for time_window in time_windows:
            data = df.loc[df["measure_date"].isin(time_window)]
            date_range = f"{str(time_window.min())[-5:]}_{str(time_window.max())[-5:]}"
            x = data["timepoint"]
            x = (x - x.min()) / np.timedelta64(1, "h")
            if detrend:
                y = data["detrended"]
            else:
                y = data["temperature"]
            series = calc_series(n, x, y)
            M_matrix.loc[f"{uuid}_{date_range}"] = series["c"]
            S_matrix = pd.DataFrame(dict(
                timepoint=data["timepoint"],
                fcs=series["S"],
                id=f"{uuid}_{date_range}"
            ))
            if detrend:
                S_matrix["detrended"] = y
            else:
                S_matrix["temperature"] = y
            S_matrices = S_matrices.append(S_matrix)
    return S_matrices, M_matrix


def tsplot(df, uuid, trend=True, detrended=True, save_image=False):
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
    fig = {
        'data': [
            {
                'x': df['timepoint'],
                'y': df["temperature"],
                'name': "Temperature",
                'opacity': 0.7,
                'mode': 'lines+markers'
            }
        ],
        'layout': {
            "title": uuid,
            'xaxis': {'title': 'Date'},
            'yaxis': {'title': "Temperature"}
        }
    }
    if trend:
        fig["data"].append(
            {
                'x': df['timepoint'],
                'y': df["trend"],
                'name': "Hodrick-Prescott filter trend",
                'opacity': 0.7
            }
        )
    if detrended:
        fig["data"][0] = {
            'x': df['timepoint'],
            'y': df["detrended"],
            'name': "De-trended",
            "yaxis": "y2",
            'opacity': 0.7
        }
        fig["layout"]["yaxis2"] = {
            "title": "Temperature",
            "overlaying": "y",
            "side": "right"
        }
    if save_image:
        offline.iplot(fig, filename=f"{uuid}", show_link=False, image="png")
    else:
        offline.iplot(fig, filename=f"{uuid}", show_link=False)
    return fig


# Read temperature data:
if TEST_METHOD:
    import matplotlib.pyplot as plt
    np.random.seed(0)
    df = pd.DataFrame()
    for uuid in ["group1_", "group2_", "group3_"]:
        x = pd.date_range('1/1/2018', periods=120, freq='H')
        xe = (x - x.min()) / np.timedelta64(1, "h")
        noise = np.random.normal(0, 1, len(x))
        s = np.std(np.cos(xe) * xe) / np.random.normal()
        y = np.cos(xe) * xe + s * noise
        df = df.append(pd.DataFrame(
            dict(timepoint=x,
                 elapsed_time=xe,
                 temperature=y,
                 id=uuid
                 )))
    plt.plot(x, y)
else:
    df = pd.read_csv("malaria.csv")
    df["timepoint"] = pd.to_datetime(df["timepoint"])
    df["id"] = df[["experiment", "case"]].apply("_".join, axis=1)
df = df.loc[:, ["timepoint", "temperature", "id"]]
ids = df["id"].unique()
df.head()

if TEST_FOUR_DAYS:
    res = pd.DataFrame()
    for uuid in ids:
        x = df[df["id"] == uuid]
        first_four_days = x["timepoint"].min() + pd.DateOffset(3)
        x = x.loc[x.timepoint <= first_four_days]
        res = res.append(x)
    df = res


# Interpolate temperature data and find when the fevers happened:
if INTERPOLATE_DATA:
    df = interpolate_data(df)

# Calculate the trend using the Hodrick-Prescott filter:
trend_list = []
for uuid in ids:
    x = df[df["id"] == uuid][["timepoint", "temperature"]]
    trend_list.append(hodrick_prescott(x, lamb=HODRICK_PRESCOTT_LAMBDA))
figs = [tsplot(x, uuid, trend=True, detrended=False, save_image=False) for x, uuid in zip(trend_list, ids) if uuid in ["E30_RKy15", "E06_RIh16", "E07B_11C166"]]

# Check the de-trended data:
if USE_FOURIER_TRANSFORM:
    M_matrix = fourier_transform(trend_list, t=TIME_WINDOW, n=FOURIER_COEF_NUM, detrend=True)
else:
    S_matrix, M_matrix = fourier_series(trend_list, t=TIME_WINDOW, n=FOURIER_COEF_NUM, detrend=True)
    dfm = S_matrix.loc[S_matrix["id"].str.startswith("E30_RKy15")]
    dfm = dfm.groupby("timepoint").mean().reset_index()
    fig = {
        'data': [
            go.Scatter(
                x=pd.to_datetime(dfm['timepoint']),
                y=dfm['detrended'],
                name="Temperature",
                opacity=0.7
            ),
            go.Scatter(
                x=pd.to_datetime(dfm["timepoint"]),
                y=dfm["fcs"],
                name="Fourier cosine series",
                opacity=0.7
            )
        ],
        'layout': {
            'xaxis': {'title': 'Date'},
            'yaxis': {'title': "Temperature"},
            'title': "E30 RKy15"
        }
    }
    offline.iplot(fig, show_link=False, image="png")
M_matrix.head()

pca = PCA(n_components=2, random_state=0)
c = pca.fit(M_matrix.T)
x, y = c.components_
print(f"Explained variance ratio: {c.explained_variance_ratio_}")
projection = pca.inverse_transform(pca.transform(M_matrix.T))
loss = ((M_matrix.T - projection) ** 2).mean()
loss.head()

obs = M_matrix.index.tolist()
data = pd.DataFrame(dict(PC1=x, PC2=y), index=obs)
data["group"] = [x.split("_")[0] for x in data.index]
fig = sns.lmplot(
    data=data, x="PC1", y="PC2", hue="group", size=10,
    fit_reg=False, scatter_kws={'alpha': 0.7}, markers=["o", "x", "s"]
)
