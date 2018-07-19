#!python3
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.offline as offline
from scipy.integrate import odeint, solve_ivp

import parameters.kcat as Kcat
import parameters.ki as Ki
import parameters.km as Km
from parameters.const import *

# offline.init_notebook_mode()


def michaelis_menten_model(t, y):
    """Compute the derivatives of substrates (from the acetylation model) at time t.

    Parameters
    ----------
        t: NumPy array
            A sequence of time points for which to solve for the substrates.
        y: list
            Current conditions on the substrates.

    Returns
    -------
        New states of the substrates.

    Constraints
    -----------
        CoA1 = max(CoA_total - AcCoA1, 0)
        CoA2 = max(CoA_total - AcCoA2, 0)
        NADH = max(NAD_total - NAD, 0)
    """
    AcCoA1, Citrate1, Citrate2, Isocitrate, AlphaKG, SucCoA, Succinate, \
        Fumarate, Malate, OXa, AcCoA2, Acetate, NAD = y
    # Ac-CoA1 = PDHA1 + ACSS1 - CS
    dAcCoA1dt = (
        (
            (Kcat.PDHA1 * PDHA1 * Pyruvate * NAD * max(CoA_total - AcCoA1, 0)) /
            (
                (Km.PDHA1_Pyruvate * (1 + (AcCoA1 / Ki.PDHA1_AcCoA)) + Pyruvate) *
                (Km.PDHA1_NAD * (1 + (AcCoA1 / Ki.PDHA1_AcCoA)) + NAD) *
                (Km.PDHA1_CoA * (1 + (AcCoA1 / Ki.PDHA1_AcCoA)) + max(CoA_total - AcCoA1, 0))
            )
        ) +
        (
            (Kcat.ACSS1 * ACSS1 * Acetate * max((CoA_total - AcCoA1), 0)) /
            (
                (Km.ACSS1_Acetate * (1 + (AcCoA1 / Ki.ACSS1_AcCoA)) + Acetate) *
                (Km.ACSS1_CoA * (1 + (AcCoA1 / Ki.ACSS1_AcCoA)) + (CoA_total - AcCoA1))
            )
        ) -
        (
            (Kcat.CS * CS * OXa * AcCoA1) /
            (
                (Km.CS_OXa * (1 + (Citrate1 / Ki.CS_Citrate)) + OXa) *
                (Km.CS_AcCoA * (1 + (Citrate1 / Ki.CS_Citrate)) + AcCoA1)
            )
        )
    )
    # citrate1  =  CS - ACO2 - SLC13A5
    dCitrate1dt = (
        (
            (Kcat.CS * CS * OXa * AcCoA1) /
            (
                (Km.CS_OXa * (1 + (Citrate1 / Ki.CS_Citrate)) + OXa) *
                (Km.CS_AcCoA * (1 + (Citrate1 / Ki.CS_Citrate)) + AcCoA1)
            )
        ) -
        (
            (
                (Kcat.ACO2_1 * ACO2 * Citrate1) /
                (Km.ACO2_Citrate * (1 + (Isocitrate / Ki.ACO2_Isocitrate)) + Citrate1)
            ) -
            (
                (Kcat.ACO2_2 * ACO2 * Isocitrate) /
                (Km.ACO2_Isocitrate * (1 + (Citrate1 / Ki.ACO2_Citrate)) + Isocitrate)
            )
        ) -
        (
            (Kcat.SLC13A5 * SLC13A5 * Citrate1) /
            (Km.SLC13A5_Citrate * (1 + (Citrate2 / Ki.SLC13A5_Citrate)) + Citrate1)
        )
    )
    # citrate2 = SLC13A5 - ACLY
    dCitrate2dt = (
        (
            (Kcat.SLC13A5 * SLC13A5 * Citrate1) /
            (Km.SLC13A5_Citrate * (1 + (Citrate2 / Ki.SLC13A5_Citrate)) + Citrate1)
        ) -
        (
            (Kcat.ACLY * ACLY * Citrate2 * max(CoA_total - AcCoA2, 0)) /
            (
                (Km.ACLY_Citrate * (1 + (AcCoA2 / Ki.ACLY_AcCoA)) + Citrate2) *
                (Km.ACLY_CoA * (1 + (AcCoA2 / Ki.ACLY_AcCoA)) + max(CoA_total - AcCoA2, 0))
            )
        )
    )
    # isocitrate  =  ACO2 - IDH2
    dIsocitratedt = (
        (
            (
                (Kcat.ACO2_1 * ACO2 * Citrate1) /
                (Km.ACO2_Citrate * (1 + (Isocitrate / Ki.ACO2_Isocitrate)) + Citrate1)
            ) -
            (
                (Kcat.ACO2_2 * ACO2 * Isocitrate) /
                (Km.ACO2_Isocitrate * (1 + (Citrate1 / Ki.ACO2_Citrate)) + Isocitrate)
            )
        ) -
        (
            (Kcat.IDH2 * IDH2 * Isocitrate * NAD) /
            (
                (Km.IDH2_Isocitrate * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + Isocitrate) *
                (Km.IDH2_NAD * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + NAD)
            )
        )
    )
    # alpha-KG  =  IDH2 - OGDH
    dAlphaKGdt = (
        (
            (Kcat.IDH2 * IDH2 * Isocitrate * NAD) /
            (
                (Km.IDH2_Isocitrate * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + Isocitrate) *
                (Km.IDH2_NAD * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + NAD)
            )
        ) -
        (
            (Kcat.OGDH * OGDH * AlphaKG * NAD) /
            (
                (Km.OGDH_AlphaKG * (1 + (SucCoA / Ki.OGDH_SucCoA)) + AlphaKG) *
                (Km.OGDH_NAD * (1 + (SucCoA / Ki.OGDH_SucCoA)) + NAD)
            )
        )
    )
    # Succinyl-CoA = OGDH - SUCLG2 (missing GDP and GTP)
    dSucCoAdt = (
        (
            (Kcat.IDH2 * IDH2 * Isocitrate * NAD) /
            (
                (Km.IDH2_Isocitrate * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + Isocitrate) *
                (Km.IDH2_NAD * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + NAD)
            )
        ) -
        (
            (Kcat.SUCLG2 * SUCLG2 * SucCoA) /
            (Km.SUCLG2_SucCoA * (1 + (Succinate / Ki.SUCLG2_Succinate)) + SucCoA)
        )
    )
    # Succinate = SUCLG2 - SDHA
    dSuccinatedt = (
        (
            (Kcat.SUCLG2 * SUCLG2 * SucCoA) /
            (Km.SUCLG2_SucCoA * (1 + (Succinate / Ki.SUCLG2_Succinate)) + SucCoA)
        ) -
        (
            (Kcat.SDHA * SDHA * Succinate) /
            (Km.SDHA_Succinate * (1 + (Fumarate / Ki.SDHA_Fumarate)) + Succinate)
        )
    )
    # Fumarate = SDHA - FH
    dFumaratedt = (
        (
            (Kcat.SDHA * SDHA * Succinate) /
            (Km.SDHA_Succinate * (1 + (Fumarate / Ki.SDHA_Fumarate)) + Succinate)
        ) -
        (
            (Kcat.FH * FH * Fumarate) /
            (Km.FH_Fumarate * (1 + (Malate / Ki.FH_Malate)) + Fumarate)
        )
    )
    # Malate = FH - MDH2
    dMalatedt = (
        (
            (Kcat.FH * FH * Fumarate) /
            (Km.FH_Fumarate * (1 + (Malate / Ki.FH_Malate)) + Fumarate)
        ) -
        (
            (Kcat.MDH2 * MDH2 * Malate * NAD) /
            (
                (Km.MDH2_Malate * (1 + (OXa / Ki.MDH2_OXa)) + Malate) *
                (Km.MDH2_NAD * (1 + (OXa / Ki.MDH2_OXa)) + NAD)
            )
        )
    )
    # OXa  =  MDH2 + PC - CS
    dOXadt = (
        (
            (Kcat.MDH2 * MDH2 * Malate * NAD) /
            (
                (Km.MDH2_Malate * (1 + (OXa / Ki.MDH2_OXa)) + Malate) *
                (Km.MDH2_NAD * (1 + (OXa / Ki.MDH2_OXa)) + NAD)
            )
        ) +
        (
            (Kcat.PC * PC * Pyruvate * HCO3) /
            (
                (Km.PC_Pyruvate * (1 + (OXa / Ki.PC_OXa)) + Pyruvate) *
                (Km.PC_HCO3 * (1 + (OXa / Ki.PC_OXa)) + HCO3)
            )
        ) -
        (
            (Kcat.CS * CS * OXa * AcCoA1) /
            (
                (Km.CS_OXa * (1 + (Citrate1 / Ki.CS_Citrate)) + OXa) *
                (Km.CS_AcCoA * (1 + (Citrate1 / Ki.CS_Citrate)) + AcCoA1)
            )
        )
    )
    # Ac-CoA2  =  ACLY + ACSS2 - ACOT12 - FASN - HMGCS1 - KAT2A - KAT2B - EP300
    # No product inhibition for reactions going into end products / pools.
    dAcCoA2dt = (
        (
            (Kcat.ACLY * ACLY * Citrate2 * max(CoA_total - AcCoA2, 0)) /
            (
                (Km.ACLY_Citrate * (1 + (AcCoA2 / Ki.ACLY_AcCoA)) + Citrate2) *
                (Km.ACLY_CoA * (1 + (AcCoA2 / Ki.ACLY_AcCoA)) + max(CoA_total - AcCoA2, 0))
            )
        ) +
        (
            (Kcat.ACSS2 * ACSS2 * Acetate * max(CoA_total - AcCoA2, 0)) /
            (
                (Km.ACSS2_Acetate * (1 + (AcCoA2 / Ki.ACSS2_AcCoA)) + Acetate) *
                (Km.ACSS2_CoA * (1 + (AcCoA2 / Ki.ACSS2_AcCoA)) + max(CoA_total - AcCoA2, 0))
            )
        ) -
        (
            (Kcat.ACOT12 * ACOT12 * AcCoA2) /
            (Km.ACOT12_AcCoA * (1 + (Acetate / Ki.ACOT12_Acetate)) * (1 + (max(CoA_total - AcCoA2, 0) / Ki.ACOT12_CoA)) + AcCoA2)
        ) -
        (
            2 * (Kcat.FASN * FASN * AcCoA2 ** 2 * HCO3 * max(NAD_total - NAD, 0) ** 2) /
            (
                (Km.FASN_AcCoA + AcCoA2) ** 2 *
                (Km.FASN_HCO3 + HCO3) *
                (Km.FASN_NADPH + max(NAD_total - NAD, 0)) ** 2
            )
        ) -
        (
            3 * (Kcat.HMGCS1 * HMGCS1 * AcCoA2 ** 3) /
            (
                (Km.ACAT2_AcCoA + AcCoA2 ** 2) *
                (Km.HMGCS1_AcCoA + AcCoA2)
            )
        ) -
        (
            (Kcat.KAT2A * KAT2A * AcCoA2) /
            (Km.KAT2A_AcCoA + AcCoA2)
        ) -
        (
            (Kcat.KAT2B * KAT2B * AcCoA2) /
            (Km.KAT2B_AcCoA + AcCoA2)
        ) -
        (
            (Kcat.EP300 * EP300 * AcCoA2) /
            (Km.EP300_AcCoA + AcCoA2)
        )
    )

    # acetate  =  SLC16A3 + HDAC1 + HDAC2 + HDAC3 + ACOT12 - ACSS1 - ACSS2
    dAcetatedt = (
        (
            Vmax_SLC16A3 *
            (
                (Acetate_blood / (Acetate_blood + Kt_Acetate)) -
                (Acetate / Acetate + Kt_Acetate)
            )
        ) +
        (Kcat.HDAC1 * HDAC1 / (Acetate / Ki.HDAC1_Acetate)) +
        (Kcat.HDAC2 * HDAC2 / (Acetate / Ki.HDAC2_Acetate)) +
        (Kcat.HDAC3 * HDAC3 / (Acetate / Ki.HDAC3_Acetate)) +
        (
            (Kcat.ACOT12 * ACOT12 * AcCoA2) /
            (Km.ACOT12_AcCoA * (1 + (Acetate / Ki.ACOT12_Acetate)) * (1 + (max(CoA_total - AcCoA2, 0) / Ki.ACOT12_CoA)) + AcCoA2)
        ) -
        (
            (Kcat.ACSS1 * ACSS1 * Acetate * max(CoA_total - AcCoA1, 0)) /
            (
                (Km.ACSS1_Acetate * (1 + (AcCoA1 / Ki.ACSS1_AcCoA)) + Acetate) *
                (Km.ACSS1_CoA * (1 + (AcCoA1 / Ki.ACSS1_AcCoA)) + max(CoA_total - AcCoA1, 0))
            )
        ) -
        (
            (Kcat.ACSS2 * ACSS2 * Acetate * max(CoA_total - AcCoA2, 0)) /
            (
                (Km.ACSS2_Acetate * (1 + (AcCoA2 / Ki.ACSS2_AcCoA)) + Acetate) *
                (Km.ACSS2_CoA * (1 + (AcCoA2 / Ki.ACSS2_AcCoA)) + max(CoA_total - AcCoA2, 0))
            )
        )
    )
    # NAD  =  v_NAD - PDHA1 - IDH2 - OGDH - MDH2
    dNADdt = (
        (
            (V_ATP_ss / n_ATP_NAD) *
            (max(NAD_total - NAD, 0) / NADH_ss)
        ) -
        (
            (Kcat.PDHA1 * PDHA1 * Pyruvate * NAD * max(CoA_total - AcCoA1, 0)) /
            (
                (Km.PDHA1_Pyruvate * (1 + (AcCoA1 / Ki.PDHA1_AcCoA)) + Pyruvate) *
                (Km.PDHA1_NAD * (1 + (AcCoA1 / Ki.PDHA1_AcCoA)) + NAD) *
                (Km.PDHA1_CoA * (1 + (AcCoA1 / Ki.PDHA1_AcCoA)) + max(CoA_total - AcCoA1, 0))
            )
        ) -
        (
            (Kcat.IDH2 * IDH2 * Isocitrate * NAD) /
            (
                (Km.IDH2_Isocitrate * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + Isocitrate) *
                (Km.IDH2_NAD * (1 + (AlphaKG / Ki.IDH2_AlphaKG)) + NAD)
            )
        ) -
        (
            (Kcat.OGDH * OGDH * AlphaKG * NAD) /
            (
                (Km.OGDH_AlphaKG * (1 + SucCoA / Ki.OGDH_SucCoA) + AlphaKG) *
                (Km.OGDH_NAD * (1 + SucCoA / Ki.OGDH_SucCoA) + NAD)
            )
        ) -
        (
            (Kcat.MDH2 * MDH2 * Malate) /
            (Km.MDH2_Malate * (1 + (OXa / Ki.MDH2_OXa)) + Malate)
        )
    )

    dxdt = [
        dAcCoA1dt, dCitrate1dt, dCitrate2dt, dIsocitratedt, dAlphaKGdt,
        dSucCoAdt, dSuccinatedt, dFumaratedt, dMalatedt, dOXadt, dAcCoA2dt,
        dAcetatedt, dNADdt
    ]
    return dxdt


if __name__ == "__main__":
    substrates = [
        "AcCoA1", "Citrate1", "Citrate2", "Isocitrate", "AlphaKG",
        "SucCoA", "Succinate", "Fumarate", "Malate", "OXa", "AcCoA2",
        "Acetate", "NAD"
    ]
    init_state = [1e-30 for _ in substrates]
    t = (0.0, 200.0)
    result = solve_ivp(michaelis_menten_model, y0=init_state, t_span=t)
    df = pd.DataFrame(result.y, index=substrates).T
    df["t"] = result.t
    data = []
    for st in substrates:
        data.append({
            "x": df["t"],
            "y": df[st],
            "name": st,
            "opacity": 0.7,
            "mode": "lines"
        })
    fig = {
        "data": data,
        "layout": {
            "title": "Histone Acetylation Pathways",
            "xaxis": {"title": "Time Step"},
            "yaxis": {"title": "Substrate"}
        }
    }
    offline.iplot(fig, filename="histone_acetylation", show_link=False)
