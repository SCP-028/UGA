import matplotlib.pyplot as plt
import seaborn as sns

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

sns.set_style("whitegrid")
sns.set_context("poster")


class kcat:
    def __init__(self):
        self.ACACA = 10.1
        self.ACAT2 = 2.1
        self.ACLY = 2.2
        self.ACO2_1 = 5.3
        self.ACO2_2 = 1.1
        self.ACOT12 = 1.48
        self.ACOT13 = 1.48
        self.ACSS1 = 1.9
        self.ACSS2 = 1.9
        self.CS = 167
        self.EP300 = 0.0383
        self.FASN = 2.7
        self.FH = 51.7
        self.HDAC1 = 2.8
        self.HDAC2 = 2
        self.HDAC3 = 1.5
        self.HMGCS1 = 0.041
        self.IDH2 = 30
        self.KAT2A = 0.012
        self.KAT2B = 0.005
        self.MDH2 = 18
        self.OGDH = 30
        self.PC = 60
        self.PDHA1 = 69
        self.SDHA = 83.3
        self.SLC13A5 = 0.02  # CTP
        self.SUCLG2 = 201


class km:
    def __init__(self):
        self.ACACA_AcCoA = 34
        self.ACACA_HCO3 = 2100
        self.ACAT2_AcCoA = 29
        self.ACLY_Citrate = 78
        self.ACLY_CoA = 14
        self.ACO2_Citrate = 480
        self.ACO2_Isocitrate = 120
        self.ACOT12_AcCoA = 47
        self.ACOT13_AcCoA = 47
        self.ACSS1_Acetate = 73
        self.ACSS1_CoA = 11
        self.ACSS2_Acetate = 73
        self.ACSS2_CoA = 11
        self.CS_AcCoA = 5
        self.CS_OXa = 5.9
        self.EP300_AcCoA = 0.28
        self.FASN_AcCoA = 7
        self.FASN_HCO3 = 2100
        self.FASN_NADPH = 5
        self.FH_Fumarate = 13
        self.HMGCS1_AcCoA = 14
        self.IDH2_Isocitrate = 320  # 2400
        self.IDH2_NAD = 80
        self.KAT2A_AcCoA = 6.7
        self.KAT2B_AcCoA = 1.1
        self.MDH2_Malate = 560
        self.MDH2_NAD = 140
        self.OGDH_AlphaKG = 13  # 4000
        self.OGDH_NAD = 80
        self.PC_Pyruvate = 220
        self.PC_HCO3 = 3200
        self.PDHA1_Pyruvate = 64.8
        self.PDHA1_NAD = 33
        self.PDHA1_CoA = 4
        self.SDHA_Succinate = 20
        self.SDHA_Fumarate = 25
        self.SLC13A5_Citrate = 600  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2913483/
        self.SUCLG2_SucCoA = 86


class ki:
    def __init__(self):
        self.ACLY_AcCoA = 10
        self.ACO2_Isocitrate = 410
        self.ACO2_Citrate = 190
        self.ACOT12_Acetate = 50  # 50 - 200
        self.ACOT12_CoA = 17
        self.ACSS1_AcCoA = 2700
        self.ACSS2_AcCoA = 10
        self.CS_Citrate = 190
        self.SLC13A5_Citrate = 190
        self.FH_Malate = 2500
        self.HDAC1_Acetate = 50  # 50 - 200
        self.HDAC2_Acetate = 50  # 50 - 200
        self.HDAC3_Acetate = 50  # 50 - 200
        self.IDH2_AlphaKG = 0.29
        self.MDH2_OXa = 9.5
        self.OGDH_SucCoA = 4
        self.PC_OXa = 9.5
        self.PDHA1_AcCoA = 20
        self.SDHA_Fumarate = 1300
        self.SUCLG2_Succinate = 14  # SUDG1, value not changed


Km = km()
Kcat = kcat()
Ki = ki()

Acetate_blood = 125
Kt_Acetate = 0.157
CoA_total = 15
HCO3 = 5000
NAD_total = 46.3
Pyruvate = 77
Vmax_SLC16A3 = 88.07
V_ATP_ss = 0.2
n_ATP_NAD = 2.5
NADH_ss = 22
# expression from LIHC
ACACA = 1.3937241
ACLY = 5.9057259
ACO2 = 12.529593
ACOT12 = 21.4535639
ACOT13 = 11.3360529
ACSS1 = 0.7508832
ACSS2 = 14.7159301
CS = 5.3787853
EP300 = 3.5370080
FASN = 26.7645829
FH = 70.3469031
HDAC1 = 9.6725321
HDAC2 = 1.3847293
HDAC3 = 5.3305112
HMGCS1 = 26.2653794
IDH2 = 92.8854856
KAT2A = 3.4616622
KAT2B = 5.8884689
MDH2 = 46.9342578
OGDH = 13.7951864
PC = 48.5013029
PDHA1 = 13.2815118
SDHA = 26.5377130
SLC13A5 = 41.3838758
SLC16A3 = 0.5819519
SUCLG2 = 54.2205369


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
    AcCoA1, Citrate1, Citrate2, Isocitrate, AlphaKG, SucCoA, Succinate, Fumarate, Malate, OXa, AcCoA2, Acetate, NAD = y
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
            (Kcat.MDH2 * MDH2 * Malate * NAD) /
            (Km.MDH2_Malate * (1 + (OXa / Ki.MDH2_OXa)) + Malate)
        )
    )

    dxdt = [
        dAcCoA1dt, dCitrate1dt, dCitrate2dt, dIsocitratedt, dAlphaKGdt,
        dSucCoAdt, dSuccinatedt, dFumaratedt, dMalatedt, dOXadt, dAcCoA2dt,
        dAcetatedt, dNADdt
    ]
    return dxdt


# Solving ODE
substrates = [
    "AcCoA1", "Citrate1", "Citrate2", "Isocitrate", "AlphaKG",
    "SucCoA", "Succinate", "Fumarate", "Malate", "OXa", "AcCoA2",
    "Acetate", "NAD"
]
init_state = pd.Series(
    [
        10, 132, 132, 34, 7.0,
        10, 20, 150, 5.22, 61, 10,
        125, 24.3
    ], index=substrates)
t = (0.0, 300.0)
result = solve_ivp(michaelis_menten_model, y0=init_state, t_span=t)
df = pd.DataFrame(result.y, index=substrates).T
df["t"] = result.t


# Plots for each substrate
for substrate in substrates:
    tmp = df[["t", substrate]]
    plt.figure(figsize=(16, 8))
    ax = sns.lineplot(x="t", y=substrate, data=tmp)
