#!python3
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def model(y, t):
    """Function that returns dx/dt."""
    Km_PDHA1_Pyruvate = 64.8
    Km_PDHA1_NAD = 33
    Kcat_PDHA1 = 69

    Km_ACSS1_Acetate = 73
    Km_ACSS1_CoA = 11
    Kcat_ACSS1 = 1.9

    Km_ACOT13_AcCoA1 = 47
    Kcat_ACOT13 = 1.48

    Km_CS_OXa = 5.9
    Km_CS_AcCoA1 = 5
    Kcat_CS = 167

    Km_ACO2_Citrate = 480
    Km_ACO2_Isocitrate = 120
    Kcat_ACO2_1 = 5.3
    Kcat_ACO2_2 = 1.1

    Km_IDH2_Isocitrate = 2400
    Km_IDH2_NAD = 80
    Kcat_IDH2 = 30

    Km_OGDH_AlphaKG = 4000
    Km_OGDH_NAD = 80
    Kcat_OGDH = 30

    Km_ACLY_Citrate = 78
    Kcat_ACLY = 2.2

    Km_ACSS2_Acetate = 73
    Km_ACSS2_CoA = 11
    Kcat_ACSS2 = 1.9

    Km_ACOT12_AcCoA2 = 47
    Kcat_ACOT12 = 1.48

    Km_ACACA_AcCoA2 = 34
    Km_ACACA_HCO3 = 2100
    Kcat_ACACA = 10.1

    Km_HMGCS1_AcCoA2 = 14
    Kcat_HMGCS1 = 0.041

    Kcat_KAT2A = 0.011
    Km_KAT2A_AcCoA2 = 6.7

    # Km_HDAC1_Ac-protein = 25,
    Kcat_HDAC1 = 2.8
    # Km_HDAC2_Ac-rotein = 32,
    Kcat_HDAC2 = 2
    # Km_HDAC3_Ac-protein = 35,
    Kcat_HDAC3 = 1.5

    Km_PC_Pyruvate = 220
    Km_PC_HCO3 = 3000
    Kcat_PC = 60
    # CONSTANTS
    Pyruvate = 77
    HCO3 = 11200
    # expression from LIHC
    ACACA = 1.3937241
    ACLY = 5.9057259
    ACO2 = 12.529593
    ACOT12 = 21.4535639
    ACOT13 = 11.3360529
    ACSS1 = 0.7508832
    ACSS2 = 14.7159301
    CS = 5.3787853
    HDAC1 = 9.6725321
    HDAC2 = 1.3847293
    HDAC3 = 5.3305112
    HMGCS1 = 26.2653794
    IDH2 = 92.8854856
    KAT2A = 3.4616622
    OGDH = 13.7951864
    PC = 48.5013029
    PDHA1 = 13.2815118
    SLC16A3 = 0.5819519
    # Ac-CoA1  =  PDHA1 + ACSS1 - ACOT13 - CS
    dAcCoA1dt = ((Kcat_PDHA1 * PDHA1 * Pyruvate * NAD) / ((Km_PDHA1_Pyruvate + Pyruvate) * (Km_PDHA1_NAD + NAD))) + ((Kcat_ACSS1 * ACSS1 * Acetate * (15 - AcCoA1)) / ((Km_ACSS1_Acetate + Acetate) * (Km_ACSS1_CoA + (15 - AcCoA1)))) - ((Kcat_ACOT13 * ACOT13 * AcCoA1) / (Km_ACOT13_AcCoA1 + AcCoA1)) - ((Kcat_CS * CS * OXa * AcCoA1) / ((Km_CS_OXa + OXa) * (Km_CS_AcCoA1 + AcCoA1)))
    # citrate  =  CS - ACO2 - ACLY
    dCitratedt = ((Kcat_CS * CS * OXa * AcCoA1) / ((Km_CS_OXa + OXa) * (Km_CS_AcCoA1 + AcCoA1))) - ((Kcat_ACO2_1 * ACO2 * Citrate) / (Km_ACO2_Citrate + Citrate) - ((Kcat_ACO2_2 * ACO2 * Isocitrate) / (Km_ACO2_Isocitrate + Isocitrate))) - ((Kcat_ACLY * ACLY * Citrate) / (Km_ACLY_Citrate + Citrate))
    # isocitrate  =  ACO2 - IDH2
    dIsocitratedt = (((Kcat_ACO2_1 * ACO2 * Citrate) / (Km_ACO2_Citrate + Citrate)) - ((Kcat_ACO2_2 * ACO2 * Isocitrate) / (Km_ACO2_Isocitrate + Isocitrate))) - ((Kcat_IDH2 * IDH2 * Isocitrate * NAD) / ((Km_IDH2_Isocitrate + Isocitrate) * (Km_IDH2_NAD + NAD)))
    # alpha-KG  =  IDH2 - OGDH
    dAlphaKGdt = ((Kcat_IDH2 * IDH2 * Isocitrate * NAD) / ((Km_IDH2_Isocitrate + Isocitrate) * (Km_IDH2_NAD + NAD))) - ((Kcat_OGDH * OGDH * AlphaKG * NAD) / ((Km_OGDH_AlphaKG + AlphaKG) * (Km_OGDH_NAD + NAD)))
    # OXa  =  OGDH + PC - CS
    dOXadt = ((Kcat_OGDH * OGDH * AlphaKG * NAD) / ((Km_OGDH_AlphaKG + AlphaKG) * (Km_OGDH_NAD + NAD))) + ((Kcat_PC * PC * Pyruvate * HCO3) / ((Km_PC_Pyruvate + Pyruvate) * (Km_PC_HCO3 + HCO3))) - ((Kcat_CS * CS * OXa * AcCoA1) / ((Km_CS_OXa + OXa) * (Km_CS_AcCoA1 + AcCoA1)))
    # Ac-CoA2  =  ACLY + ACSS2 - ACACA - HMGCS1 - KAT2A - ACOT12
    dAcCoA2dt = ((Kcat_ACLY * ACLY * Citrate) / (Km_ACLY_Citrate + Citrate)) + ((Kcat_ACSS2 * ACSS2 * Acetate * (15 - AcCoA2)) / ((Km_ACSS2_Acetate + Acetate) * (Km_ACSS2_CoA + (15 - AcCoA2)))) - ((2 * Kcat_ACACA * AcCoA2 * ACACA * HCO3) / ((Km_ACACA_AcCoA2 + AcCoA2) * (Km_ACACA_HCO3 + HCO3))) - ((Kcat_HMGCS1 * HMGCS1 * AcCoA2 * AcCoA2) / ((Km_HMGCS1_AcCoA2 + AcCoA2) * (Km_HMGCS1_AcCoA2 + AcCoA2))) - ((Kcat_KAT2A * KAT2A * AcCoA2) / (Km_KAT2A_AcCoA2 + AcCoA2)) - ((Kcat_ACOT12 * ACOT12 * AcCoA2) / (Km_ACOT12_AcCoA2 + AcCoA2))
    # acetate  =  SLC16A3 + HDAC1 + HDAC2 + HDAC3 + ACOT13 + ACOT12 - ACSS1 - ACSS2
    dAcetatedt = 88.07 + (Kcat_HDAC1 * HDAC1) + (Kcat_HDAC2 * HDAC2) + (Kcat_HDAC3 * HDAC3) + ((Kcat_ACOT13 * ACOT13 * AcCoA1) / (Km_ACOT13_AcCoA1 + AcCoA1)) + ((Kcat_ACOT12 * ACOT12 * AcCoA2) / (Km_ACOT12_AcCoA2 + AcCoA2)) - ((Kcat_ACSS1 * ACSS1 * Acetate * (15 - AcCoA1)) / ((Km_ACSS1_Acetate + Acetate) * (Km_ACSS1_CoA + (15 - AcCoA1)))) - ((Kcat_ACSS2 * ACSS2 * Acetate * (15 - AcCoA2)) / ((Km_ACSS2_Acetate + Acetate) * (Km_ACSS2_CoA + (15 - AcCoA2))))
    # NAD  =  v_NAD - PDHA1 - IDH2 - OGDH
    dNADdt = ((0.20844 - (0.0012 * NAD)) / 85) - ((Kcat_PDHA1 * PDHA1 * Pyruvate * NAD) / ((Km_PDHA1_Pyruvate + Pyruvate) * (Km_PDHA1_NAD + NAD))) - ((Kcat_IDH2 * IDH2 * Isocitrate * NAD) / ((Km_IDH2_Isocitrate + Isocitrate) * (Km_IDH2_NAD + NAD))) - ((Kcat_OGDH * OGDH * AlphaKG * NAD) / ((Km_OGDH_AlphaKG + AlphaKG) * (Km_OGDH_NAD + NAD)))

    dxdt = [
        dAcCoA1dt, dCitratedt, dIsocitratedt, dAlphaKGdt, dOXadt,
        dAcCoA2dt, dAcetatedt, dNADdt
    ]
    return dxdt


# initial condition
x0 = [1e-30, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30, 1e-30]
t = np.linspace(0, 100, 1000)

# solve ODE
y = odeint(model, x0, t)