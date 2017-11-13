# ODE solver
library(deSolve)
source("./parameters.R")


michaelis.menten <- function(t, y, parameters) {
    with(as.list(c(y, parameters)), {
        # Ac-CoA1 = PDHA1 + ACSS1 - ACOT13 - CS
        dAcCoA1 <- ((Kcat_PDHA1 * PDHA1 * Pyruvate * NAD) /((Km_PDHA1_Pyruvate + Pyruvate) * (Km_PDHA1_NAD + NAD))) + ((Kcat_ACSS1 * ACSS1 * Acetate * (15 - AcCoA1)) / ((Km_ACSS1_Acetate + Acetate) * (Km_ACSS1_CoA + (15 - AcCoA1)))) - ((Kcat_ACOT13 * ACOT13 * AcCoA1) / (Km_ACOT13_AcCoA1 + AcCoA1)) - ((Kcat_CS * CS * OXa * AcCoA1) / ((Km_CS_OXa + OXa) * (Km_CS_AcCoA1 + AcCoA1)))
        # citrate = CS - ACO2 - ACLY
        dCitrate <- ((Kcat_CS * CS * OXa * AcCoA1) / ((Km_CS_OXa + OXa) * (Km_CS_AcCoA1 + AcCoA1))) - ((Kcat_ACO2_1 * ACO2 * Citrate) / (Km_ACO2_Citrate + Citrate) - ((Kcat_ACO2_2 * ACO2 * Isocitrate) / (Km_ACO2_Isocitrate + Isocitrate))) - ((Kcat_ACLY * ACLY * Citrate) / (Km_ACLY_Citrate + Citrate))
        # isocitrate = ACO2 - IDH2
        dIsocitrate <- (((Kcat_ACO2_1 * ACO2 * Citrate) / (Km_ACO2_Citrate + Citrate)) - ((Kcat_ACO2_2 * ACO2 * Isocitrate) / (Km_ACO2_Isocitrate + Isocitrate))) - ((Kcat_IDH2 * IDH2 * Isocitrate * NAD) / ((Km_IDH2_Isocitrate + Isocitrate) * (Km_IDH2_NAD + NAD)))
        # alpha-KG = IDH2 - OGDH
        dAlphaKG <- ((Kcat_IDH2 * IDH2 * Isocitrate * NAD) / ((Km_IDH2_Isocitrate + Isocitrate) * (Km_IDH2_NAD + NAD))) - ((Kcat_OGDH * OGDH * AlphaKG * NAD) / ((Km_OGDH_AlphaKG + AlphaKG) * (Km_OGDH_NAD + NAD)))
        # OXa = OGDH + PC - CS
        dOXa <- ((Kcat_OGDH * OGDH * AlphaKG * NAD) / ((Km_OGDH_AlphaKG + AlphaKG) * (Km_OGDH_NAD + NAD))) + ((Kcat_PC * PC * Pyruvate * HCO3) / ((Km_PC_Pyruvate + Pyruvate) * (Km_PC_HCO3 + HCO3))) - ((Kcat_CS * CS * OXa * AcCoA1) / ((Km_CS_OXa + OXa) * (Km_CS_AcCoA1 + AcCoA1)))
        # Ac-CoA2 = ACLY + ACSS2 - ACACA - HMGCS1 - KAT2A - ACOT12
        dAcCoA2 <- ((Kcat_ACLY * ACLY * Citrate) / (Km_ACLY_Citrate + Citrate)) + ((Kcat_ACSS2 * ACSS2 * Acetate * (15 - AcCoA2)) / ((Km_ACSS2_Acetate + Acetate) * (Km_ACSS2_CoA + (15 - AcCoA2)))) - ((2 * Kcat_ACACA * AcCoA2 * ACACA * HCO3) / ((Km_ACACA_AcCoA2 + AcCoA2) * (Km_ACACA_HCO3 + HCO3))) - ((Kcat_HMGCS1 * HMGCS1 * AcCoA2 * AcCoA2) / ((Km_HMGCS1_AcCoA2 + AcCoA2) * (Km_HMGCS1_AcCoA2 + AcCoA2))) - ((Kcat_KAT2A * KAT2A * AcCoA2) / (Km_KAT2A_AcCoA2 + AcCoA2)) - ((Kcat_ACOT12 * ACOT12 * AcCoA2) / (Km_ACOT12_AcCoA2 + AcCoA2))
        # acetate = SLC16A3 + HDAC1 + HDAC2 + HDAC3 + ACOT13 + ACOT12 - ACSS1 - ACSS2
        dAcetate <- 88.07 + (Kcat_HDAC1 * HDAC1) + (Kcat_HDAC2 * HDAC2) + (Kcat_HDAC3 * HDAC3) + ((Kcat_ACOT13 * ACOT13 * AcCoA1) / (Km_ACOT13_AcCoA1 + AcCoA1)) + ((Kcat_ACOT12 * ACOT12 * AcCoA2) / (Km_ACOT12_AcCoA2 + AcCoA2)) - ((Kcat_ACSS1 * ACSS1 * Acetate * (15 - AcCoA1)) / ((Km_ACSS1_Acetate + Acetate) * (Km_ACSS1_CoA + (15 - AcCoA1)))) - ((Kcat_ACSS2 * ACSS2 * Acetate * (15 - AcCoA2)) / ((Km_ACSS2_Acetate + Acetate) * (Km_ACSS2_CoA + (15 - AcCoA2))))
        # NAD = v_NAD - PDHA1 - IDH2 - OGDH
        dNAD <- ((0.20844 - (0.0012 * NAD)) / 85) - ((Kcat_PDHA1 * PDHA1 * Pyruvate * NAD) / ((Km_PDHA1_Pyruvate + Pyruvate) * (Km_PDHA1_NAD + NAD))) - ((Kcat_IDH2 * IDH2 * Isocitrate * NAD) / ((Km_IDH2_Isocitrate + Isocitrate) * (Km_IDH2_NAD + NAD))) - ((Kcat_OGDH * OGDH * AlphaKG * NAD) / ((Km_OGDH_AlphaKG + AlphaKG) * (Km_OGDH_NAD + NAD)))

        list(c(dAcCoA1, dCitrate, dIsocitrate, dAlphaKG, dOXa, dAcCoA2, dAcetate, dNAD))
        }
    )
}

times <- seq(0, 1000, by = 0.01)
out <- ode(y = initial_state, parms = parameters,
           func = michaelis.menten, times = times)
