v_PDHA1 <- (
    (Kcat_PDHA1 * PDHA1 * Pyruvate * NAD) /
    (
        (Km_PDHA1_Pyruvate + Pyruvate) *
        (Km_PDHA1_NAD + NAD)
    ))
v_ACSS1 <- (
    (Kcat_ACSS1 * ACSS1 * Acetate * (15 - AcCoA1)) /
    (
        (Km_ACSS1_Acetate + Acetate) *
        (Km_ACSS1_CoA + (15 - AcCoA1))
    ))
v_ACOT13 <- (
    (Kcat_ACOT13 * ACOT13 * AcCoA1) /
    (Km_ACOT13_AcCoA1 + AcCoA1))
v_CS <- (
    (Kcat_CS * CS * OXa * AcCoA1) /
    (
        (Km_CS_OXa + OXa) *
        (Km_CS_AcCoA1 + AcCoA1)
    ))
v_ACO2 <- (
    (
        (Kcat_ACO2_1 * ACO2 * Citrate) /
        (Km_ACO2_Citrate + Citrate)
    ) -
    (
        (Kcat_ACO2_2 * ACO2 * Isocitrate) /
        (Km_ACO2_Isocitrate + Isocitrate)
    ))
v_IDH2 <- (
    (Kcat_IDH2 * IDH2 * Isocitrate * NAD) /
    (
        (Km_IDH2_Isocitrate + Isocitrate) *
        (Km_IDH2_NAD + NAD)
    ))
v_OGDH <- (
    (Kcat_OGDH * OGDH * AlphaKG * NAD) /
    (
        (Km_OGDH_AlphaKG + AlphaKG) *
        (Km_OGDH_NAD + NAD)
    ))
v_ACLY <- (
    (Kcat_ACLY * ACLY * Citrate) /
    (Km_ACLY_Citrate + Citrate))
v_ACSS2 <- (
    (Kcat_ACSS2 * ACSS2 * Acetate * (15 - AcCoA2)) /
    (
        (Km_ACSS2_Acetate + Acetate) *
        (Km_ACSS2_CoA + (15 - AcCoA2))
    ))
v_ACOT12 <- (
    (Kcat_ACOT12 * ACOT12 * AcCoA2) /
    (Km_ACOT12_AcCoA2 + AcCoA2))
v_ACACA <- (
    (2 * Kcat_ACACA * AcCoA2 * ACACA * HCO3) /
    (
        (Km_ACACA_AcCoA2 + AcCoA2) *
        (Km_ACACA_HCO3 + HCO3)
    ))
# HMGCS1 unsure
v_HMGCS1 <- (
    (Kcat_HMGCS1 * HMGCS1 * AcCoA2 * AcCoA2) /
    (
        (Km_HMGCS1_AcCoA2 + AcCoA2) *
        (Km_HMGCS1_AcCoA2 + AcCoA2)
    ))
v_KAT2A <- (
    (Kcat_KAT2A * KAT2A * AcCoA2) / 
    (Km_KAT2A_AcCoA2 + AcCoA2))
v_HDAC1 <- (Kcat_HDAC1 * HDAC1)
v_HDAC2 <- (Kcat_HDAC2 * HDAC2)
v_HDAC3 <- (Kcat_HDAC3 * HDAC3)
v_PC <- (
    (Kcat_PC * PC * Pyruvate * HCO3) /
    (
        (Km_PC_Pyruvate + Pyruvate) *
        (Km_PC_HCO3 + HCO3)
    ))
# v_SLC16A3 = Vmax_SLC * (LAC / (LAC + Kt_LAC))
v_SLC16A3 <- 88.07
# v_NAD = (0.4 * v_ATP_ss * NADH) / 85
v_NAD <- ((0.20844 - (0.0012 * NAD)) / 85))
# CoA1 = 15 - AcCoA1
# CoA2 = 15 - AcCoA2
# ACO2 reaction goes both ways