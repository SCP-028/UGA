v_PDHA1 <- (
    (Kcat_PDHA1 * c_PDHA1 * c_Pyruvate * c_NAD) /
    (
        (Km_PDHA1_Pyruvate + c_Pyruvate) *
        (Km_PDHA1_NAD + c_NAD) *
        (Km_PDHA1_CoA1 + c_CoA1)
    )
)
v_ACSS1 <- (
    (Kcat_ACSS1 * c_ACSS1 * c_Acetate * c_CoA1) /
    (
        (Km_ACSS1_Acetate + c_Acetate) *
        (Km_ACSS1_CoA1 + c_CoA1)
    )
)
v_CS <- (
    (Kcat_CS * c_CS * c_OXa * c_AcCoA1) /
    (
        (Km_CS_OXa + c_OXa) *
        (Km_CS_AcCoA1 + c_AcCoA1)
    )
)
v_ACO2 <- (
    (
        (Kcat_ACO2_1 * c_ACO2 * c_Citrate) /
        (Km_ACO2_Citrate + c_Citrate)
    ) -
    (
        (Kcat_ACO2_2 * c_ACO2 * c_Isocitrate) /
        (Km_ACO2_Isocitrate + c_Isocitrate)
    )
)
v_ACLY <- (
    (Kcat_ACLY * c_ACLY * c_Citrate * c_CoA2) /
    (
        (Km_ACLY_Citrate + c_Citrate) *
        (Km_ACLY_CoA2 + c_CoA2)
    )
)
v_IDH2 <- (
    (Kcat_IDH2 * c_IDH2 * c_Isocitrate * c_NAD) /
    (
        (Km_IDH2_Isocitrate + c_Isocitrate) *
        (Km_IDH2_NAD + c_NAD)
    )
)
v_OGDH <- (
    (Kcat_OGDH * c_OGDH * c_AlphaKG * c_NAD) /
    (
        (Km_OGDH_AlphaKG + c_AlphaKG) *
        (Km_OGDH_NAD + c_NAD)
    )
)
v_PC <- (
    (Kcat_PC * c_PC * c_Pyruvate * c_HCO3) /
    (
        (Km_PC_Pyruvate + c_Pyruvate) *
        (Km_PC_HCO3 + c_HCO3)
    )
)
v_ACSS2 <- (
    (Kcat_ACSS2 * c_ACSS2 * c_Acetate * c_CoA2) /
    (
        (Km_ACSS2_Acetate + c_Acetate) *
        (Km_ACSS2_CoA2 + c_CoA2)
    )
)
v_ACACA <- (
    (Kcat_ACACA * c_ACACA * c_AcCoA2 * c_HCO3) /
    (
        (Km_ACACA_AcCoA2 + c_AcCoA2) *
        (Km_ACACA_HCO3 + c_HCO3)
    )
)
v_FASN <- (
    (2 * Kcat_FASN * c_FASN * c_AcCoA2 ^ 2 * c_HCO3 * c_NADPH ^ 2) /
    (
        (Km_FASN_AcCoA2 + c_AcCoA2) ^ 2 *
        (Km_FASN_HCO3 + c_HCO3) *
        (Km_FASN_NADPH + c_NADPH) ^ 2
    )
)
v_HMGCS1 <- (
    (3 * Kcat_HMGCS1 * c_HMGCS1 * c_AcCoA2 ^ 3) /
    (
        (Km_ACAT2_AcCoA2 + c_AcCoA2) ^ 2 *
        (Km_HMGCS1_AcCoA2 + c_AcCoA2)
    )
)
v_KAT2A <- (
    (Kcat_KAT2A * c_KAT2A * c_AcCoA2) / 
    (Km_KAT2A_AcCoA2 + c_AcCoA2)
)
v_ACOT12 <- (
    (Kcat_ACOT12 * c_ACOT12 * c_AcCoA2) /
    (Km_ACOT12_AcCoA2 + c_AcCoA2)
)
v_HDAC1 <- (Kcat_HDAC1 * c_HDAC1)
v_HDAC2 <- (Kcat_HDAC2 * c_HDAC2)
v_HDAC3 <- (Kcat_HDAC3 * c_HDAC3)

v_SLC16A3 <- (
    v_max_SLC16A3 *
    (
        (   
            c_Acetate_blood /
            (c_Acetate_blood + k_T_Acetate)
        ) -
        (
            c_Acetate /
            (c_Acetate + k_T_Acetate)
        )
    )
)
v_NAD <- (
    (v_ATP_ss / n_ATP_NAD) *
    (c_NADH / c_NADH_ss)
)
