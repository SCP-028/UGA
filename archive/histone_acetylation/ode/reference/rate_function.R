v_PDHA1 <- (
    (Kcat_PDHA1 * c_PDHA1 * c_Pyruvate * c_NAD * c_CoA1) /
    (
        (Km_PDHA1_Pyruvate + c_Pyruvate) *
        (Km_PDHA1_NAD + c_NAD) *
        (Km_PDHA1_CoA1 + c_CoA1) *
        (1 + (c_AcCoA1 / Ki))^3
    )
)
v_ACSS1 <- (
    (Kcat_ACSS1 * c_ACSS1 * c_Acetate * c_CoA1) /
    (
        (Km_ACSS1_Acetate + c_Acetate) *
        (Km_ACSS1_CoA1 + c_CoA1) *
        (1 + (c_Citrate1 / Ki))
    )
)
v_CS <- (
    (Kcat_CS * c_CS * c_OXa * c_AcCoA1) /
    (
        (Km_CS_OXa + c_OXa) *
        (Km_CS_AcCoA1 + c_AcCoA1) *
        (1 + (c_Citrate1 / Ki))^2
    )
)
v_ACO2 <- (
    (
        (Kcat_ACO2_1 * c_ACO2 * c_Citrate) /
        (
            (Km_ACO2_Citrate + c_Citrate) *
            (1 + (c_Isocitrate / Ki))
        )
    ) -
    (
        (Kcat_ACO2_2 * c_ACO2 * c_Isocitrate) /
        (Km_ACO2_Isocitrate + c_Isocitrate) *
        (1 + (c_Citrate1 / Ki))
    )
)
v_CTP <- (
    (Kcat_CTP * c_CTP * c_Citrate1) /
    (
        (Km_CTP_Citrate + c_Citrate1) *
        (1 + (c_Citrate2 / Ki))
    )
)
v_ACLY <- (
    (Kcat_ACLY * c_ACLY * c_Citrate2 * c_CoA2 * c_CoA2) /
    (
        (Km_ACLY_Citrate + c_Citrate) *
        (Km_ACLY_CoA2 + c_CoA2) *
        (1 + (c_AcCoA2 / Ki))
    )
)
v_IDH2 <- (
    (Kcat_IDH2 * c_IDH2 * c_Isocitrate * c_NAD) /
    (
        (Km_IDH2_Isocitrate + c_Isocitrate) *
        (Km_IDH2_NAD + c_NAD) *
        (1 + (c_AlphaKG / Ki))^2
    )
)
v_OGDH <- (
    (Kcat_OGDH * c_OGDH * c_AlphaKG * c_NAD) /
    (
        (Km_OGDH_AlphaKG + c_AlphaKG) *
        (Km_OGDH_NAD + c_NAD) *
        (1 + (c_SuccinylCoA / Ki))^2
    )
)
v_SUDG1 <- (
    (Kcat_OGDH * c_SUDG1  * c_SuccinylCoA) /
    (
        (Km_SUDG1_SuccinylCoA + c_SuccinylCoA) *
        (1 + (c_Succinate / Ki))
    )
)
v_SDHA <- (
    (Kcat_SDHA * c_SDHA * c_Succinate) /
    (
        (Km_SDHA_Succinate + c_Succinate) *
        (1 + (c_Fumarate / Ki))
    )
)
v_FH <- (
    (Kcat_FH * c_FH * c_Fumarate) /
    (
        (Km_FH_Fumarate + c_Fumarate) *
        (1 + (c_Malate / Ki))
    )
)
v_MDH2 <- (
    (Kcat_MDH2 * c_MDH2 * c_Malate) /
    (
        (Km_MDH2_Malate + c_Malate) *
        (1 + (c_OXa / Ki))
    )
)
v_PC <- (
    (Kcat_PC * c_PC * c_Pyruvate * c_HCO3) /
    (
        (Km_PC_Pyruvate + c_Pyruvate) *
        (Km_PC_HCO3 + c_HCO3) *
        (1 + (c_OXa / Ki))^2
    )
)
v_ACSS2 <- (
    (Kcat_ACSS2 * c_ACSS2 * c_Acetate * c_CoA2) /
    (
        (Km_ACSS2_Acetate + c_Acetate) *
        (Km_ACSS2_CoA2 + c_CoA2) *
        (1 + (c_AcCoA2 / Ki))^2
    )
)
v_ACOT12 <- (
    (Kcat_ACOT12 * c_ACOT12 * c_AcCoA2) /
    (
        (Km_ACOT12_AcCoA2 + c_AcCoA2) *
        (1 + (c_Acetate / Ki)) *
        (1 + (c_CoA2 / Ki))
    )
)
v_FASN <- (
    (2 * Kcat_FASN * c_FASN * c_AcCoA2^2 * c_HCO3 * c_NADPH^2) /
    (
        (Km_FASN_AcCoA2 + c_AcCoA2)^2 *
        (Km_FASN_HCO3 + c_HCO3) *
        (Km_FASN_NADPH + c_NADPH)^2
    )
)
v_HMGCS1 <- (
    (3 * Kcat_HMGCS1 * c_HMGCS1 * c_AcCoA2^3) /
    (
        (Km_ACAT2_AcCoA2 + c_AcCoA2) ^ 2 *
        (Km_HMGCS1_AcCoA2 + c_AcCoA2)
    )
)
v_KAT2A <- (
    (Kcat_KAT2A * c_KAT2A * c_AcCoA2) /
    (Km_KAT2A_AcCoA2 + c_AcCoA2)
)
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
v_HDAC1 <- (
    (Kcat_HDAC1 * c_HDAC1) /
    (c_Acetate / Ki)
)
v_HDAC2 <- (
    (Kcat_HDAC2 * c_HDAC2) /
    (c_Acetate / Ki)
)
v_HDAC3 <- (
    (Kcat_HDAC3 * c_HDAC3) /
    (c_Acetate / Ki)
)
v_NAD <- (
    (v_ATP_ss / n_ATP_NAD) *
    (c_NADH / c_NADH_ss)
)
