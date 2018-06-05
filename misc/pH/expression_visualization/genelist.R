load_aa_data <- function() {
    arginine <- c(
        "SLC22A2", "PRMT1", "PRMT2", "PRMT3", "PRMT5", "PRMT6", "PRMT7", "PRMT8", "PRMT9",
        "DDAH1", "PYCR1", "AGMAT", "ASS1", "ASL", "AZIN2", "SRM", "SMS", "RAI1", "SRMS",
        "ALDH18A1", "OAT", "ODC1", "SLC7A1", "SAT1", "SLC26A1", "SMO", "SMOX", "SLC3A2",
        "ALDH4A1", "SLC22A1", "GATM", "SLC22A3", "SLC7A8"
    )
    serine <- c("PHGDH", "PSAT1", "PSPH")
    proline <- c(
        "GPI", "PFKM", "PFKP", "PFKL", "ALDOA", "ALDOB", "ALDOC",
        "TPI1", "GAPDH", "ALDH18A1", "PYCR1"
    )
    tryptophan <- c("TDO2", "IDO1", "IDO2", "AFMID", "NQO1", "KMO", "KYNU", "HAAO")
    lysine <- c(
        "AASS", "ALDH7A1", "AADAT", "OGDH", "OGDHL", "DLST", "GCDH", "ECHS1", "EHHADH",
        "HADHA", "HADH", "ACAT1", "ACAT2"
    )
    glutaminolysis <- c(
        "SLC1A5", "SLC6A14", "SLC6A19", "SLC7A5", "SLC7A6", "SLC7A7", "SLC7A8", "SLC7A9",
        "SLC38A1", "SLC38A2", "SLC38A3", "SLC38A5", "SLC38A7", "SLC38A8", "GLS", "GLS2",
        "GLUD1", "GLUD2", "GPT", "GOT1", "GOT2", "OGDH", "DLD", "DLST", "SUCLG1", "SUCLG2",
        "SUCLA2", "SDHA", "SDHB", "SDHC", "SDHD", "FH", "MDH1", "MDH2", "CS", "ME1", "ME2",
        "ME3", "ACLY", "FASN"
    )
    # https://doi.org/10.1046/j.1440-1746.2000.02205.x - Fig. 1
    bcaa <- c(
        "BCAT1", "BCAT2", "BCKDHA", "BCKDHB", "DBT", "DLD", "IVD", "SLC7A5", "SLC7A6",
        "SLC7A7", "SLC7A8", "SLC43A1", "SLC43A2", "SLC43A3", "SLC16A1", "SLC16A2",
        "SLC16A4", "SLC16A5", "SLC16A6", "SLC16A7"
    )
    genelist <- list(arginine, serine, proline, tryptophan, lysine, glutaminolysis, bcaa)
    names(genelist) <- c(
        "arginine", "serine", "proline", "tryptophan",
        "lysine", "glutaminolysis", "bcaa"
    )
    return(genelist)
}

