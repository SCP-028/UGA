from cobra import Model, Reaction, Metabolite
import cobra.io

model = Model("Lysine")

# Metabolites use entry ids in KEGG
C00047 = Metabolite(
    id="Lysine",
    name="L-Lysine"
)
C00449 = Metabolite(
    id="Saccharopine",
    name="Saccharopine"
)
C04076 = Metabolite(
    id="Alpha_AASA",
    name="L-2-Aminoadipate-6-semialdehyde"
)
C00956 = Metabolite(
    id="Aminoadipate",
    name="L-2-Aminoadipate"
)
C00322 = Metabolite(
    id="Oxoadipate",
    name="2-Oxoadipate"
)
C06157 = Metabolite(
    id="Glutaryl_dihydrolipoamide",
    name="S-Glutaryl-dihydrolipoamide"
)
C00527 = Metabolite(
    id="Glutaryl_CoA",
    name="Glutaryl-CoA"
)
C00877 = Metabolite(
    id="Crotonoyl_CoA",
    name="Crotonoyl-CoA"
)
C01144 = Metabolite(
    id="Hydroxy_butanoyl_CoA",
    name="(S)-3-Hydroxy-butanoyl-CoA"
)
C00332 = Metabolite(
    id="Acetoacetyl_CoA",
    name="Acetoacetyl-CoA"
)
C00024 = Metabolite(
    id="Acetyl_CoA",
    name="Acetyl-CoA"
)

# Side metabolites
C00026 = Metabolite(
    id="Oxoglutarate",
    name="2-Oxoglutarate"
)
C00005 = Metabolite(
    id="NADPH",
    name="NADPH"
)
C00080 = Metabolite(
    id="H+",
    name="H+"
)
C00006 = Metabolite(
    id="NADP",
    name="NADP"
)
C00001 = Metabolite(
    id="H2O",
    name="H2O"
)
C00003 = Metabolite(
    id="NAD",
    name="NAD+"
)
C00025 = Metabolite(
    id="Glutamate",
    name="L-Glutamate"
)
C00004 = Metabolite(
    id="NADH",
    name="NADH"
)
C15972 = Metabolite(
    id="Lipoamide_E",
    name="Enzyme N6-(lipoyl)lysine"
)
C00011 = Metabolite(
    id="CO2",
    name="CO2"
)
C15973 = Metabolite(
    id="Dihydrolipoamide_E",
    name="Enzyme N6-(dihydrolipoyl)lysine"
)
C00010 = Metabolite(
    id="CoA",
    name="Coenzyme A"
)
C04253 = Metabolite(
    id="Electron_transferring_flavoprotein",
    name="Electron_transferring_flavoprotein"
)
C04570 = Metabolite(
    id="Reduced_electron_transferring_flavoprotein",
    name="Reduced_electron_transferring_flavoprotein"
)


# Reactions
R00716 = Reaction(
    id="EC1.5.1.8",
    name="N6-(L-1,3-Dicarboxypropyl)-L-lysine:NADP+ oxidoreductase"
)
R00716.add_metabolites({
    C00047: -1.0,
    C00026: -1.0,
    C00005: -1.0,
    C00080: -1.0,
    C00449: 1.0,
    C00006: 1.0,
    C00001: 1.0
}
)
R00716.gene_reaction_rule = "AASS"

R02315 = Reaction(
    id="EC1.5.1.9",
    name="N6-(L-1,3-Dicarboxypropyl)-L-lysine:NADP+ oxidoreductase"
)
R02315.add_metabolites({
    C00449: -1.0,
    C00006: -1.0,
    C00001: -1.0,
    C00025: 1.0,
    C04076: 1.0,
    C00005: 1.0,
    C00080: 1.0
})
R02315.gene_reaction_rule = "AASS"

R03103 = Reaction(
    id="EC1.2.1.31",
    name="L-2-aminoadipate-6-semialdehyde:NADP+ 6-oxidoreductase"
)
R03103.add_metabolites({
    C04076: -1.0,
    C00006: -1.0,
    C00001: -1.0,
    C00956: 1.0,
    C00005: 1.0,
    C00080: 1.0
})
R03103.gene_reaction_rule = "ALDH7A1"

R01939 = Reaction(
    id="EC2.6.1.39",
    name="L-2-aminoadipate:2-oxoglutarate aminotransferase"
)
R01939.add_metabolites({
    C00956: -1.0,
    C00026: -1.0,
    C00322: 1.0,
    C00025: 1.0
})
R01939.gene_reaction_rule = "AADAT"

R01940 = Reaction(
    id="EC1.2.4.2",
    name="2-Oxoadipate:lipoamde 2-oxidoreductase"
)
R01940.add_metabolites({
    C00322: -1.0,
    C15972: -1.0,
    C06157: 1.0,
    C00011: 1.0
})
R01940.gene_reaction_rule = "( OGDH or OGDHL )"

R02571 = Reaction(
    id="EC2.3.1.61",
    name="Glutaryl-CoA:dihydrolipoamide S-succinyltransferase"
)
R02571.add_metabolites({
    C00527: 1.0,
    C15973: 1.0,
    C00010: -1.0,
    C06157: -1.0
})
R02571.gene_reaction_rule = "DLST"

R02488 = Reaction(
    id="EC1.3.8.6",
    name="glutaryl-CoA:electron-transfer flavoprotein 2,3-oxidoreductase "
)
R02488.add_metabolites({
    C00527: -1.0,
    C04253: -1.0,
    C00877: 1.0,
    C04570: 1.0,
    C00011: 1.0
})
R02488.gene_reaction_rule = "GCDH"

R03026 = Reaction(
    id="EC4.2.1.17",
    name="(S)-3-hydroxybutanoyl-CoA hydro-lyase"
)
R03026.add_metabolites({
    C01144: 1.0,
    C00877: -1.0,
    C00001: -1.0
})
R03026.gene_reaction_rule = "( ECHS1 and EHHADH and HADHA )"

R01975 = Reaction(
    id="EC1.1.1.35",
    name="(S)-3-Hydroxybutanoyl-CoA:NAD+ oxidoreductase"
)
R01975.add_metabolites({
    C01144: -1.0,
    C00003: -1.0,
    C00332: 1.0,
    C00004: 1.0,
    C00080: 1.0
})
R01975.gene_reaction_rule = " ( EHHADH and HADH )"

R00238 = Reaction(
    id="EC2.3.1.9",
    name="Acetyl-CoA:acetyl-CoA C-acetyltransferase"
)
R00238.add_metabolites({
    C00024: 2.0,
    C00010: -1.0,
    C00332: -1.0
})
R00238.gene_reaction_rule = "( ACAT1 or ACAT2 )"


model.add_reactions([R00716, R02315, R03103, R01939, R01940,
                     R02571, R02488, R03026, R01975, R00238])

print("Reactions\n---------")
for x in model.reactions:
    print("%s : %s" % (x.id, x.reaction))

print("\nGenes\n-----")
for x in model.genes:
    associated_ids = (i.id for i in x.reactions)
    print("%s is associated with reactions: %s" %
          (x.id, "{" + ", ".join(associated_ids) + "}"))
cobra.io.save_json_model(model, filename="lysine_degradation.json")

# Requirements
# certifi==2018.1.18
# cobra==0.11.2
# future==0.16.0
# mpmath==1.0.0
# numpy==1.14.1
# optlang==1.3.0
# pandas==0.22.0
# python-dateutil==2.6.1
# pytz==2018.3
# ruamel.yaml==0.14.12
# six==1.11.0
# swiglpk==1.4.4
# sympy==1.1.1
# tabulate==0.8.2
# wincertstore==0.2
