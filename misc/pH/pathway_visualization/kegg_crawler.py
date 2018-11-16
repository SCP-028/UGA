""" Use KEGG API to get entrez IDs for a certain customized pathway.
REACTANT and PRODUCT should have the same lengths.
Sample output:
{
    "putrescine": {
        "R00018": {
            "equation": "2 C00134 <=> C06366 + C00014",
            "enzyme": [
                "2.5.1.44"
            ]
        },
        "R00670": {
            "equation": "C00077 <=> C00134 + C00011",
            "enzyme": [
                "4953"
            ]
        }
    },
    "spermidine": {
        "R01918": {
            "equation": "C05730 + C00001 <=> C00051 + C00315",
            "enzyme": [
                "3.5.1.78",
                "3.5.1.-"
            ]
        }
    }
}
"""
import re
import json

import requests

# PATHWAYNAME = 'Arginine'
# REACTANT = ["L-arginine", "ADMA", "L-Citrulline", "argininosuccinate",
#             "L-arginine", "agmatine", "putrescine", "spermidine",
#             "spermine", "diacetylspermine", "putrescine", "proline",
#             "1-Pyrroline-5-carboxylate", "L-Glutamate 5-semialdehyde"]
# PRODUCT = ["ADMA", "L-Citrulline", "argininosuccinate", "L-arginine",
#            "agmatine", "putrescine", "spermidine", "spermine",
#            "diacetylspermine", "N8-acetylspermidine", "ornithine",
#            "1-Pyrroline-5-carboxylate", "L-Glutamate 5-semialdehyde", "ornithine"]
PATHWAYNAME = "Tryptophan"
REACTANT = [
    "tryptophan", "formylkynurenine", "kynurenine", "kynurenine",
    "3-Hydroxykynurenine", "3-Hydroxyanthranilate"
]
PRODUCT = [
    "formylkynurenine", "kynurenine", "kynurenic acid", "3-Hydroxykynurenine",
    "3-Hydroxyanthranilate", "Quinolinate"
]
METABOLITE = set(REACTANT + PRODUCT)

SS = requests.Session()
URL = "http://rest.kegg.jp"
compounds = {}
reactions = {}
pathway = {}
reactions_regex = re.compile("REACTION\s+([R\d\s\n]*)\n[A-Z]{2,10}")
strip_regex = re.compile("\n|\s+")
ec_regex = re.compile("ENZYME\s+([^A-Z]+)\n")
equation_regex = re.compile("EQUATION\s+([^ABD-Z]+)\n")
genes_regex = re.compile("GENES\s+HSA:\s([^:]+[\n\s]?)+")

# Get KEGG compounds IDs of the metabolites, and corresponding reaction IDs
for c in METABOLITE:
    r = SS.get(f"{URL}/find/compound/{c}")
    cpd = r.text.split('\t', 1)[0]  # get first returned item
    compounds[c] = cpd
    r = SS.get(f"{URL}/get/{cpd}")
    rxn = reactions_regex.findall(r.text)
    if rxn:
        rxn = rxn[0]
        # store as dict for adding entrez later
        reactions[c] = dict.fromkeys(strip_regex.sub(" ", rxn).split(" "))

# Find enzymes and equations for all reactions
for cpd, rxns in reactions.items():
    for rxn in rxns.keys():
        rxns[rxn] = {
            'equation': None,
            'enzyme': []
        }
        r = SS.get(f"{URL}/get/{rxn}")
        eqn = equation_regex.findall(r.text)
        if eqn:
            eqn = eqn[0]
            rxns[rxn]['equation'] = strip_regex.sub(" ", eqn)
        ec = ec_regex.findall(r.text)
        if ec:
            ec = ec[0]
            ec = strip_regex.sub(",", ec).split(",")
            for enzyme in ec:
                r = SS.get(f"{URL}/get/ec:{enzyme}")
                genes = genes_regex.findall(r.text)
                if genes:
                    genes = genes[0]
                    genes = re.findall("(\d+)\(", genes)
                    rxns[rxn]['enzyme'].append(genes)
                else:
                    rxns[rxn]['enzyme'].append([enzyme])
            rxns[rxn]['enzyme'] = [item for sublist in rxns[rxn]['enzyme'] for item in sublist if item != "///"]
    reactions[cpd] = rxns

with open(f"{PATHWAYNAME}_reactions.json", "w") as f:
    json.dump(reactions, f)
with open(f"{PATHWAYNAME}_compounds.json", "w") as f:
    json.dump(compounds, f)
