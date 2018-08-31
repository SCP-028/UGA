import re

import numpy as np
import pandas as pd


def parse_ID(seq):
    """Return the EC number inside the input string.

    Example
    -------
        "ID\t1.1.1.1\n*****************\n" -> "1.1.1.1"
    """
    seq = seq.split("\n")[0]
    return seq[2:].strip()


def parse_hash_id(seq):
    """Return a list of protein numbers within the given string.

    Example
    -------
        "#10,41,43,150#" -> ["10", "41", "43", "150"]
    """
    ans = seq.strip("#").split(",")
    return ans


def parse_protein(seq):
    """Return a list of protein numbers that are human proteins.

    Example
    -------
        "PROTEIN\nPR\t#8# Homo sapiens\n<10,11,273>\nPR\t#9# Rattus norvegicus   <10,49,50,51,142,167,224,228>"
        -> ["8"]
    """
    seq = seq.split("\nPR\t")  # PR is the identifier for proteins
    seq = [x.split(" ")[0] for x in seq if "Homo sapiens" in x]  # extract protein hash ID
    seq = [parse_hash_id(x) for x in seq]
    seq = [x for item in seq for x in item]
    return seq


def parse_pH(seq, protein, condition=["optimum", "range"]):
    """Extract the pH optimum / range values of human proteins.

    Parameters
    ----------
        seq: str
        protein: set(str)
        condition: str, "optimum" or "range"

    Returns
    -------
        list[list[str], str] human protein IDs and the pH optimum.
        Note that the "optimum" could still be a range.
    """
    if condition == "optimum":
        seq = seq.split("\nPHO\t")
    elif condition == "range":
        seq = seq.split("\nPHR\t")
    else:
        raise ValueError("The parameter <condition> must be 'optimum' or 'range'.")
    seq = [x.split(" ", maxsplit=2) for x in seq]  # protein hash ID, pH & comment
    seq = seq[1:]  # remove the "PH_OPTIMUM at the beginning"
    seq = [[parse_hash_id(x[0]), x[1], x[2]] for x in seq]
    seq = [x for x in seq if set(x[0]).intersection(protein)]
    return seq


def populate_pH_optimum(seq, uuid, protein):
    """Generate a [pandas.DataFrame] with pH optimum information.

    Parameters
    ----------
        seq: list[list[str], str]
            Output from `parse_pH_optimum`.
        uuid: str
        protein: set(str)
    """
    ans = pd.DataFrame(columns=["uuid", "protein", "pH", "activity", "pH_range_comment"])
    for entry in seq:
        entry[0] = set(entry[0]).intersection(protein)
        for p in entry[0]:
            row = pd.Series({
                "uuid": uuid,
                "protein": p,
                "pH": entry[1],
                "activity": 1,
                "pH_range_comment": entry[2]
            })
            ans = ans.append(row, ignore_index=True)
    return ans


def populate_pH_range(seq, uuid, protein):
    # TODO: extract information from comment
    ans = pd.DataFrame(columns=["uuid", "protein", "pH", "activity", "pH_range_comment"])
    for entry in seq:
        entry[0] = set(entry[0]).intersection(protein)
        for p in entry[0]:
            row = pd.Series({
                "uuid": uuid,
                "protein": p,
                "pH": entry[1],
                "activity": np.nan,
                "pH_range_comment": entry[2]
            })
            ans = ans.append(row, ignore_index=True)
    return ans


with open("/data/annotation/brenda/brenda_download.txt", "r", encoding="utf-8") as f:
    annot = [line.strip() for line in f.readlines()]
# ECs are separated with "///"
annotation = "\n".join(annot).split("\n\n///\n")
ans = pd.DataFrame(columns=["uuid", "protein", "pH", "activity", "pH_range_comment"])
for i, enzyme in enumerate(annotation):
    enzyme = "".join(enzyme).split("\n\n")
    info = [x for x in enzyme if x.startswith(("ID", "PROTEIN", "PH_OPTIMUM", "PH_RANGE"))]
    # TODO: figure out how to filter the proteins with missing values. Approach now is too strict.
    if (len(info) == 4):
        uuid = parse_ID(info[0])  # TODO: "()" not in uuid
        protein = parse_protein(info[1])
        pH_optimum = parse_pH(info[2], protein, condition="optimum")
        pH_range = parse_pH(info[3], protein, condition="range")
        ans = ans.append(populate_pH_optimum(pH_optimum, uuid, protein))
        ans = ans.append(populate_pH_range(pH_range, uuid, protein))
