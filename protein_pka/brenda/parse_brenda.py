def parse_ID(seq):
    """Return the EC number inside the input string.

    Example
    -------
        "ID\t1.1.1.1\n*****************\n" -> "1.1.1.1"
    """
    seq = seq.split("\n")[0]
    return seq[2:].strip()


def parse_hash_id(ids):
    """Return a list of protein numbers within the given string.

    Example
    -------
        "#10,41,43,150#" -> ["10", "41", "43", "150"]
    """
    ans = [seq.strip("#").split(",") for seq in ids]
    ans = [x for item in ans for x in item]
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
    seq = parse_hash_id(seq)
    return seq


with open("/data/annotation/brenda/brenda_download.txt", "r", encoding="utf-8") as f:
    annot = [line.strip() for line in f.readlines()]
# ECs are separated with "///"
annotation = "\n".join(annot).split("\n\n///\n")
for enzyme in annotation:
    enzyme = "".join(enzyme).split("\n\n")
    info = [x for x in enzyme if x.startswith(("ID", "PROTEIN", "PH_OPTIMUM", "PH_RANGE"))]

enzyme = annotation[0]
enzyme = "".join(enzyme).split("\n\n")
info = [x for x in enzyme if x.startswith(("ID", "PROTEIN", "PH_OPTIMUM", "PH_RANGE"))]
# TODO: figure out how to filter the proteins with missing values
