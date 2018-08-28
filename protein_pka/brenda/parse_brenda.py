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


def parse_ph_optimum(seq, human):
    """Extract the pH optimum values of human proteins.

    Parameters
    ----------
        seq: str
        human: list[str] or set(str)

    Returns
    -------
        list[list[str], str] human protein IDs and the pH optimum.
        Note that the "optimum" could still be a range.
    """
    human = set(human)
    seq = seq.split("\nPHO\t")
    seq = [x.split(" ")[:2] for x in seq]  # protein hash ID & pH
    seq = seq[1:]  # remove the "PH_OPTIMUM at the beginning"
    seq = [[parse_hash_id(x[0]), x[1]] for x in seq]
    seq = [x for x in seq if set(x[0]).intersection(human)]
    return seq


with open("/data/annotation/brenda/brenda_download.txt", "r", encoding="utf-8") as f:
    annot = [line.strip() for line in f.readlines()]
# ECs are separated with "///"
annotation = "\n".join(annot).split("\n\n///\n")
ans = {}
for enzyme in annotation:
    enzyme = "".join(enzyme).split("\n\n")
    info = [x for x in enzyme if x.startswith(("ID", "PROTEIN", "PH_OPTIMUM", "PH_RANGE"))]
    if (len(info) == 4):
        uuid = parse_ID(info[0])
        protein = parse_protein(info[1])
        pH_optimum = parse_ph_optimum(info[2], protein)
        pH_range = parse_ph_range(info[3], protein)
        ans[uuid] = {
            "protein": protein,
            "pH_optimum": pH_optimum,
            "pH_range": pH_range
        }

# TODO: figure out how to filter the proteins with missing values
