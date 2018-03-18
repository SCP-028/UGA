"""Calculate protein net charge at pH=7."""
import glob
import os
import re
import sys

import pandas as pd


class pka:
    def __init__(self):
        self.acidic = ["ASP", "GLU"]
        self.alkaline = ["ARG", "LYS", "HIS"]
        self.pH = 7.0

    def clean_raw_files(self):
        """Clean data and save in the "cleaned" directory.

        Returns
        -------
            A generator object containing all the cleaned files.
        """
        cleaned_files = glob.iglob(os.path.join(ROOTPATH, "cleaned", "*.pka"))
        cleaned_files = [x[-8:] for x in cleaned_files]
        raw_files = set(glob.glob("*.pka")) - set(cleaned_files)

        strip_regex = re.compile("\s+")
        for file in raw_files:
            df = ["Protein,Residue,Position,pK,slope,1000Chi2"]
            with open(file, "r") as f:
                next(f)  # Skip header line
                for line in f:
                    line = line.strip()
                    if not line.endswith("range"):
                        line = strip_regex.sub(",", line)
                        df.append(f"{file[:4]}," + line)
            with open(f"./cleaned/{file}", "w") as f:
                f.write("\n".join(df))
        os.chdir(os.path.join(ROOTPATH, "cleaned"))
        filenames = glob.iglob("*.pka")
        return filenames

    def read_df(self, filenames, savefile=True):
        """Read all files into one long dataframe.

        Parameters
        ----------
            filenames: iterable
                All the cleaned files from clean_raw_files.
            savefile: boolean, optional

        Returns
        -------
            A pandas.DataFrame with columns:
            Protein, Residue, Position, pK, slope, 1000Chi2
        """
        df = pd.concat((pd.read_csv(f) for f in filenames),
                       ignore_index=True)
        df = df.loc[df.Protein.str.len() == 4]  # Remove error rows
        if savefile:
            df.to_csv(os.path.join(ROOTPATH, "result", "pka_dataframe.csv"), index=False)
        return df


if __name__ == "__main__":
    # Set working directory
    if 1 != len(sys.argv):
        ROOTPATH = os.path.dirname(os.path.realpath(sys.argv[1]))
    else:
        ROOTPATH = os.path.dirname(os.path.realpath(sys.argv[0]))
    os.chdir(ROOTPATH)
    try:
        os.makedirs(os.path.join(ROOTPATH, "cleaned"))
        os.makedirs(os.path.join(ROOTPATH, "result"))
    except OSError as e:
        pass

    x = pka()
    filenames = x.clean_raw_files()
    df = x.read_df(filenames=filenames)
