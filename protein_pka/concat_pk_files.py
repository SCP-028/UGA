"""Concatenate protein pK files."""
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

    def read_df(self, filenames, savefile=True):
        """Read all files into one long data frame.

        Parameters
        ----------
            filenames: iterable
                All the cleaned files from clean_raw_files.
            savefile: Boolean, optional

        Returns
        -------
            A pandas.DataFrame with columns:
            Protein, Residue, Position, pK, slope, 1000Chi2
        """
        df = pd.concat((pd.read_csv(f, header=None) for f in filenames),
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
    filenames = glob.iglob(f"{ROOTPATH}/cleaned/*.pka")
    df = x.read_df(filenames=filenames)
