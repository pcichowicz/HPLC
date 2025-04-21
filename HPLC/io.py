# File I/O: loading chromatograms, etc.
import numpy as np
import pandas as pd

def load_chromatogram(fname, cols, delimiter=',', dropna=False):
    """
    Load and parse file containing chromatogram and returns as Pandas Dataframe.

    Parameters
    ----------
    :param fname: `str`
        The path to chromatogram file, must be text file (i.e. not `.xlsx`).
    :param cols: `list` or `dict`
        The columns present in the text file. Dict will allow for renaming.
        with `key` -> `value`.
    :param delimiter: `str`
        Delimiter character separating columns (i.e. `,`, `\t`)
    :param dropna: `bool`
        If True, drops NaN's from chromatogram.

    Returns
    -------
    df : `pandas.core.frame.DataFrame`
        The chromatograph loaded as a Pandas Dataframe
    """

    if type(cols) == dict:
        _colnames = list(cols.keys())
    else:
        _colnames = cols
    skip = 0
    num = 0

    if len(_colnames) != 0:
        with open(fname, 'r') as f:
            _lines = f.readlines()
            halted = False
            for line in _lines:
                if np.array([nom.lower() in line.lower() for nom in _colnames]).all():
                    halted = True
                    num +=1
                else:
                    if num == 0:
                        skip +=1
            if not halted:
                raise ValueError(
                    "Column name(s) not found in file provided"
                )
    if num >1:
        raise RuntimeError(
            "More than one chromatogram is not supported. Provide file with only one chromatogram."
        )

    # Finally load in dataframe with proper parameters given
    df = pd.read_csv(fname, skiprows=skip, delimiter=delimiter)
    if type(cols) == dict:
        df.rename(columns = cols, inplace = True)
        _colnames = list(cols.values())
    if dropna:
        df.dropna(inplace=True)

    df = df[_colnames]
    return df