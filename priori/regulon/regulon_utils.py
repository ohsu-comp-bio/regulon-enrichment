import warnings
warnings.simplefilter("ignore", UserWarning)
import pandas as pd
import dill as pickle
import functools
import os
from sklearn.feature_selection import f_regression, mutual_info_regression
from scipy.stats import spearmanr, pearsonr
import scipy.stats as st
import numpy as np
from tqdm import tqdm
import timeit


def write_pickle(obj, relnm):
    """ Serialize object to pickle and write to disk at relnm

    Args:
        obj (`:obj:`) : Python object to be pickled
        relnm (str) : Relative name/path to pickle on disk

    Returns:
         'Serialized object to disk at {}'.format(relnm)

    """
    with open(relnm, 'wb') as f:
        pickle.dump(obj, f, protocol = -1)
    return 'Serialized object to disk at {}'.format(relnm)


def read_pickle(relnm):
    """ Read serialized object from pickle on disk at relnm

    Args:
        relnm (str) : Relative name/path to pickled object

    Returns:
        obj (`:obj: unpickled object`)

    """

    with open(relnm, 'rb') as f:
        obj = pickle.load(f)
    print('Loaded network from {}'.format(relnm))
    return obj

def read_custom_network(relnm):
    """ Read custom network at relnm

    Args:
        relnm (str) : Relative name/path to tsv object

    Returns:
        obj (`:obj: unpickled object`)

    """

    try:

        df = pd.read_csv(relnm, sep='\t')

        # Check if "regulator" columns exist
        if 'regulator' not in df.columns:
            raise ValueError("The 'regulator' column was not in the input network. Please re-format your network and try again.")

        # Check if "target" columns exist
        if 'target' not in df.columns:
            raise ValueError("The 'target' column was not in the input network. Please re-format your network and try again.")

        # Check for duplicate rows
        if df.duplicated().any():
            dup_rows = df[df.duplicated()]
            # raise ValueError("Duplicate rows were found in your input network. Please re-format your input network.")
            raise ValueError(f"Duplicate rows were found in your input network:\n{dup_rows}. \nPlease re-format your input network.")

    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {relnm}")

    obj = pd.DataFrame({
        'UpGene': df['regulator'],
        'Type': 'controls-expression-of',
        'DownGene': df['target']
    })

    print('Loaded network from {}'.format(relnm))
    return obj


def ensure_dir(relnm):
    """ Accept relative filepath string, create it if it doesnt already exist
        return filepath string

    Args:
        relnm (str) : Relative name/path

    Returns:
        relnm (str)

    """

    d = os.path.join(os.getcwd(), relnm)
    if not os.path.exists(d):
        print('--- path does not exist : {} ---'.format(d))
        print('--- constructing path : {} ---'.format(d))
        os.makedirs(d)

    return relnm