import warnings
warnings.simplefilter("ignore", UserWarning)
import pandas as pd
import dill as pickle
import functools
import os
from sklearn.feature_selection import f_regression, mutual_info_regression
from sklearn.mixture import BayesianGaussianMixture as GMM
from scipy.stats import spearmanr, pearsonr
import scipy.stats as st
import numpy as np
from tqdm import tqdm
import timeit



def load_sif():
    return pd.read_csv(sif_file, names = ['UpGene', 'Type', 'DownGene'], sep = '\t', header = None)


def filter_sif(sif, intx_type = 'controls-expression-of'):
    return sif[(sif['Type'] == intx_type)]


def load_secondary_itx_sif():
    """ Load precompiled secondary interaction sif

    Returns:
        (pandas.DataFrame): pandas.DataFrame obj of length: n interactions and
            columns: ['UpGene','Type',DownGene']

    """
    return pd.read_csv(sec_intx_file, names = ['UpGene', 'Type', 'DownGene'], sep = '\t', header = None)


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
    print('Loaded object from disk at {}'.format(relnm))
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


def traverse_interactions(regulator, filt_sif):
    """ Parse interaction network and add secondary interactions on a per regulator basis

    Args:
        regulator (str): Regulator to expand interaction network
        filt_sif (pandas.DataFrame): pandas.DataFrame obj of length: n interactions and
            columns: ['UpGene','Type',DownGene']

    Returns:
        comb_idx (pandas.DataFrame):pandas.DataFrame obj of length: n interactions + secondary interactions and
            columns: ['UpGene','Type',DownGene']

    """
    sub_reg = filt_sif[(filt_sif.UpGene == regulator)]
    down_genes = sub_reg.DownGene.unique()
    secondary_itx = filt_sif[(filt_sif.UpGene.isin(down_genes))]
    secondary_itx.UpGene = regulator
    comb_idx = pd.concat([sub_reg, secondary_itx])
    comb_idx.Type = 'controls-expression-of'
    comb_idx = comb_idx.drop_duplicates()
    comb_idx = comb_idx[(comb_idx.DownGene != regulator)]

    return comb_idx


def generate_expanded_regulon():
    """ Generates an expanded Pathway Commons regulon with secondary down-stream interactions for
        regulators that control the expression of other regulators

    Returns:
        Nothing - Generates a pickled pandas dataframe for future reference/use

    """
    print('--- Generating regulon with primary and secondary interactions ---')
    sif = load_sif()
    filt_sif = filter_sif(sif)
    regulators = filt_sif.UpGene.unique()

    regulon_list = list(map(functools.partial(traverse_interactions, filt_sif = filt_sif), regulators))

    regulon = pd.concat(regulon_list)
    regulon.set_index('UpGene', inplace = True)
    regulon.reset_index(inplace=True)
    print('---- Regulon constructed ---')

    write_pickle(regulon, '../data/secondary_intx_regulon.pkl')


def bgm_moa(regul_weights):
    """ Fits regulon mode of activation weights to a bayesian gaussian mixture model with three components and computes
        the probability of the three distributions (repression, non-influential, activation) for each regulator

    Args:
        regul_weights (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions between regulator and
            downstream members of its regulon of shape [len(Target), ['Regulator','Target','MoA','likelihood']

    Returns:

    """

    g = GMM(n_components = 3, max_iter = 1000)

    sub_reg = regul_weights.copy()
    sub_reg_vals = sub_reg.MoA.values.reshape(-1, 1)
    g.fit(sub_reg_vals)

    mu = g.means_.flatten()
    sigma = np.sqrt(g.covariances_).flatten()
    fit = sorted(list(zip(mu, sigma)))

    activation = (st.norm.cdf(sub_reg_vals, fit[2][0], fit[2][1]))
    repression = 1 - st.norm.cdf(sub_reg_vals, fit[0][0], fit[0][1])
    total_lower = 1 - st.norm.cdf(sub_reg_vals, fit[1][0], fit[1][1])
    total_upper = (st.norm.cdf(sub_reg_vals, fit[1][0], fit[1][1]))

    copy_target = sub_reg.copy()
    copy_target['up'] = 0
    copy_target['down'] = 0
    copy_target.loc[(copy_target.MoA >= 0), 'up'] = 1
    copy_target.loc[(copy_target.MoA <= 0), 'down'] = 1

    up_moa = copy_target.up.values.reshape(copy_target.shape[0], 1)
    down_moa = copy_target.down.values.reshape(copy_target.shape[0], 1)

    Mode = (activation / (repression + total_lower + activation) * up_moa) -\
           (repression / (repression + total_upper + activation) * down_moa)

    return Mode


def prune_regulon(expr, regulon, regulon_size):
    """ Prunes regulon with secondary interactions that do not meet the necessary number of downstream interactions
    metric {regulon_size}

    Args:
        expr (pandas DataFrame obj): pandas DataFrame of [n_samples, n_features]
        regulon (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions between regulator and
            downstream members of its regulon of shape [len(Target), ['Regulator','Target','MoA','likelihood']
        regulon_size (int) : number of downstream interactions required for a given regulator in order to calculate
            enrichment score

    Returns:
        filtered_regulon (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions between regulator
            and downstream members of its regulon of shape [len(Target), ['Regulator','Target','MoA','likelihood']

    """

    expr_filtered_regulon = regulon[((regulon.UpGene.isin(expr.columns)) & (regulon.DownGene.isin(expr.columns)))]
    expr_filtered_regulon.set_index('UpGene', inplace=True)
    idx = (expr_filtered_regulon.index.value_counts() >= regulon_size)
    filt_idx = idx[idx==True]
    filtered_regulon = expr_filtered_regulon.loc[filt_idx.index]
    filtered_regulon.reset_index(inplace=True)

    return filtered_regulon


def regulon_weight_assignment(regulator, expr, filtered_regulon):
    """ Assigns probability and weights for regulator - target interactions

    Args:
        regulator (str): Regulator to expand interaction network
        expr (:obj: `pandas DataFrame`) : pandas DataFrame containing scaled expression data of
            shape [n_samples, n_features]
        filtered_regulon (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions between regulator
            and downstream members of its regulon of shape [len(Target), ['Regulator','Target','MoA','likelihood']

    Returns:
        regul_weights (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions between regulator and
            downstream members of its regulon of shape [len(Target), ['Regulator','Target','MoA','likelihood']

    """

    sub_reg = filtered_regulon[(filtered_regulon['UpGene'] == regulator)]

    X = expr.reindex(sub_reg.DownGene.values, axis = 1).dropna(axis = 1)
    y = expr.reindex([regulator], axis = 1)

    spr_results = X.apply(lambda col: spearmanr(col, y.iloc[:, 0]), axis = 0).apply(pd.Series)

    spr_result = spr_results[0]
    spr_pvalues = spr_results[1]

    f_test, _ = f_regression(X, y.values.ravel())
    weights = f_test

    weights_spr = weights + abs(spr_result)

    regul_weights = (weights_spr * np.sign(spr_result)).to_frame()
    regul_weights.columns = ['MoA']
    regul_weights.index.name = 'Target'
    regul_weights.reset_index(inplace = True)
    regul_weights['Regulator'] = regulator
    regul_weights['likelihood'] = spr_pvalues.values
    regul_weights = regul_weights.reindex(['Regulator', 'Target', 'MoA', 'likelihood'], axis = 1)
    regul_weights.set_index('Regulator', inplace = True)
    regul_weights = regul_weights[~np.isinf(regul_weights.MoA)]

    return regul_weights


def structure_weights(regulator, pruned_regulon, f_statistics, r_frame, p_frame):
    """ Calculates weights associated with regulators. Weights are the summation of the F-statistic and absolute
        spearman correlation coefficient. The weight retains the sign of the spearman correlation coefficient.

    Args:
        regulator (str): A feature to assign weights to downstream interactions
        pruned_regulon (:obj:`pd.DataFrame`, shape = [n_interactions, 3]
        f_statistics (dict) : Dictionary with key:{regulator} key and
        r_frame (:obj:`pd.DataFrame`), shape = [n_features, n_features]
        p_frame (:obj:`pd.DataFrame`), shape = [n_features, n_features]

    Returns:
        weights_ordered (:obj:`pd.DataFrame`), shape = [n_interactions, 3]

    """

    sub_regul = pruned_regulon[(pruned_regulon['UpGene'] == regulator)]
    targs = sub_regul.DownGene
    p_ = p_frame.loc[targs, regulator]
    p_.name = 'likelihood'
    f_ = f_statistics[regulator][0]
    r_ = r_frame.loc[targs, regulator]
    w_ = (f_ + abs(r_)) * np.sign(r_)
    w_.index.name = 'Target'
    w_.name = 'MoA'
    weights = w_.to_frame()
    weights['likelihood'] = p_
    weights['Regulator'] = regulator
    weights_ordered = weights.reset_index().reindex(['Regulator', 'Target', 'MoA', 'likelihood'],
                                                    axis = 1).set_index('Regulator')

    return weights_ordered


def generate_bolstered_regulon(expr, cohort, regulon_size=15):
    """ Calculate weights for PC regulon and a dataset using mutual information, f-statistic to test for linear
        relationships, and the spearman correlation coefficient to determine the mode of regulation

    Args:
        expr (:obj: `pandas DataFrame`) : pandas DataFrame containing scaled expression data of
            shape [n_samples, n_features]
        cohort (str) : name of cohort to associate with compiled regulon
        regulon_size (int) : required number of downstream interactions for a give regulator

    Returns:
        regul_weights (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions between regulator and
            downstream members of its regulon of shape [len(Target), ['Regulator','Target','MoA','likelihood']

    """
    bolstered_relnm = os.path.join(dirname, '../experiments/{0}/data/{0}_bolstered_regulon.pkl'.format(cohort))

    # Check to see if bolstered regulon exists
    if os.path.isfile(bolstered_relnm):
        print('--- loading context specific regulon ---')
        total_regulon = read_pickle(bolstered_relnm)

    else:
        if os.path.isfile(sec_intx_file):
            print('--- loading unfiltered regulon ---')
            regulon = read_pickle(sec_intx_file)
        else:
            generate_expanded_regulon()
            regulon = read_pickle(sec_intx_file)

        print('--- pruning regulon ---')
        filtered_regulon = prune_regulon(expr, regulon, regulon_size)
        regulators = filtered_regulon.UpGene.unique()

        print('--- compiling regulon of {} regulators and {} interactions with a minimum of {} interactions ---'.
              format(len(regulators), filtered_regulon.shape[0], regulon_size))

        regulon_list = list(map(functools.partial(regulon_weight_assignment, expr=expr,
                                                  filtered_regulon = filtered_regulon), tqdm(regulators)))
        total_regulon = pd.concat(regulon_list)

        relnm = os.path.join(dirname, '../experiments/{0}/data'.format(cohort))

        ensure_dir(relnm)
        write_pickle(total_regulon, os.path.join(relnm, '{}_bolstered_regulon.pkl'.format(cohort)))

    return total_regulon
