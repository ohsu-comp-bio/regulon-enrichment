import datetime
import pandas as pd
import scipy.stats as st
import numpy as np

def subset_regulon(regulator, regulon, expr):
    """ Subset expression frame by regulator targets expressed in expression frame and by mode of regulation

    Args:
        regulator (str) : Regulator to subset expression frame by
        regulon (:obj: `pandas DataFrame`) : pandas DataFrame of Regulator-Target interactions
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_samps, n_feats]

    Returns:
        down_reg_sub (:obj: `pandas DataFrame`) : pandas DataFrame of down regulated targets regulator normed expression
        values
        up_reg_sub (:obj: `pandas DataFrame`) : pandas DataFrame of up regulated targets regulator normed expression
        values

    """
    # Subset regulon by transcription factor
    sub_regul = regulon.loc[regulator]
    sub_expr = expr.reindex(sub_regul.Target.values, axis = 1)

    # Identify up-regulated and down-regulated targets
    down_reg_sub = sub_expr.loc[:, (sub_regul.MoA < 0.0).values]
    up_reg_sub = sub_expr.loc[:, (sub_regul.MoA > 0.0).values]

    return down_reg_sub, up_reg_sub


def rank_and_order_total(expr_sub, regulator, regulon, ascending, expr):
    """ Rank and order transcription factor targets expression frame

    Args:
        expr_sub (:obj: `pandas DataFrame`) : pandas DataFrame of regulated targets regulator normed expression
        regulator (str) : Regulator to subset expression frame by
        regulon (:obj: `pandas DataFrame`): pandas DataFrame of regulon returned by compile_regulon
            with columns ['Target', 'MoA', 'likelihood']
        ascending (bool): Boolean flag to rank regulon gene set via ascending/descending manner
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_samps, n_feats]

    Returns:
        rank_ordered (:obj: `pandas DataFrame`) : pandas DataFrame of regulated targets regulator normed expression

    """
    # Rank target genes (either up-regulated or down-regulated)
    total_ranked = expr.rank(method = 'max', ascending = ascending, axis = 1)
    moa_frame = regulon.loc[regulator, ].loc[regulon.loc[regulator, ].Target.isin(expr_sub.columns),
                                             ['Target', 'MoA', 'likelihood']].reset_index()

    # Normalize weights by the likelihood (adjusted p-value)
    moa_frame.index = moa_frame.Target
    moa_frame.likelihood = moa_frame.likelihood / moa_frame.likelihood.max()
    moa_frame['weights'] = moa_frame.MoA * moa_frame.likelihood
    moa_frame = moa_frame.loc[:, 'MoA'].to_frame().T
    ranks = total_ranked.loc[:, expr_sub.columns]

    # Align
    aligned = ranks.align(moa_frame, axis = 1)

    # Multiply weights and ranks
    weighted_ranks = np.multiply(aligned[0], np.asarray(aligned[1])).sum(axis = 1).to_frame()

    # Store minimum and maximum ranks for samples
    rank_min = weighted_ranks.min().values[0]
    rank_max = weighted_ranks.max().values[0]

    weighted_ranks['min'] = rank_min
    weighted_ranks['max'] = rank_max

    return weighted_ranks


def format_nes_frame(down_reg_ordered, up_reg_ordered, regulator):
    """ Function to concatenate and sum down and up regulated z-score rankings

    Args:
        down_reg_ordered (:obj: `pandas DataFrame`) : pandas DataFrame of z-scores for down regulated
        targets of regulator
        up_reg_ordered (:obj: `pandas DataFrame`) : pandas DataFrame of z-scores for up regulated
        targets of regulator
        regulator (str) : Regulator that controls the activity of a regulon

    Returns:
        zframe (:obj: `pandas DataFrame`) : pandas DataFrame of average z-scores for up and down-regulated targets

    """

    # Identify up-regulated and down-regulated targets
    down_normed = pd.DataFrame(down_reg_ordered.loc[:, 0].values, columns = ['down-regulated-targets'],
                               index = down_reg_ordered.index)
    down_normed = down_normed.fillna(0.0)

    up_normed = pd.DataFrame(up_reg_ordered.loc[:, 0].values, columns = ['up-regulated-targets'],
                             index = up_reg_ordered.index)
    up_normed = up_normed.fillna(0.0)

    # Combine into single matrix
    join_r = pd.concat([down_normed, up_normed], axis = 1)

    # Combine values, weighting the downregulated targets by -1
    join_r.columns = ['down-regulated-targets', 'up-regulated-targets']
    zframe = ((join_r['down-regulated-targets'] * -1) + join_r['up-regulated-targets']).to_frame()

    # z-transform by each transcription factor
    zframe.columns = [regulator]
    zframe[regulator] = st.zscore(zframe[regulator])
    zframe = (zframe - zframe.median()) / zframe.std()

    return zframe

def score_enrichment(regulator, expr, regulon):
    """ Function to subset and generate regulator activity scores based
        on rank ordering of up-regulated and down-regulated targets

    Args:
        regulator (str) : Regulator to subset expression frame by
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_samps, n_feats]
        regulon (:obj: `pandas DataFrame`): pandas DataFrame of regulon returned by compile_regulon
            with columns ['Target', 'MoA', 'likelihood']
    Return:
        enrichment_score (:obj: `pandas DataFrame`): pandas DataFrame of activity scores for specified regulator

    """
    print(regulator)

    # Subset network and expression matrix by the regulon
    down_reg_sub, up_reg_sub = subset_regulon(regulator, regulon, expr)

    # Rank up and down regulated targets by z-scores. Sum rank values across rows
    # (Compute numerical data ranks [1 through n] along axis) and sort samples lowest to highest summed rank score.
    down_reg_ordered = rank_and_order_total(down_reg_sub, regulator, regulon, ascending=False, expr=expr)
    up_reg_ordered = rank_and_order_total(up_reg_sub, regulator, regulon, ascending=True, expr=expr)

    # Concatenate and sum up-regulated and down-regulated z-score rankings
    zframe = format_nes_frame(down_reg_ordered, up_reg_ordered, regulator)

    return zframe


def logger(**kwargs):
    """ Generates a log file of arguments passed to EnrichR.py

    Args:
        **kwargs: paired key word arguments

    Returns:
        None
    """
    cohort = kwargs['cohort']
    relnm = os.path.join(dirname, '../experiments/{0}/data'.format(cohort))
    now = datetime.datetime.now()
    ensure_dir(relnm)
    out_f = open(os.path.join(relnm, '{}_kwargs.txt'.format(cohort)), 'w')
    out_f.write("EnrichR generated regulon, enrichment scores and scaled expression data-set compiled on "
                "{} with the following **kwargs \n".
                format(now.strftime("%Y-%m-%d %H:%M")))
    for k, v in kwargs.items():
        out_f.write('* {} : {} \n'.format(k, v))
    out_f.close()
