import datetime
import pandas as pd
import scipy.stats as st
import numpy as np


def assign_weights(lh):
    """ Generate normalized likelihood weights and assigns those weights to the absolute gene expression signature

    Args:
        lh (pandas DataFrame): sparse DataFrame indicating likelihood for transcription factors

    Returns:
        nes_wt (pandas series): weight associated for each regulator based on their absolute gene expression signature
        wts (pandas DataFrame): sparse DataFrame indicating likelihood for regulators
    """
    # Generate normalized likelihood weights
    wts = lh.T.divide(lh.max(axis = 1)).T
    wts = wts.fillna(0.0)

    # Absolute gene expression signature
    nes_wt = pd.DataFrame((wts ** 2).sum(axis=1)**.5)

    wts = (wts.T / (wts.T.sum())).T

    return nes_wt, wts


def ensure_overlap(lh, mor, expr):
    """ ensures label overlap with weights and expression matrix

    Args:
        lh (pandas DataFrame): sparse DataFrame indicating likelihood for regulators
        mor (pandas DataFrame): sparse DataFrame indicating mode or regulation for transcription factors
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_feats, n_samps]

    Returns:
        expression, mor, lh (pandas DataFrame): expression, mode of regulation, and likelihood frames, respectfully,
            re-indexed to be concordant with weights associated with regulon
        """

    expression = expr.reindex(mor.columns)
    mor = mor.reindex(expression.index, axis=1)
    lh = lh.reindex(expression.index, axis=1)

    return expression, mor, lh


def deconstruct_regulon(regulon):
    """Returns a pandas DataFrame for likelihood and mode of regulation inferred by compile_total_regulon
    Args:
        regulon (pandas DataFrame): Aracne-AP Regulon file in four column format i.e. Regulator,Target,MoA,likelihood

    Returns:
        lh, mor (pandas DataFrame): sparse DataFrame indicating likelihood and mode or
            regulation for transcription factors
    """

    lh = regulon.pivot(index = 'Regulator', columns = 'Target', values = 'likelihood').fillna(0.0)
    mor = regulon.pivot(index = 'Regulator', columns = 'Target', values = 'MoA').fillna(0.0)

    return lh, mor


def quantile_nes_score(regulon, expr):
    """ Generates quantile transformed rank position enrichment scores

    Args:
        regulon (:obj: `pandas DataFrame`): Aracne-AP Regulon file in four column format i.e.
           Regulator,Target,MoA,likelihood
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_feats, n_samps]

    Returns:
        nes (:obj: `pandas DataFrame`): normalized enrichment scores per regulator
    """

    reset_regulon = regulon.reset_index()
    lh, mor = deconstruct_regulon(reset_regulon)
    expression, mor, lh = ensure_overlap(lh, mor, expr)

    nes_wt, wts = assign_weights(lh)

    pos = expression.index.get_indexer_for(lh.columns)
    t2 = expression.rank() / (expression.shape[0] + 1)
    t2 = pd.DataFrame(st.norm.ppf(t2.iloc[pos, ], loc=0, scale=1), columns=t2.columns, index=t2.index)
    sum1 = (mor * wts).dot(t2)
    nes = pd.DataFrame(sum1.values * nes_wt.values, columns=sum1.columns, index=sum1.index)

    return nes


def load_quantile(regulon, expr, cohort):
    """ return pandas series of quantile enrichment scores for a given regulator

    Args:
        regulon (:obj: `pandas DataFrame`): Aracne-AP Regulon file in four column format i.e.
           Regulator,Target,MoA,likelihood
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_feats, n_samps]
        cohort (str) : name of cohort to associate with compiled regulon

    Returns:
        nes (obj: pandas series): series of quantile enrichment scores for a give regulator
    """
    quantile_nes = os.path.join(dirname, '../experiments/{0}/data/{0}_quantile_ranks.pkl'.format(cohort))
    if os.path.isfile(quantile_nes):
        print('--- Loading quantile normalization scores ---')
        nes = read_pickle(quantile_nes)
    else:
        print('--- Generating quantile normalization scores ---')
        nes = quantile_nes_score(regulon, expr.T)

    return nes


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
    sub_regul = regulon.loc[regulator]
    sub_expr = expr.reindex(sub_regul.Target.values, axis = 1)

    # Subset expression values by up and down regulated targets
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

    total_ranked = expr.rank(method = 'max', ascending = ascending, axis = 1)
    moa_frame = regulon.loc[regulator, ].loc[regulon.loc[regulator, ].Target.isin(expr_sub.columns),
                                             ['Target', 'MoA', 'likelihood']].reset_index()

    moa_frame.index = moa_frame.Target
    moa_frame.likelihood = moa_frame.likelihood / moa_frame.likelihood.max()

    moa_frame['weights'] = moa_frame.MoA * moa_frame.likelihood

    moa_frame = moa_frame.loc[:, 'MoA'].to_frame().T

    ranks = total_ranked.loc[:, expr_sub.columns]

    weighted_ranks = pd.np.multiply(ranks, moa_frame).sum(axis = 1).to_frame()

    # Store minimum rank for samples - this value is equivalent to the total number of targets in expr_sub i.e. if all
    # genes for a particular sample rank first the rank_min = 1.0 * #genes
    # rank_min = 1.0 * expr_sub.shape[1]

    rank_min = weighted_ranks.min().values[0]

    # Store maximum rank for samples - this value is equivalent to the total number of targets in expr_sub i.e. if all
    # genes for a particular sample rank last the rank_min = #samples * #genes

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

    down_normed = pd.DataFrame(down_reg_ordered.loc[:, 0].values, columns = ['down-regulated-targets'],
                               index = down_reg_ordered.index)
    down_normed = down_normed.fillna(0.0)

    up_normed = pd.DataFrame(up_reg_ordered.loc[:, 0].values, columns = ['up-regulated-targets'],
                             index = up_reg_ordered.index)
    up_normed = up_normed.fillna(0.0)

    join_r = pd.concat([down_normed, up_normed], axis = 1)

    join_r.columns = ['down-regulated-targets', 'up-regulated-targets']
    zframe = ((join_r['down-regulated-targets'] * -1) + join_r['up-regulated-targets']).to_frame()

    zframe.columns = [regulator]
    zframe[regulator] = st.zscore(zframe[regulator])
    zframe = (zframe - zframe.median()) / zframe.std()

    return zframe


def format_delta(down_reg_sub, up_reg_sub):
    """ Take the mean difference from up/down - regulated targets and scale

    Args:
        down_reg_sub (:obj: `pandas DataFrame`) : pandas DataFrame of down regulated targets regulator normed expression
        up_reg_sub (:obj: `pandas DataFrame`) : pandas DataFrame of up regulated targets regulator normed expression

    Returns:
        delta (:obj: `pandas DataFrame`) : pandas DataFrame of the mean delta difference of regulator targets
        per sample

    """

    up_reg_sub = pd.DataFrame(st.zscore(up_reg_sub, axis = 1), columns = up_reg_sub.columns, index = up_reg_sub.index)
    up_reg_normed = up_reg_sub - up_reg_sub.median()
    up_reg_normed = (up_reg_normed.median(axis = 1) / up_reg_normed.median(axis = 1).max()).\
        fillna(0.0).replace([np.inf, -np.inf], 0.0)

    down_reg_sub = pd.DataFrame(st.zscore(down_reg_sub), columns = down_reg_sub.columns, index = down_reg_sub.index)
    down_reg_normed = down_reg_sub + down_reg_sub.median()
    down_reg_normed = (down_reg_normed.median(axis = 1) / down_reg_normed.median(axis = 1).max()).\
        fillna(0.0).replace([np.inf, -np.inf], 0.0)

    delta = up_reg_normed - down_reg_normed
    delta = delta.to_frame()
    delta.columns = ['Delta']
    delta['Delta'] = st.zscore(delta)

    return delta


def score_enrichment(regulator, expr, regulon, quant_nes):
    """ Function to subset and generate regulator activity scores based
        on rank ordering of up-regulated and down-regulated targets

    Args:
        regulator (str) : Regulator to subset expression frame by
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_samps, n_feats]
        regulon (:obj: `pandas DataFrame`): pandas DataFrame of regulon returned by compile_regulon
            with columns ['Target', 'MoA', 'likelihood']
        quant_nes (obj: `pandas DataFrame`): quantile enrichment scores for regulators
    Return:
        enrichment_score (:obj: `pandas DataFrame`): pandas DataFrame of activity scores for specified regulator

    """
    print(regulator)
    down_reg_sub, up_reg_sub = subset_regulon(regulator, regulon, expr)

    # Rank up and down regulated targets by z-scores. Sum rank values across rows
    # (Compute numerical data ranks [1 through n] along axis)
    # and sort samples lowest to highest summed rank score.

    down_reg_ordered = rank_and_order_total(down_reg_sub, regulator, regulon, ascending=False, expr=expr)
    up_reg_ordered = rank_and_order_total(up_reg_sub, regulator, regulon, ascending=True, expr=expr)

    zframe = format_nes_frame(down_reg_ordered, up_reg_ordered, regulator)
    delta = format_delta(down_reg_sub, up_reg_sub)

    zframe[regulator] = zframe.values + delta.values

    enrichment_score = zframe[regulator] + quant_nes.loc[regulator]

    return enrichment_score


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


def generate_enrichment_scores(expr_f, cohort, norm_type = 'robust', feature = True, sample = False,
                               thresh_filter = 0.4, scale = True, regulon_size = 15):
    """ Runs expression and regulon_utils functions to generate cohort specific regulon and enrichment scores

    Args:
        expr_f (str): absolute path to tab delimited expression file of shape = [n_features, n_samples]
        cohort (str) : name of cohort to associate with compiled regulon and enrichment scores
        norm_type (str): Scaler to normalized features/samples by: standard | robust | minmax | quant
        feature (bool): Scale expression data by features
        sample (bool): Scale expression data by both features and samples
        thresh_filter (float): Prior to normalization remove features that do not have the mean unit of
            a feature (i.e. 1 tpm) is greater than {thresh_filter}
        scale (bool): optional arg to avoid scaling dataset if data set has been normalized prior to analysis
        regulon_size (int) : required number of downstream interactions for a give regulator

    Returns:
        None

    """
    input_args = locals()
    total_enrichment_nes = os.path.join(dirname, '../experiments/{0}/data/{0}_total_enrichment.pkl'.format(cohort))
    if os.path.isfile(total_enrichment_nes):
        print('--- Regulon enrichment scores pre-computed ---')
        print(total_enrichment_nes)

    else:
        non_scaled_expr = load_expr(expr_f)

        expr = load_scaled_expr(non_scaled_expr, cohort = cohort, norm_type = norm_type, feature = feature,
                                sample = sample, thresh_filter = thresh_filter, scale = scale)

        regulon = generate_bolstered_regulon(expr, cohort, regulon_size = regulon_size)
        quant_nes = load_quantile(regulon, expr, cohort)
        regulators = regulon.index.unique()

        print('--- Calculating regulon enrichment scores ---')
        nes_list = list(map(functools.partial(score_enrichment, expr=expr, regulon = regulon, quant_nes=quant_nes),
                            tqdm(regulators)))
        total_enrichment = pd.concat(nes_list, axis=1)

        relnm = os.path.join(dirname, '../experiments/{0}/data'.format(cohort))

        ensure_dir(relnm)
        write_pickle(total_enrichment, os.path.join(relnm, '{}_total_enrichment.pkl'.format(cohort)))
        logger(**input_args)
