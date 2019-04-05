import pandas as pd
import os
import sys
bcell_project_dir = os.path.abspath('../..')
sys.path.append(bcell_project_dir)

from enrichment.regulon_enrichment import *
pri_intx_file = DATA_PATH + '/primary_intx_regulon.pkl'


def load_entrez_map():
    """ Load entrez to gene symbol map

    Returns:
        map_frame (`obj:pandas DataFrame`): pandas dataframe used to translate from entrez to gene symbol

    """

    data_path = '/Users/estabroj/PycharmProjects/regulon_enrichment/experiments/Bcell/data'
    map_frame = pd.read_csv(os.path.join(data_path,'entrez_to_gene.txt'),sep='\t',index_col=0)
    map_frame.index = map_frame.index.map(str)

    return map_frame


def load_array_expression(dataset):
    data_path = '/Users/estabroj/PycharmProjects/regulon_enrichment/experiments/Bcell/data'

    entrez_map = load_entrez_map()

    knockdown_frame = pd.read_csv(os.path.join(data_path,dataset),sep='\t',index_col=0)
    knockdown_frame.index.name='Entrez'
    knockdown_frame.index = knockdown_frame.index.map(str)
    if '.' in knockdown_frame.index[0]:
        idx_split = '.'
    else:
        idx_split = '_'

    knockdown_frame.reset_index(inplace=True)
    knockdown_frame.Entrez = knockdown_frame.Entrez.str.split(idx_split,expand=True)[0]
    knockdown_frame.set_index('Entrez',inplace=True)
    knockdown_joined = knockdown_frame.join(entrez_map)
    knockdown_joined.set_index('entrez2gene', inplace = True)
    knockdown_joined = knockdown_joined.groupby('entrez2gene').median()
    knockdown_joined.index.name = 'Gene'

    return knockdown_joined


def regulon_weight_assignment_remove_mi_no_gm_no_scale(regulator, expr, filtered_regulon):
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
    # regul_weights = traverse_interactions_(regul_weights = regul_weights, filtered_regulon = filtered_regulon, regulator = regulator, expr = expr, intx_thresh=300)
    regul_weights = regul_weights[~np.isinf(regul_weights.MoA)]

    return regul_weights


def rwa_no_traverse(regulator, expr, filtered_regulon):
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

    sub_reg = filtered_regulon

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
    print(regul_weights.shape)

    return regul_weights


def traverse_interactions_(regul_weights, filtered_regulon, regulator, expr, intx_thresh=300):
    """ Parse interaction network and add secondary interactions on a per regulator basis

    Args:
        regulator (str): Regulator to expand interaction network
        filt_sif (pandas.DataFrame): pandas.DataFrame obj of length: n interactions and
            columns: ['UpGene','Type',DownGene']

    Returns:
        comb_idx (pandas.DataFrame):pandas.DataFrame obj of length: n interactions + secondary interactions and
            columns: ['UpGene','Type',DownGene']

    """

    targets = regul_weights.loc[:,['Target','MoA']].set_index('Target').MoA.abs().sort_values(ascending=False)
    filt_targets = targets.index[targets.index.isin(filtered_regulon.UpGene)].tolist()
    concat_frame = pd.DataFrame(columns = filtered_regulon.columns)
    if len(filt_targets) > 1 and regul_weights.shape[0] < intx_thresh:
        n = 0
        for t in filt_targets:
            if n < intx_thresh:
                sub_regul = filtered_regulon[(filtered_regulon['UpGene']==t)].copy()
                sub_regul['UpGene'] = regulator
                filt_sub = sub_regul[~sub_regul.DownGene.isin(regul_weights.Target)]
                concat_frame = pd.concat([concat_frame, filt_sub])
                n = concat_frame.shape[0] + regul_weights.shape[0]

        sub_regul_weights = rwa_no_traverse(regulator, expr, concat_frame)
        regul_weights = pd.concat([regul_weights, sub_regul_weights])
        regul_weights = regul_weights[~regul_weights.Target.duplicated()]

    return regul_weights


def generate_regulon():
    """ Generates an expanded Pathway Commons regulon with secondary down-stream interactions for
        regulators that control the expression of other regulators

    Returns:
        Nothing - Generates a pickled pandas dataframe for future reference/use

    """
    print('--- Generating regulon with primary interactions ---')
    sif = load_sif()
    filt_sif = filter_sif(sif)
    filt_sif.set_index('UpGene', inplace = True)
    filt_sif.reset_index(inplace=True)
    write_pickle(filt_sif, '../data/primary_intx_regulon.pkl')


def generate_bolstered_regulon(expr, cohort, regulon_size=5):
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
        filtered_regulon = prune_regulon(expr = expr, regulon = regulon, regulon_size = regulon_size)
        regulators = filtered_regulon.UpGene.unique()

        print('--- compiling regulon of {} regulators and {} interactions with a minimum of {} interactions ---'.
              format(len(regulators), filtered_regulon.shape[0], regulon_size))

        regulon_list = list(map(functools.partial(regulon_weight_assignment_remove_mi_no_gm_no_scale, expr = expr,
                                                  filtered_regulon = filtered_regulon), tqdm(regulators)))
        total_regulon = pd.concat(regulon_list)

        relnm = os.path.join(dirname, '../experiments/{0}/data'.format(cohort))

        ensure_dir(relnm)
        write_pickle(total_regulon, os.path.join(relnm, '{}_bolstered_regulon.pkl'.format(cohort)))

    return total_regulon


def generate_enrichment_scores_array(non_scaled_expr, cohort, norm_type = 'robust', feature = True, sample = False,
                               thresh_filter = 0.4, scale = True, regulon_size = 5):
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

        expr = load_scaled_expr(non_scaled_expr, cohort = cohort, norm_type = norm_type, feature = feature,
                                sample = sample, thresh_filter = thresh_filter, scale = scale)


        regulon = generate_bolstered_regulon(expr, cohort, regulon_size = regulon_size)

        # ###
        # idx = regulon.MoA.abs() > 1.0
        # regulon = regulon.loc[idx]
        # regulon = regulon[(regulon.index.value_counts() >= regulon_size)]
        # ###

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


def subset_by_treatment_scale(expr, string_match = 'ctrl', axis = 1):
    if axis == 1:
        idx = expr.columns.str.contains(string_match)
        expr = expr.loc[:, idx]
    else:
        idx = expr.index.str.contains(string_match)
        expr = expr.loc[idx, :]
    return expr


def scale_dataset(expr):
    """

    Args:
        expr:

    Returns:

    """
    # print(expr.iloc[:,3])
    # print('--- Scaling variance of controls ---')
    # print('--- Cohort : {} ---'.format(cohort))
    # scaled_expr = scale_dataset(expr=expr.T)
    # expr = scaled_expr.T

    mean_control = subset_by_treatment_scale(expr=expr).mean(axis = 1)
    mean_difference = expr.T - mean_control
    control = subset_by_treatment_scale(expr)
    square_difference = ((control.T - mean_control) ** 2).T
    dot_prod = square_difference.dot(np.repeat(1,square_difference.shape[1]))
    variance = dot_prod.divide(control.shape[1]-1)
    std = np.sqrt(variance)
    signature = (mean_difference.divide(std)).T

    return signature


def generate_bcell_enrichment():
    datasets = {'MEF2B': 'p3hr1_mef2b_bcell.txt', 'FOXM1': 'st486_foxm1_bcell.txt', 'MYB': 'st486_myb_bcell.txt',
                'BCL6_ly7': 'ly7_bcl6_bcell.txt', 'BCL6_pfeiffer': 'pfeiffer_bcl6_bcell.txt',
                'STAT3': 'snb19_stat3_bcell.txt'}
    data_path = '/Users/estabroj/PycharmProjects/regulon_enrichment/experiments/Bcell/data'
    for cohort in datasets:
        print(cohort)
        expr_f = os.path.join(data_path,datasets[cohort])
        print(expr_f)
        non_scaled_expr = load_array_expression(expr_f)
        generate_enrichment_scores_array(non_scaled_expr, cohort, norm_type = 'robust', feature = True, sample = False,
                                       thresh_filter = 0.4, scale = True, regulon_size = 5)


def subset_by_treatment(expr, string_match='ctrl'):
    idx = expr.index.str.contains(string_match)
    control = expr.loc[idx,:]
    trx = expr.loc[~idx,:]

    return control, trx


def generate_rankings_expr():

    enrichment_files = {'MEF2B': 'MEF2B_robust_True_False_0.4_True_frame.pkl',
                'FOXM1': 'FOXM1_robust_True_False_0.4_True_frame.pkl',
                'MYB': 'MYB_robust_True_False_0.4_True_frame.pkl',
                'BCL6_ly7': 'BCL6_ly7_robust_True_False_0.4_True_frame.pkl',
                'BCL6_pfeiffer': 'BCL6_pfeiffer_robust_True_False_0.4_True_frame.pkl',
                'STAT3': 'STAT3_robust_True_False_0.4_True_frame.pkl'}

    bcell_data_path = '/Users/estabroj/PycharmProjects/regulon_enrichment/experiments/Bcell'
    stats_dict = {}
    for file in enrichment_files:
        enrichment_file = os.path.join(bcell_data_path, file, 'data', enrichment_files[file])
        reg_enrichment = read_pickle(enrichment_file)
        control, trx = subset_by_treatment(reg_enrichment)
        ttest_results = st.ttest_ind(control, trx)
        stats_frame = pd.DataFrame(np.array([ttest_results.pvalue, ttest_results.statistic]),
                                   index = ['pvalue','t-statistic'], columns = trx.columns).T
        stats_frame.sort_values(by = 'pvalue', ascending = True, inplace = True)
        if file not in stats_dict:
            stats_dict[file] = stats_frame
    out_p = os.path.join(bcell_data_path,'data','expression_stats_dictionary.pkl')
    write_pickle(stats_dict,out_p)


def generate_rankings(path=None):
    if path == None:
        path = os.getcwd()
    enrichment_files = {'MEF2B': 'MEF2B_total_enrichment.pkl',
                        'FOXM1': 'FOXM1_total_enrichment.pkl',
                        'MYB': 'MYB_total_enrichment.pkl',
                        'BCL6_ly7': 'BCL6_ly7_total_enrichment.pkl',
                        'BCL6_pfeiffer': 'BCL6_pfeiffer_total_enrichment.pkl',
                        'STAT3': 'STAT3_total_enrichment.pkl'}
    stats_dict = {}
    for file in enrichment_files:
        enrichment_file = os.path.join(path, file, 'data', enrichment_files[file])
        reg_enrichment = read_pickle(enrichment_file)
        control, trx = subset_by_treatment(reg_enrichment)
        ttest_results = st.ttest_ind(control, trx)
        stats_frame = pd.DataFrame(np.array([ttest_results.pvalue, ttest_results.statistic]),
                                   index = ['pvalue', 't-statistic'], columns = trx.columns).T
        stats_frame.sort_values(by = 'pvalue', ascending = True, inplace = True)
        if file not in stats_dict:
            stats_dict[file] = stats_frame
    out_p = os.path.join(path, 'data', 'enrichment_stats_dictionary.pkl')

    write_pickle(stats_dict, out_p)
    print_loc(stats_dict)


def print_loc(stats_dict):
    for k in stats_dict:
        if '_' in k:
            j = k.split('_')[0]
            print(k)
            print(stats_dict[k].index.get_loc(j))
        else:
            print(k)
            print(stats_dict[k].index.get_loc(k))


def join_regulon_frames(frame1, frame2):
    intx = list(set(frame1.Target).intersection(frame2.Target))

    f1 = frame1.set_index('Target')
    f1.columns = ['F1_moa','F1_lh']

    f2 = frame2.set_index('Target')
    f2.columns = ['F2_moa','F2_lh']

    joined_frames = f1.join(f2)
    filt_joined = joined_frames.loc[intx]
    return filt_joined


def structure_dir():
    dirs = ['BCL6_ly7', 'BCL6_pfeiffer', 'FOXM1', 'MEF2B', 'MYB', 'STAT3']
    path = os.getcwd()
    for dir in dirs:
        out_path = os.path.join(path, dir, 'data')
        if not os.path.exists(out_path):
            os.makedirs(out_path)


def restructure_dir(fh_out):
    import shutil
    dirs = ['BCL6_ly7', 'BCL6_pfeiffer', 'FOXM1', 'MEF2B', 'MYB', 'STAT3']
    path = os.getcwd()
    out_path = os.path.join(path,'Bcell',fh_out)
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    for d in dirs:
        fh = os.path.join(path,d)
        shutil.move(fh, out_path)


def fraction_random_regulon(regulon, fraction, random_state=1):
    """ Fraction of total regulon to shuffle to determine Accuracy by fraction of randomized regulon
    Accuracy : ((xn - xr) / xn) *100.0  - where xn is the total number of regulators evaluated and
    xr is the ranked position of the siRNA regulator. Accuracy percent reported

    Args:
        regulon:
        fraction:
        random_state:

    Returns:

    """
    copy_regul = regulon.reset_index()
    copy_regul.set_index(['Regulator', 'Target'],inplace=True)
    rand_regul = mi_regulon.sample(frac = fraction, random_state=random_state)
    sample_regul = rand_regul.sample(frac = 1, random_state=random_state)
    sample_regul.index = rand_regul.index
    copy_regul.update(sample_regul)
    copy_regul.reset_index().set_index('Regulator',inplace=True)

    return copy_regul


def fraction_regulon_size(regulon, fraction, random_state=1):
    """ Decimate regulon based on percentage, ensure that the regulon is of size x, decimate from 100% to the equivalent
     of size x in %.
    This functions will be utilized in conjunction with the subset regulon

    Args:
        regulon:
        fraction:
        random_state:

    Returns:

    """
    copy_regul = regulon.sample(frac = fraction, random_state=random_state)

    return copy_regul


def gaussian_noise(expr, _std):
    """ Adding random gaussian noise centered at zero with a set standard deviation

    Args:
        expr:
        _std:

    Returns:

    """
    noise = np.random.normal(0, _std, expr.shape)
    joined_w_noise = expr + noise

    return joined_w_noise


def fraction_signature_size(expr, fraction, random_state=1):
    """ Randomly removes genes from signature

    Args:
        expr:
        fraction:
        random_state:

    Returns:

    """

    copy_expr = expr.sample(frac = fraction, random_state=random_state)

    return copy_expr