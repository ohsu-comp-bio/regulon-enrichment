import numpy as np
import os
from sklearn.preprocessing import RobustScaler, MinMaxScaler, StandardScaler, QuantileTransformer
import pandas as pd


def load_expr(expr_f):
    """

    Args:
        expr_f (str): absolute path to tab delimited expression file of shape = [n_features, n_samples]

    Returns:
        expr (:obj:`np.array` of :obj:`float`,
          shape = [n_samples, n_features])

    """
    expr = pd.read_csv(expr_f, sep='\t', index_col=0)

    return expr


def log_norm(expr):
    """Log-normalizes a dataset, usually RNA-seq expression.
    Puts a matrix of continuous values into log-space after adding
    a constant derived from the smallest non-zero value.
    Args:
        expr (:obj:`np.array` of :obj:`float`,
                  shape = [n_samples, n_features])
    Returns:
        expr (:obj:`pandas DataFrame`): shape = [n_features, n_samples]

    Examples:
        >>> norm_expr = log_norm(np.array([[1.0, 0], [2.0, 8.0]]))
        >>> print(norm_expr)
                [[ 0.5849625 , -1.],
                 [ 1.32192809,  3.08746284]]

    """

    log_add = np.nanmin(expr[expr > 0]) * 0.5
    norm_mat = np.log2(expr + log_add)
    norm_mat = norm_mat.dropna(axis=1)

    return norm_mat


def fit_and_transform_array(expr, norm_type='robust', feature = True, sample = False, thresh_filter = 0.4, scale=True):
    """ Fit and scale expression data based on a specified data scaler algorithm

    Args:
        expr (pandas DataFrame obj): pandas DataFrame of [n_features, n_samples]
        norm_type (str): Scaler to normalized features/samples by: standard | robust | minmax | quant
        feature (bool): Scale expression data by features
        sample (bool): Scale expression data by both features and samples
        thresh_filter (float): Prior to normalization remove features that do not have the mean unit of
            a feature (i.e. 1 tpm) is greater than {thresh_filter}
        scale (bool): optional arg to avoid scaling dataset if data set has been normalized prior to analysis

    Returns:
        scaled_frame (:obj: `pandas DataFrame`) : pandas DataFrame containing scaled expression data of
            shape [n_samples, n_features]

    """
    scaler_opt = {'standard': StandardScaler(), 'robust': RobustScaler(), 'minmax': MinMaxScaler(),
                  'quant': QuantileTransformer()}
    print('--- setting {} as scaler to normalize features|samples by ---'.format(norm_type))
    scaler = scaler_opt[norm_type]

    if scale:
        # Transpose frame to correctly orient frame for scaling and machine learning algorithms
        expr = expr.groupby(expr.index).mean()
        expr = expr[(expr.mean(axis=1) > thresh_filter)].T
        expr = log_norm(expr)
        print('--- log2 normalization ---')
        if feature and sample:
            # scale by both feature and sample
            print('--- scaling by feature and sample ---')

            scaler_s = scaler_opt[norm_type]
            f_scaled_expr = pd.DataFrame(scaler.fit_transform(expr), index = expr.index, columns = expr.columns).T
            scaled_frame = pd.DataFrame(scaler_s.fit_transform(f_scaled_expr), index = f_scaled_expr.index,
                                        columns = f_scaled_expr.columns).T
        elif feature and not sample:
            # scale by feature
            print('--- scaling by feature ---')

            scaled_frame = pd.DataFrame(scaler.fit_transform(expr), index = expr.index, columns = expr.columns)

        else:
            print('--- expression dataset will not be scaled ---')
            scaled_frame = expr
    else:
        print('--- expression dataset will not be scaled ---')
        scaled_frame = expr

    return scaled_frame


def load_scaled_expr(expr, cohort, norm_type='robust', feature = True, sample = False, thresh_filter = 0.4, scale=True):
    """ Checks if expression file has been normalized, if not fits and scales expression data based on a specified data
        scaler algorithm, else loads the pickled object
    Args:
        expr (pandas DataFrame obj): pandas DataFrame of [n_features, n_samples]
        cohort (str) : name of cohort to associate with compiled regulon
        norm_type (str): Scaler to normalized features/samples by: standard | robust | minmax | quant
        feature (bool): Scale expression data by features
        sample (bool): Scale expression data by both features and samples
        thresh_filter (float): Prior to normalization remove features that do not have the mean unit of
            a feature (i.e. 1 tpm) is greater than {thresh_filter}
        scale (bool): optional arg to avoid scaling dataset if data set has been normalized prior to analysis

    Returns:
        scaled_frame (:obj: `pandas DataFrame`) : pandas DataFrame containing scaled expression data of
            shape [n_samples, n_features]
    """
    scaled_expr = os.path.join(dirname, '../experiments/{cohort}/data/{cohort}_{norm_type}_{feature}_{sample}_{thresh_filter}_{scale}_frame.pkl'.
                               format(cohort = cohort, norm_type = norm_type, feature = feature, sample = sample, thresh_filter = thresh_filter, scale = scale))

    if os.path.isfile(scaled_expr):
        print('--- Loading scaled expression data ---')
        nes = read_pickle(scaled_expr)
    else:
        print('--- Generating scaled expression data ---')
        nes = fit_and_transform_array(expr = expr, norm_type=norm_type, feature = feature,
                                      sample = sample, thresh_filter = thresh_filter, scale = scale)
        write_pickle(nes, scaled_expr)

    return nes
