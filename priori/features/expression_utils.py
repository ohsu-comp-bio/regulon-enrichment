import numpy as np
import os
from sklearn.preprocessing import RobustScaler, MinMaxScaler, StandardScaler, QuantileTransformer
import pandas as pd
import re


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

def format_expr(df):
    """
    Re-format the input expression matrix so the gene names are stored in the index.

    """
    if 'features' in df.columns:
        df = df.set_index('features') # If "features" column already exists, set index as features
        return df

    # If the features column doesn't exist, assign features from the rownames
    else:
        print("The 'features' column wasn't found. Assigning features from the rownames.")
        df = df.rename_axis('features')
        return df

def check_for_duplicates(df):
    """
    Check for duplicate rows or columns in the input expression matrix.

    """
    # Check for duplicate rows
    if df.index.duplicated().any():
        dup_features = df.index[df.index.duplicated()]
        raise ValueError(f"Duplicate genes ({', '.join(map(str, dup_features))}) were found in your input expression matrix. Please re-format your input expression matrix.")

def check_for_ensembl_id_format(df):
    """
    Check if at least one gene name in the list follows the Ensembl ID format.

    """
    ensembl_id_pattern = re.compile(r'^ENSG')

    for name in df.index:
        if ensembl_id_pattern.match(name):
            raise ValueError(f"Error: '{name}' is in Ensembl ID format. Please provide features as gene symbols.")

def check_expr(df):
    """
    Re-format the input expression matrix so the gene names are stored in a column called "features".

    Args:
        - df (data frame): the input expression matrix

    Returns:
        - df or None: Returns an error message if the input matrix does not have gene names stored 
            in a column called "features" or in the rownames.
    """

    # Format the expression matrix
    df = format_expr(df)

    # Check for duplicate features or sample names
    check_for_duplicates(df)

    # Verify that genes are listed as gene symbols, not as ensemble IDs
    check_for_ensembl_id_format(df)

    return df


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
