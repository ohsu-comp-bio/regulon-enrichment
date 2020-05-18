from sklearn.utils.validation import check_array
from .features.expression_utils import log_norm
import warnings
import enricher.regulon.regulon_enrichment as regulon_enrichment
import enricher.features.expression_utils as expression_utils
import enricher.regulon.regulon_utils as regulon_utils
from tqdm import tqdm
import os
import pandas as pd
import numpy as np
import scipy.stats as st
import functools

if __name__ == '__main__':
    DATA_PATH = os.path.join(os.getcwd(), 'data')
else:
    dirname = os.path.dirname(__file__)

    DATA_PATH = os.path.join(dirname, 'data')

sif_file = DATA_PATH + '/PathwayCommons9.All.hgnc.sif.gz'
sec_intx_file = DATA_PATH + '/secondary_intx_regulon.pkl'


def score_enrichment(regulator, expr, regulon, quant_nes):
    """ Function to subset and generate regulator activity scores based on rank ordering of up-regulated
     and down-regulated targets

    Args:
        regulator (str) : Regulator to subset expression frame by
        expr (:obj: `pandas DataFrame`): pandas DataFrame of shape [n_samps, n_feats]
        regulon (:obj: `pandas DataFrame`): pandas DataFrame of regulon returned by compile_regulon
            with columns ['Target', 'MoA', 'likelihood']
        quant_nes (obj: `pandas DataFrame`): quantile enrichment scores for regulators
    Return:
        enrichment_score (:obj: `pandas DataFrame`): pandas DataFrame of activity scores for specified regulator

    """

    down_reg_sub, up_reg_sub = subset_regulon(regulator, regulon, expr)

    # Rank up and down regulated targets by z-scores. Sum rank values across rows
    # (Compute numerical data ranks [1 through n] along axis)
    # and sort samples lowest to highest summed rank score.

    down_reg_ordered = rank_and_order_total(down_reg_sub, regulator, regulon, ascending = False, expr = expr)
    up_reg_ordered = rank_and_order_total(up_reg_sub, regulator, regulon, ascending = True, expr = expr)

    zframe = format_nes_frame(down_reg_ordered, up_reg_ordered, regulator)
    delta = format_delta(down_reg_sub, up_reg_sub)

    zframe[regulator] = zframe.values + delta.values

    enrichment_score = zframe[regulator] + quant_nes.loc[regulator]

    return enrichment_score


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


class Error(Exception):
    """Base class for other exceptions"""
    pass


class OmicError(Error):
    """Raised when duplications in omic features or samples are detected"""
    pass


class Enrichment(object):
    """Base enrichment class for predicting regulon enrichment from -omic datasets.

    Args:
        cohort :
        expr (:obj:`pd.DataFrame`, shape = [n_feats, n_samps])
        regulon (:obj: `pandas DataFrame`)
        regulon_size (int): Minimum number of edges for a given regulator.
        sec_intx_file (str): Path to pre-compiled secondary interaction network.

    """

    def __init__(self, cohort, expr, regulon = None, regulon_size = 15, sec_intx = sec_intx_file, thresh_filter=0.1):
        if not isinstance(expr, pd.DataFrame):
            raise TypeError("`expr` must be a pandas DataFrame, found "
                            "{} instead!".format(type(expr)))

        if len(set(expr.index)) != expr.shape[0]:
            print(len(set(expr.index)))
            print(expr.shape)
            raise OmicError("Duplicate feature names in {cohort} dataset!".format(cohort = cohort))

        if len(set(expr.columns)) != expr.shape[1]:
            raise OmicError("Duplicate sample names in {cohort} dataset!".format(cohort = cohort))

        self.cohort = cohort
        self.expr = expr

        if regulon is None:
            self.regulon = regulon_utils.read_pickle(sec_intx)

        else:
            self.regulon = regulon

        self.scaler_type = None
        self.scaled = False
        self.regulon_size = regulon_size
        self.regulon_weights = None
        self.thresh_filter = thresh_filter
        self.total_enrichment = None
        self.regulators = None
        self.quant_nes = None

    def __str__(self):
        return """------\nCohort: {}\nn-features: {}\nn-samples: {}\nscaler: {}\nscaled: {}\nregulon threshold: {}\
            \nregulon nodes: {}\nregulon edges: {}\n------\n""".format(self.cohort, self.expr.shape[0],
                self.expr.shape[1], self.scaler_type, self.scaled,self.regulon_size, len(self.regulon.UpGene.unique()),
                self.regulon.shape[0])

    def __repr__(self):
        return """------\nCohort: {}\nn-features: {}\nn-samples: {}\nscaler: {}\nscaled: {}\nregulon threshold: {}\
        \nregulon nodes: {}\nregulon edges: {}\n------\n""".\
            format(self.cohort, self.expr.shape[0], self.expr.shape[1], self.scaler_type, self.scaled,
                   self.regulon_size, len(self.regulon.UpGene.unique()), self.regulon.shape[0])

    @staticmethod
    def _preprocess_data(expr, scaler_type = 'robust', thresh_filter = 0.1):
        """ Centers expression data based on a specified data scaler algorithm

        Args:
            expr (pandas DataFrame obj): pandas DataFrame of [n_features, n_samples]
            scaler_type (str): Scaler to normalized features/samples by: standard | robust | minmax | quant
            thresh_filter (float): Prior to normalization remove features that have a standard deviation per feature
            less than {thresh_filter}

        Returns:
            scaled_frame (:obj: `pandas DataFrame`) : pandas DataFrame containing scaled expression data of
                shape [n_samples, n_features]

        """

        # By default, the input is checked to be a non-empty 2D array containing
        # only finite values.
        _ = check_array(expr)

        scaler_opt = {'standard': expression_utils.StandardScaler(), 'robust': expression_utils.RobustScaler(),
                      'minmax': expression_utils.MinMaxScaler(), 'quant': expression_utils.QuantileTransformer()}

        if scaler_type not in scaler_opt:
            raise KeyError('{scaler_type} not supported scaler_type! Supported types include: {keys}'.format(
                scaler_type = scaler_type, keys = ' | '.join(scaler_opt.keys())))

        scaler = scaler_opt[scaler_type]

        # Transpose frame to correctly orient frame for scaling and machine learning algorithms
        print('--- log2 normalization ---')

        expr_t = expr[(expr.std(axis = 1) > thresh_filter)].T
        expr_lt = expression_utils.log_norm(expr_t)

        print('--- Centering features with {} scaler ---'.format(scaler_type))
        scaled_frame = pd.DataFrame(scaler.fit_transform(expr_lt), index = expr_lt.index, columns = expr_lt.columns)

        return scaled_frame

    @staticmethod
    def _prune_regulon(expr, regulon, regulon_size):
        """ Prunes regulon with secondary interactions that do not meet the necessary number of downstream interactions
        metric {regulon_size}

        Args:
            expr (pandas DataFrame obj): pandas DataFrame of [n_samples, n_features]
            regulon (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions between regulator and
                downstream members of its regulon of shape [len(Target), ['Regulator','Target','MoA','likelihood']
            regulon_size (int) : number of downstream interactions required for a given regulator in order to calculate
                enrichment score

        Returns:
            filtered_regulon (:obj: `pandas DataFrame`) : pandas DataFrame containing weight interactions
                between regulator and downstream members of its regulon of shape :
                [len(Target), ['Regulator','Target','MoA','likelihood']

        """

        expr_filtered_regulon = regulon[
            ((regulon.UpGene.isin(expr.columns)) & (regulon.DownGene.isin(expr.columns)))].set_index('UpGene')
        idx = (expr_filtered_regulon.index.value_counts() >= regulon_size)

        filtered_regulon = expr_filtered_regulon.loc[idx[idx == True].index].reset_index()

        return filtered_regulon

    @staticmethod
    def _structure_weights(regulator, pruned_regulon, f_statistics, r_frame, p_frame):
        """ Calculates weights associated with regulators. Weights are the summation of the F-statistic and
            absolute spearman correlation coefficient. The weight retains the sign of the spearman
            correlation coefficient.

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

    def scale(self, scaler_type = 'robust', thresh_filter = 0.1):
        self.scaler_type = scaler_type
        self.expr = self._preprocess_data(self.expr, self.scaler_type, thresh_filter)
        self.scaled = True

    def assign_weights(self):
        if not self.scaled:
            warnings.warn('Assigning interaction weights without scaling dataset!')

        pruned_regulon = self._prune_regulon(self.expr, self.regulon, self.regulon_size)
        # noinspection PyTypeChecker
        r, p = regulon_utils.spearmanr(self.expr)

        r_frame = pd.DataFrame(r, columns = self.expr.columns, index = self.expr.columns)
        p_frame = pd.DataFrame(p, columns = self.expr.columns, index = self.expr.columns)

        F_statistics = {regulator: regulon_utils.f_regression(self.expr.reindex(frame.DownGene, axis = 1),
                        self.expr.reindex([regulator], axis = 1).values.ravel())
                        for regulator, frame in pruned_regulon.groupby('UpGene')}

        weights = pd.concat([self._structure_weights(regulator, pruned_regulon, F_statistics, r_frame, p_frame)
                            for regulator in F_statistics])

        self.regulon_weights = weights[~np.isinf(weights.MoA)]

    def calculate_enrichment(self):
        if self.regulon_weights is None:
            raise TypeError("`regulon_weights` must be assigned prior to enrichment calculation, found "
                            "{} instead!".format(type(self.regulon_weights)))

        quant_nes = regulon_enrichment.quantile_nes_score(self.regulon_weights, self.expr.T)
        self.quant_nes = quant_nes
        self.regulators = self.regulon_weights.index.unique()

        print('--- Calculating regulon enrichment scores ---')
        nes_list = list(map(functools.partial(regulon_enrichment.score_enrichment, expr =self.expr,
                        regulon = self.regulon_weights, quant_nes =quant_nes), tqdm(self.regulators)))

        self.total_enrichment = pd.concat(nes_list, axis =1)
