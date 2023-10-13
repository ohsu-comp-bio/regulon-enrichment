import warnings 
import os
import functools
import pandas as pd
import numpy as np
import statsmodels.api as sm
from tqdm import tqdm
from sklearn.utils.validation import check_array
import priori.regulon.regulon_enrichment as regulon_enrichment
import priori.features.expression_utils as expression_utils
import priori.regulon.regulon_utils as regulon_utils
import argparse
import re
warnings.simplefilter("ignore", UserWarning)

# Set global data path
if __name__ == '__main__':
    DATA_PATH = os.path.join(os.getcwd(), 'data')
else:
    dirname = os.path.dirname(__file__)

    DATA_PATH = os.path.join(dirname, 'data')

# Set path to local Pathway Commons relationships
pathway_commons = DATA_PATH + '/primary_intx_regulon.pkl'

class Error(Exception):
    """Base class for other exceptions"""


class OmicError(Error):
    """Raised when duplications in omic features or samples are detected"""


class Priori(object):
    """Base enrichment class for predicting regulon enrichment from -omic datasets.

    Args:
        expr (:obj:`pd.DataFrame`, shape = [n_feats, n_samps])
        regulon (:obj: `pandas DataFrame`)
        regulon_size (int): Minimum number of edges for a given regulator.

    """

    def __init__(self, expr, regulon=None, regulon_size=15,
                 thresh_filter=0.1):
        if not isinstance(expr, pd.DataFrame):
            raise TypeError("`expr` must be a pandas DataFrame, found "
                            "{} instead!".format(type(expr)))

        if len(set(expr.index)) != expr.shape[0]:
            print(len(set(expr.index)))
            print(expr.shape)
            raise OmicError("Duplicate feature names in the dataset!")

        if len(set(expr.columns)) != expr.shape[1]:
            raise OmicError("Duplicate sample names in the dataset!")

        # Assign object attributes defined by arguments
        self.expr = expr
        self.regulon_size = regulon_size
        self.thresh_filter = thresh_filter

        if regulon is None:
            self.regulon = regulon_utils.read_pickle(pathway_commons)

        else:
            self.regulon = regulon_utils.read_custom_network(regulon)

        # Initialize attributes to be defined
        self.scaler_type = None
        self.scaled = False
        self.regulon_weights = None
        self.enrichment = None
        self.regulators = None

        self.correlations = None
        self.p_values = None
        self.adjusted_p_values = None

    def __str__(self):
        return """------\nn-features: {}\nn-samples: {}\nscaler: {}\nscaled:\
        {}\nregulon threshold: {}\nregulon nodes: {}\nregulon edges: {}\n------\n""".\
            format(self.expr.shape[0],
                   self.expr.shape[1],
                   self.scaler_type,
                   self.scaled, self.regulon_size,
                   len(self.regulon.UpGene.unique()),
                   self.regulon.shape[0])

    def __repr__(self):
        return """------\nn-features: {}\nn-samples: {}\nscaler: {}\nscaled: {}\
        \nregulon threshold: {}\nregulon nodes: {}\nregulon edges: {}\n------\n""".\
            format(self.expr.shape[0],
                   self.expr.shape[1],
                   self.scaler_type,
                   self.scaled,
                   self.regulon_size,
                   len(self.regulon.UpGene.unique()),
                   self.regulon.shape[0])

    @staticmethod
    def _preprocess_data(expr, scaler_type='robust', thresh_filter=0.1):
        """ Centers expression data based on a specified data scaler algorithm

        Args:
            expr (pandas DataFrame obj): pandas DataFrame of [n_features, n_samples]
            scaler_type (str): Scaler to normalized features/samples by:
                standard | robust | minmax | quant
            thresh_filter (float): Prior to normalization remove features that have
                a standard deviation per feature less than {thresh_filter}

        Returns:
            scaled_frame (:obj: `pandas DataFrame`) : pandas DataFrame containing
                scaled expression data of shape [n_samples, n_features]

        """
        # By default, the input is checked to be a non-empty 2D array containing
        # only finite values.
        _ = check_array(expr)

        scaler_opt = {'standard': expression_utils.StandardScaler(),
                      'robust': expression_utils.RobustScaler(),
                      'minmax': expression_utils.MinMaxScaler(),
                      'quant': expression_utils.QuantileTransformer()}

        if scaler_type not in scaler_opt:
            raise KeyError('{scaler_type} not supported scaler_type!'
                           ' Supported types include: {keys}'.format(
                               scaler_type=scaler_type, keys=' | '.join(scaler_opt.keys())))

        scaler = scaler_opt[scaler_type]

        # Transpose frame to correctly orient frame for scaling and machine learning algorithms
        print('--- log2 normalization ---')
        expr_t = expr[(expr.std(axis=1) > thresh_filter)].T

        # Shift expression matrix so values are non-negative
        min_expr = expr_t.min().min()

        if min_expr <= 0:
            expr_t = expr_t - min_expr + 1e-12

        # Log2 transform
        expr_lt = expression_utils.log_norm(expr_t+1e-12)

        print('--- Centering features with {} scaler ---'.format(scaler_type))
        scaled_frame = pd.DataFrame(scaler.fit_transform(expr_lt),
                                    index=expr_lt.index,
                                    columns=expr_lt.columns)

        return scaled_frame

    @staticmethod
    def _prune_regulon(expr, regulon, regulon_size):
        """ Prunes regulon with secondary interactions that do not meet
            the necessary number of downstream interactions metric {regulon_size}

        Args:
            expr (pandas DataFrame obj): pandas DataFrame of [n_samples, n_features]
            regulon (:obj: `pandas DataFrame`) : pandas DataFrame containing weight
                interactions between regulator and downstream members of its regulon
                of shape [len(Target), ['Regulator','Target','MoA','likelihood']
            regulon_size (int) : number of downstream interactions required for a
                given regulator in order to calculate enrichment score

        Returns:
            filtered_regulon (:obj: `pandas DataFrame`) : pandas DataFrame containing weight
                interactions between regulator and downstream members of its regulon of shape :
                [len(Target), ['Regulator','Target','MoA','likelihood']

        """

        # Identify transcription factors that are present in the expression data as well as the prior network
        expr_filtered_regulon = regulon[
            ((regulon.UpGene.isin(expr.columns)) & (regulon.DownGene.isin(expr.columns)))].\
            set_index('UpGene')
        idx = (expr_filtered_regulon.index.value_counts() >= regulon_size)

        # Filter regulon and expression datasets
        filtered_regulon = expr_filtered_regulon.loc[idx[idx == True].index].reset_index()
        filtered_regulon=filtered_regulon.rename(columns={"index":"UpGene"})
        edges = list(set(filtered_regulon.UpGene) | set(filtered_regulon.DownGene))
        sub_expr = expr.loc[:,edges]

        return filtered_regulon, sub_expr

    @staticmethod
    def _structure_weights(regulator, pruned_regulon, f_statistics, r_frame, p_frame):
        """ Calculates weights associated with regulators. Weights are the summation of
            the F-statistic and absolute spearman correlation coefficient. The weight
            retains the sign of the spearman correlation coefficient.

        Args:
            regulator (str): A feature to assign weights to downstream interactions
            pruned_regulon (:obj:`pd.DataFrame`, shape = [n_interactions, 3]
            f_statistics (dict) : Dictionary with key:{regulator} key and
            r_frame (:obj:`pd.DataFrame`), shape = [n_features, n_features]
            p_frame (:obj:`pd.DataFrame`), shape = [n_features, n_features]

        Returns:
            weights_ordered (:obj:`pd.DataFrame`), shape = [n_interactions, 3]

        """

        # Identify the transcription factors and their targets
        sub_regul = pruned_regulon[(pruned_regulon['UpGene'] == regulator)]
        targs = sub_regul.DownGene

        # Format the spearman correlation ranks and p-values as well as the F regression values
        p_ = p_frame.loc[targs, regulator]
        p_.name = 'likelihood'
        f_ = f_statistics[regulator][0]
        r_ = r_frame.loc[targs, regulator]

        # Calculate transcription factor gene weights
        w_ = (f_ * abs(r_)) * np.sign(r_)

        # Reformat weights, assigning the likelihood as the spearman correlation p value
        w_.index.name = 'Target'
        w_.name = 'MoA'
        weights = w_.to_frame()
        weights['likelihood'] = p_
        weights['Regulator'] = regulator
        weights_ordered = weights.reset_index().\
            reindex(['Regulator', 'Target', 'MoA', 'likelihood'], axis=1)\
            .set_index('Regulator')

        return weights_ordered

    def scale(self, scaler_type='robust', thresh_filter=0.1):
        """ Fit and scale expression data based on a specified data scaler algorithm

        Args:
            scaler_type (str): Scaler to normalized features/samples by:
                standard | robust | minmax | quant
            thresh_filter (float): Prior to normalization remove features that do not have
                the mean unit of a feature (i.e. 1 tpm) is greater than {thresh_filter}

        """
        self.scaler_type = scaler_type
        if scaler_type == None:
            warnings.warn('Proceeding without scaling dataset!')
            self.expr = self.expr.T
        else:
            self.expr = self._preprocess_data(self.expr, self.scaler_type, thresh_filter)
            self.scaled = True

    def assign_weights(self):
        """
        Generate normalized likelihood weights and assigns those weights to the absolute gene
            expression signature
        """
        if not self.scaled:
            warnings.warn('Assigning interaction weights without scaling dataset!')

        # Filter regulon and expression by transcription factors present in each
        pruned_regulon, sub_expr = self._prune_regulon(self.expr, self.regulon, self.regulon_size)
        self.expr = sub_expr

        # Calculate spearman correlation
        r, p = regulon_utils.spearmanr(self.expr)

        # FDR-adjust p-values using the Benjamini-Hochberg method
        bool, p_adj, sidak, bonf = sm.stats.multipletests(p.flatten(), method = "fdr_bh")
        p_adj = p_adj.reshape(p.shape)
        
        # Create data frames
        r_frame = pd.DataFrame(r, columns=self.expr.columns, index=self.expr.columns)
        p_frame = pd.DataFrame(p, columns=self.expr.columns, index=self.expr.columns)
        p_adj_frame = pd.DataFrame(p_adj, columns=self.expr.columns, index=self.expr.columns)

        # Calculate F regression to correct for collinearity
        F_statistics = {regulator: regulon_utils.f_regression(
            self.expr.reindex(frame.DownGene, axis=1),
            self.expr.reindex([regulator], axis=1).values.ravel())
                        for regulator, frame in pruned_regulon.groupby('UpGene')}

        # Assign spearman correlation coefficient, adjusted p-values, and F-statistics to weights
        weights = pd.concat([self._structure_weights(regulator,
                                                     pruned_regulon,
                                                     F_statistics,
                                                     r_frame,
                                                     p_adj_frame)
                             for regulator in F_statistics])

        self.regulon_weights = weights[~np.isinf(weights.MoA)]
        self.correlations = r_frame
        self.p_values = p_frame
        self.adjusted_p_values = p_adj_frame

    def calculate_enrichment(self):
        """
        Subset and generate regulator activity scores based on rank ordering of up-regulated
            and down-regulated targets

        """
        if self.regulon_weights is None:
            raise TypeError("`regulon_weights` must be assigned prior to enrichment calculation,"
                            " found {} instead!".format(type(self.regulon_weights)))

        self.regulators = self.regulon_weights.index.unique()

        print('--- Calculating regulon enrichment scores ---')
        enrich_list = list(map(functools.partial(regulon_enrichment.score_enrichment,
                                              expr=self.expr, regulon=self.regulon_weights),
                            tqdm(self.regulators)))

        self.enrichment = pd.concat(enrich_list, axis=1)

def main():
    parser = argparse.ArgumentParser(
        "Infer transcription factor activity from gene expression data utilizing pathway and molecular interactions "
        "and mechanisms available through Pathway Commons."
    )

    parser.add_argument('expr', type=str, help="which tab delimited expression matrix to use "
                                               "shape : [n_features, n_samples]"
                                               "units : TPM, RPKM")
    parser.add_argument('out_dir', type=str, help="output directory")
    parser.add_argument('--regulon', type=str, help="optional regulon containing weight interactions between "
                                                  "regulator and downstream members of its regulon"
                                                  "shape : [len(Target), ['Regulator','Target','MoA','likelihood']",
                                                  default=None)
    parser.add_argument('--regulon_size', type=int, help="number of downstream interactions required for a given "
                                                       "regulator in order to calculate enrichment score", default=15)
    parser.add_argument('--scaler_type', type=str, help="Scaler to normalized features/samples by: "
                                                      "standard | robust | minmax | quant", default='robust')
    parser.add_argument('--thresh_filter', type=float, help="Prior to normalization remove features that have a standard "
                                                        "deviation per feature less than {thresh_filter}",
                                                        default=0.1)
    # Parse command line arguments
    args = parser.parse_args()

    # Read and check expression matrix
    expr_matrix = pd.read_table(args.expr)
    expr_matrix = expression_utils.check_expr(expr_matrix)
    print('Loaded input expression data from {}'.format(args.expr))

    # Create Priori object
    enr_obj = Priori(expr=expr_matrix, 
        regulon=args.regulon,
        regulon_size=args.regulon_size, 
        thresh_filter=args.thresh_filter)

    print(enr_obj)
    print('\nScaling data...\n')
    enr_obj.scale(scaler_type=args.scaler_type, thresh_filter=args.thresh_filter)
    print('\nData scaled!\n')

    print('\nAssigning weights...\n')
    enr_obj.assign_weights()
    print('\nWeights assigned!\n')

    print('\nCalculating enrichment...\n')
    enr_obj.calculate_enrichment()
    print('\nPriori activity scores calculated!\n')

    regulon_utils.ensure_dir(args.out_dir)
    regulon_utils.write_pickle(enr_obj, os.path.join(args.out_dir,'priori_object.pkl'))
    enr_obj.enrichment.to_csv(os.path.join(args.out_dir,'priori_activity_scores.tsv'),sep='\t')
    enr_obj.regulon_weights.to_csv(os.path.join(args.out_dir,'priori_activity_score_weights.tsv'),sep='\t')
    print('Complete')

if __name__ == "__main__":
    main()