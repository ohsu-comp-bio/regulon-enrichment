"""Unit tests for regulon enrichment and regulon_utils modules.

See Also:

    :class:`..regulon.regulon_enrichment`:
    :class:`..regulon.regulon_utils`:

Author: Joey Estabrook <estabroj@ohsu.edu>
"""

import os
import sys
import unittest
import warnings
import pandas as pd
from enricher.enrich import regulon_utils, regulon_enrichment
warnings.simplefilter("ignore", UserWarning)

base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, "resources")
sys.path.extend([os.path.join(base_dir, '../..')])


def load_test_sif(sif='test.sif'):
    return pd.read_table(os.path.join(data_dir, sif), index_col=0)


def load_test_expr(expr='test_expr.tsv'):
    return pd.read_csv(os.path.join(data_dir, expr), index_col=0, sep = '\t')


class RegulonUtilsTestCase(unittest.TestCase):
    def test_load_test_sif(self):
        sif = load_test_sif()
        self.assertSequenceEqual(sif.shape, (1302, 3))

    def test_filter_sif(self):
        sif = load_test_sif()
        self.assertEqual(regulon_utils.filter_sif(sif).shape[0], 800)
        self.assertEqual(regulon_utils.filter_sif(sif)['Type'].unique().tolist()[0], 'controls-expression-of')

    def test_load_test_expr(self):
        expr = load_test_expr().T
        self.assertSequenceEqual(expr.shape, (6, 16297))

    def test_prune_regulon(self):
        sif = load_test_sif()
        expr = load_test_expr().T
        filt_sif = regulon_utils.filter_sif(sif)
        filtered_regulon = regulon_utils.prune_regulon(expr, filt_sif, 15)
        self.assertSequenceEqual(filtered_regulon.shape, (606, 3))

    def test_regulon_weight_assignment(self):
        sif = load_test_sif()
        expr = load_test_expr().T
        filt_sif = regulon_utils.filter_sif(sif)
        filtered_regulon = regulon_utils.prune_regulon(expr, filt_sif, 15)
        regul_weights = regulon_utils.regulon_weight_assignment('TP53', expr, filtered_regulon)
        self.assertSequenceEqual(regul_weights.shape, (606, 3))
        self.assertSequenceEqual(regul_weights.columns.tolist(), ['Target', 'MoA', 'likelihood'])
        self.assertEqual(regul_weights.iloc[0, :].tolist()[0], 'AARS')
        self.assertAlmostEqual(regul_weights.iloc[0, :].tolist()[1], 0.1724812122096268)
        self.assertAlmostEqual(regul_weights.iloc[0, :].tolist()[2], 0.8717434402332361)

    def test_quantile_nes_score(self):
        sif = load_test_sif()
        expr = load_test_expr().T
        filt_sif = regulon_utils.filter_sif(sif)
        filtered_regulon = regulon_utils.prune_regulon(expr, filt_sif, 15)
        regul_weights = regulon_utils.regulon_weight_assignment('TP53', expr, filtered_regulon)
        nes = regulon_enrichment.quantile_nes_score(regul_weights, expr.T)
        self.assertSequenceEqual(nes.columns.tolist(),
                                 ['Test_A1', 'Test_A2', 'Test_A3', 'Test_D1', 'Test_D2', 'Test_D3'])
        self.assertAlmostEqual(nes.values.mean(), -0.8396900886423703)

    def test_score_enrichment(self):
        sif = load_test_sif()
        expr = load_test_expr().T
        filt_sif = regulon_utils.filter_sif(sif)
        filtered_regulon = regulon_utils.prune_regulon(expr, filt_sif, 15)
        regul_weights = regulon_utils.regulon_weight_assignment('TP53', expr, filtered_regulon)
        nes = regulon_enrichment.quantile_nes_score(regul_weights, expr.T)
        enrichment = regulon_enrichment.score_enrichment('TP53', expr, regul_weights, nes)
        self.assertAlmostEqual(enrichment.values.mean(), -0.7888152886796247)


if __name__ == '__main__':
    unittest.main()
