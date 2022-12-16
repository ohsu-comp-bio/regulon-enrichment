"""Unit tests for enrich class.

See Also:

    :class:`..enrich`:

Author: Joey Estabrook <estabroj@ohsu.edu>
"""

import os
import sys
base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, "tests/resources")
#sys.path.extend([os.path.join(base_dir, '../..')])

from sklearn.utils.validation import check_array
from priori import Priori
import regulon.regulon_enrichment as regulon_enrichment
import features.expression_utils as expression_utils
import regulon.regulon_utils as regulon_utils
from features.expression_utils import log_norm

from tqdm import tqdm
import warnings
import unittest
import pandas as pd
import numpy as np
import scipy.stats as st
import functools




def load_test_sif(sif='test.sif'):
    return pd.read_table(os.path.join(data_dir, sif), index_col=0)


def load_test_expr(expr='test_expr.tsv'):
    return pd.read_csv(os.path.join(data_dir, expr), index_col=0, sep = '\t')


class EnrichTestCase(unittest.TestCase):
    def test_load_test_sif(self):
        sif = load_test_sif()
        self.assertSequenceEqual(sif.shape, (1302, 3))

    def test_load_test_expr(self):
        expr = load_test_expr()
        self.assertSequenceEqual(expr.shape, (8723, 6))

    def test_enrichment(self):
        sif = load_test_sif()
        expr = load_test_expr()
        filt_sif = regulon_utils.filter_sif(sif)
        enr = Priori(expr=expr, regulon=filt_sif)
        self.assertSequenceEqual(enr.expr.shape, expr.shape)
        self.assertSequenceEqual(enr.regulon.shape, filt_sif.shape)
        self.assertEqual(enr.scaled, False)
        self.assertEqual(enr.regulators, None)
        self.assertEqual(enr.regulon_size, 15)
        self.assertEqual(enr.regulon_weights, None)
        self.assertEqual(enr.thresh_filter, 0.1)
        self.assertEqual(enr.total_enrichment, None)
        self.assertEqual(enr.quant_nes, None)

        enr.scale()
        self.assertEqual(enr.scaled, True)

        enr.assign_weights()
        self.assertSequenceEqual(enr.regulon_weights.shape, (433, 3))
        self.assertAlmostEqual(enr.regulon_weights.MoA.mean(), 1.1555032640512617)

        enr.calculate_enrichment()
        self.assertSequenceEqual(enr.regulators.tolist(), ['TP53'])
        self.assertSequenceEqual(enr.total_enrichment.shape, (6, 1))


if __name__ == '__main__':
    unittest.main()
