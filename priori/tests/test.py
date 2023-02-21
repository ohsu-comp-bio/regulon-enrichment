"""Unit tests for priori class.

See Also:

    :class:`..enrich`:

Author: Joey Estabrook <estabroj@ohsu.edu> & Will Yashar <yashar@ohsu.edu>
"""

import os
import sys
base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, "resources")
sys.path.extend([os.path.join(base_dir, '../..')])

from sklearn.utils.validation import check_array
from priori import priori
import priori.regulon.regulon_enrichment as regulon_enrichment
import priori.regulon.regulon_utils as regulon_utils

from tqdm import tqdm
import warnings
import unittest
import pandas as pd
import numpy as np
import scipy.stats as st
import functools

def load_test_network(net='test.pkl'):
    return regulon_utils.read_pickle(os.path.join(data_dir, net))

def load_test_expr(expr='test_expr.tsv'):
    return pd.read_csv(os.path.join(data_dir, expr), index_col=0, sep = '\t')

class PrioriTestCase(unittest.TestCase):

    def test_enrichment(self):

        expr = load_test_expr()
        net = load_test_network()

        self.assertSequenceEqual(net.shape, (122785, 3))
        self.assertSequenceEqual(expr.shape, (8723, 6))

        enr = priori.Priori(expr=expr, regulon=net)

        self.assertSequenceEqual(enr.expr.shape, expr.shape)
        self.assertSequenceEqual(enr.regulon.shape, net.shape)
        self.assertEqual(enr.scaled, False)
        self.assertEqual(enr.regulators, None)
        self.assertEqual(enr.regulon_size, 15)
        self.assertEqual(enr.regulon_weights, None)
        self.assertEqual(enr.thresh_filter, 0.1)
        self.assertEqual(enr.enrichment, None)

        enr.scale()
        self.assertEqual(enr.scaled, True)

        enr.assign_weights()
        self.assertSequenceEqual(enr.regulon_weights.shape, (26146, 3))
        self.assertAlmostEqual(enr.regulon_weights.MoA.mean(), 2.4612614776160298)

        enr.calculate_enrichment()
        self.assertSequenceEqual(enr.regulators.tolist(), ['AKT1', 'APP', 'AR', 'ARNT', 'ATF1', 'ATF2', 'ATF4', 'BACH1', 'BPTF', 'BRCA1', 'CARM1', 'CAT', 'CBFA2T3', 'CEBPD', 'CHD9', 'CREB1', 'CREBBP', 'CUX1', 'DBP', 'DDIT3', 'E2F3', 'E4F1', 'ELF1', 'ELF2', 'ELK1', 'EP300', 'ESRRA', 'ETS2', 'FOXM1', 'FOXO3', 'GABPB1', 'GATA2', 'GFI1', 'GTF2A2', 'GTF3A', 'HIF1A', 'HMGA1', 'HSF1', 'HSF2', 'IGF1R', 'IRF1', 'IRF2', 'IRF7', 'IRF8', 'LMO2', 'MAFG', 'MAPK1', 'MAX', 'MAZ', 'MED1', 'MIF', 'MYB', 'MYC', 'NCOA1', 'NCOA2', 'NCOA6', 'NF1', 'NFE2L1', 'NFE2L2', 'NR1H3', 'NR3C1', 'NRF1', 'PAX8', 'PCBP1', 'POU2F1', 'PPARA', 'RB1', 'RELA', 'REPIN1', 'RFX1', 'RREB1', 'RUNX1', 'RUNX2', 'RXRB', 'SF1', 'SIRT1', 'SMAD1', 'SMAD3', 'SMAD4', 'SMARCD3', 'SOD2', 'SP1', 'SP3', 'SPI1', 'SREBF1', 'STAT1', 'STAT2', 'STAT3', 'STAT5A', 'STAT5B', 'STAT6', 'TBL1X', 'TBL1XR1', 'TBP', 'TCF12', 'TCF3', 'TFAP2A', 'TFAP4', 'TFCP2', 'TFDP1', 'TFDP2', 'TGFB1', 'TGIF1', 'TGS1', 'TP53', 'TRAF4', 'UBP1', 'XBP1', 'YY1', 'ZBTB14'])
        self.assertSequenceEqual(enr.enrichment.shape, (6, 110))


if __name__ == '__main__':
    unittest.main()
