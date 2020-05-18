"""Unit tests for regulon enrichment and regulon_utils modules.

See Also:

    :class:`..features.expression_utils`:


Author: Joey Estabrook <estabroj@ohsu.edu>
"""

import os
import sys
import unittest

base_dir = os.path.dirname(__file__)
data_dir = os.path.join(base_dir, "resources")
sys.path.extend([os.path.join(base_dir, '../..')])


from enricher.enrich import expression_utils

def load_test_expr(expr='test_expr.tsv'):
    return pd.read_csv(os.path.join(data_dir,expr), index_col=0, sep = '\t')


class ExpressionUtilsTestCase(unittest.TestCase):

    def test_load_test_expr(self):
        expr = load_test_expr().T
        self.assertSequenceEqual(expr.shape, (6, 16297))

    def test_log_norm(self):
        expr = load_test_expr().T
        log_expr = expression_utils.log_norm(expr)
        self.assertAlmostEqual(log_expr['A1BG'].mean(), 6.845200696748805)

    def test_fit_and_transform(self):
        expr = load_test_expr().T
        normed_expr = expression_utils.fit_and_transform_array(expr)
        self.assertSequenceEqual(normed_expr.shape,(6, 16258))
        self.assertAlmostEqual(expr_normed.mean().mean(), -0.13697546259217588)

if __name__ == '__main__':
    unittest.main()