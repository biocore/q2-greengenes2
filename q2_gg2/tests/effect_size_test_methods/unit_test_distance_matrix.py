import unittest
import pandas as pd
from skbio import DistanceMatrix
import numpy as np
from q2_gg2._methods import (_split_distance_matrix)


class UtilTests(unittest.TestCase):

    def setUp(self):
        self.data = np.array([[0.0, 0.5, .75, .2, 1.0],
                              [0.5, 0.0, .4, .5, .2],
                              [.75, .4, 0.0, .9, .2],
                              [.2, .5, .9, 0.0, .25],
                              [1.0, .2, .2, .25, 0.0]])
        self.ids = ["A1", "A2", "A3", "A4", "A5"]
        self.dm = DistanceMatrix(self.data, self.ids)
        self.col_names = ['paired_sample', 'preparation']
        self.ind = pd.Index(['A1', 'A2', 'A3', 'A4', 'A5'])

    def test_split_distance_matrix_strict(self):

        table = pd.DataFrame([['A2', '16S'], ['A1', 'WGS'],
                              ['A5', '16S'], ['', '16S'],
                              ['A3', 'WGS']],
                             index=self.ind,
                             columns=self.col_names)
        expected_one = self.dm.filter(["A1", "A3"])
        expected_two = self.dm.filter(["A2", "A5"])
        self.assertEqual(_split_distance_matrix(self.dm,
                                                table,
                                                strict=True),
                                               (expected_one, expected_two))

    def test_split_distance_matrix_not_strict(self):

        table = pd.DataFrame([['A2', '16S'], ['A1', 'WGS'],
                              ['A5', '16S'], ['', '16S'],
                              ['A3', 'WGS']],
                             index=self.ind,
                             columns=self.col_names)
        expected_one = self.dm.filter(["A1", "A3", "A4"])
        expected_two = self.dm.filter(["A2", "A5"])
        self.assertEqual(_split_distance_matrix(self.dm, table),
                                               (expected_one, expected_two))

    def test_split_distance_matrix_raises_preparation(self):

        col_names = ['paired_sample']
        table = pd.DataFrame([['A2'], ['A1'], ['A5'],
                              [''], ['A3']], index=self.ind,
                             columns=col_names)
        with self.assertRaisesRegex(KeyError,
                                    'Metadata must include preparation'):
            _split_distance_matrix(self.dm, table)

    def test_split_distance_matrix_raises_paired_sample(self):

        col_names = ['preparation']
        table = pd.DataFrame([['16S'], ['WGS'], ['16S'],
                              ['16S'], ['WGS']], index=self.ind,
                             columns=col_names)
        with self.assertRaisesRegex(KeyError,
                                    'Metadata must include paired_sample'):
            _split_distance_matrix(self.dm, table)

    def test_split_distance_matrix_raises_missing_preparation(self):

        table = pd.DataFrame([['A2', '16S'], ['A1', 'WGS'], ['A5', None],
                              ['', '16S'], ['A3', 'WGS']], index=self.ind,
                             columns=self.col_names)
        with self.assertRaisesRegex(ValueError, 'Missing preparations'):
            _split_distance_matrix(self.dm, table)

    def test_split_distance_matrix_raises_too_many_preparations(self):

        table = pd.DataFrame([['A2', '16S'], ['A1', 'WGS'], ['A5', '16S'],
                              [None, '16S'], ['A3', 'RR']], index=self.ind,
                             columns=self.col_names)
        with self.assertRaisesRegex(ValueError,
                                    'Must have exactly 2 unique prepartions'):
            _split_distance_matrix(self.dm, table)

    def test_split_distance_matrix_raises_missing_ids(self):
        self.ind = pd.Index(['A1', 'A2', 'A3', 'A4', 'A6'])
        table = pd.DataFrame([['A2', '16S'], ['A1', 'WGS'],
                             ['A5', '16S'], [None, '16S'],
                             ['A3', 'WGS']], index=self.ind,
                             columns=self.col_names)
        with self.assertRaisesRegex(ValueError,
                                    'Sample data missing IDs '
                                    'in distance matrix'):
            _split_distance_matrix(self.dm, table)


if __name__ == '__main__':
    # run our unit tests!
    unittest.main()
