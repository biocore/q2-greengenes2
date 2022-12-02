import qiime2
from qiime2.plugins import evident
from skbio import DistanceMatrix
from q2_gg2._methods import (_split_distance_matrix,
                             _compute_effect_size)
from pandas.testing import assert_frame_equal
import pandas as pd
import numpy as np
import unittest


class UtilTests(unittest.TestCase):
    def setUp(self):

        # construct distance matrix
        mat = np.random.uniform(0, 1, size=(100, 100))
        mat += mat.T
        np.fill_diagonal(mat, 0)
        dm = DistanceMatrix(mat,
                            ["A%d" % i for i in range(1, len(mat) + 1)])

        # import metadata
        table = pd.read_csv('samp.tsv', sep='\t')
        self.table = table.set_index('sampleid')

        # split distance
        self.split = _split_distance_matrix(dm, self.table)

        # convert objects into Artifacts
        self.table = qiime2.metadata.Metadata.load('samp.tsv')
        self.table_16s = self.table.filter_ids(self.split[0].ids)
        self.table_wgs = self.table.filter_ids(self.split[1].ids)

    def test_compute_effect_size_no_columns(self):
        _16s = qiime2.Artifact.import_data('DistanceMatrix', self.split[0])
        wgs = qiime2.Artifact.import_data('DistanceMatrix', self.split[1])
        self.table.filter_columns(column_type='categorical')
        columns = list(self.table.to_dataframe().columns.values)
        columns.remove('paired_sample')
        columns.remove('preparation')
        self.table_16s = self.table_16s.to_dataframe()[columns]
        self.table_wgs = self.table_wgs.to_dataframe()[columns]
        self.table_16s = qiime2.Metadata(self.table_16s)
        self.table_wgs = qiime2.Metadata(self.table_wgs)
        _16s_eff_size, = evident.methods.multivariate_effect_size_by_category(
                            data=_16s,
                            sample_metadata=self.table_16s,
                            group_columns=columns)
        _wgs_eff_size, = evident.methods.multivariate_effect_size_by_category(
                            data=wgs,
                            sample_metadata=self.table_wgs,
                            group_columns=columns)
        effect_size_16 = _16s_eff_size.view(pd.DataFrame)
        effect_size_wgs = _wgs_eff_size.view(pd.DataFrame)
        effect_size_16.rename(columns={'effect_size': 'effect_size_16s',
                                       'metric': 'metric_16s'},
                              inplace=True)
        effect_size_wgs.rename(columns={'effect_size': 'effect_size_wgs',
                                        'metric': 'metric_wgs'},
                               inplace=True)
        effect_sizes = pd.merge(effect_size_16, effect_size_wgs,
                                how='outer', on='column')
        assert_frame_equal(_compute_effect_size(self.split[0], self.split[1],
                                                self.table), effect_sizes)

    def test_compute_effect_size_with_columns(self):
        _16s = qiime2.Artifact.import_data('DistanceMatrix', self.split[0])
        wgs = qiime2.Artifact.import_data('DistanceMatrix', self.split[1])
        columns = ['category1', 'category2']
        self.table_16s = self.table_16s.to_dataframe()[columns]
        self.table_wgs = self.table_wgs.to_dataframe()[columns]
        self.table_16s = qiime2.Metadata(self.table_16s)
        self.table_wgs = qiime2.Metadata(self.table_wgs)
        _16s_eff_size, = evident.methods.multivariate_effect_size_by_category(
                            data=_16s,
                            sample_metadata=self.table_16s,
                            group_columns=columns)
        _wgs_eff_size, = evident.methods.multivariate_effect_size_by_category(
                            data=wgs,
                            sample_metadata=self.table_wgs,
                            group_columns=columns)
        effect_size_16 = _16s_eff_size.view(pd.DataFrame)
        effect_size_wgs = _wgs_eff_size.view(pd.DataFrame)
        effect_size_16.rename(columns={'effect_size': 'effect_size_16s',
                                       'metric': 'metric_16s'},
                              inplace=True)
        effect_size_wgs.rename(columns={'effect_size': 'effect_size_wgs',
                                        'metric': 'metric_wgs'},
                               inplace=True)
        effect_sizes = pd.merge(effect_size_16, effect_size_wgs,
                                how='outer', on='column')
        assert_frame_equal(_compute_effect_size(self.split[0], self.split[1],
                                                self.table,
                                                ['category1', 'category2']),
                           effect_sizes)

    def test_compute_effect_size_raises_columns_missing(self):
        with self.assertRaisesRegex(KeyError,
                                    'Must include at least one '
                                    'categorical column'):
            _compute_effect_size(self.split[0], self.split[1],
                                 self.table, columns=[])

    def test_compute_effect_size_raises_matrix_dimension(self):
        _16s_data = np.array([[0]])
        id = ['A1']
        dm = DistanceMatrix(_16s_data, id)
        with self.assertRaisesRegex(ValueError,
                                    'Distance matrix dimensions '
                                    'must be larger than 1x1'):
            _compute_effect_size(dm, self.split[1], self.table)

    def test_compute_effect_size_raises_columns_invalid(self):
        with self.assertRaisesRegex(KeyError,
                                    'Columns are not defined in metadata'):
            _compute_effect_size(self.split[0], self.split[1],
                                 self.table, columns=['a', 'b'])


if __name__ == '__main__':
    # run our unit tests
    unittest.main()
