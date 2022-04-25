# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import io
import tempfile

import biom
import numpy as np
import pandas as pd
import pandas.testing as pdt

from q2_gg2._methods import (_load_tree_and_cache,
                             _fetch_taxonomy,
                             _fetch_unclassified,
                             _classify,
                             filter_features,
                             taxonomy_from_features,
                             taxonomy_from_table,
                             relabel)


class GG2MethodTests(unittest.TestCase):
    def setUp(self):
        self.treedata = io.StringIO("(((a,b)i1,(c,d)i2)i3,(e,f)i4)i5;")
        self.seqdata = (">a\nA\n"
                        ">b\nA\n"
                        ">c\nA\n"
                        ">d\nA\n"
                        ">e\nA\n"
                        ">f\nA\n"
                        ">X\nA\n"
                        ">Y\nA\n"
                        ">Z\nA\n")

    def test_load_tree_and_cache(self):
        obs = _load_tree_and_cache(self.treedata)
        for tip in obs.tips():
            ancestors = [n.name for n in tip.ancestors()]
            self.assertEqual(ancestors, tip.ancestor_cache)

    def test_fetch_taxonomy(self):
        tree = _load_tree_and_cache(self.treedata)
        exp = pd.DataFrame([['a', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['b', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['c', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['d', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['e', 'i5; i4; c__; o__; f__; g__; s__', 1.0],
                            ['f', 'i5; i4; c__; o__; f__; g__; s__', 1.0]],
                           columns=['Feature ID', 'Taxon', 'Confidence'])
        obs = _fetch_taxonomy(tree, 'abcdef')
        pdt.assert_frame_equal(obs, exp)

    def test_fetch_unclassified(self):
        exp = pd.DataFrame([['a', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['b', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['c', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['d', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['e', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['f', 'd__; p__; c__; o__; f__; g__; s__', 1.0]],
                           columns=['Feature ID', 'Taxon', 'Confidence'])
        obs = _fetch_unclassified('abcdef')
        pdt.assert_frame_equal(obs, exp)

    def test_classify(self):
        tree = _load_tree_and_cache(self.treedata)
        exp = pd.DataFrame([['a', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['b', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['c', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['d', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['e', 'i5; i4; c__; o__; f__; g__; s__', 1.0],
                            ['f', 'i5; i4; c__; o__; f__; g__; s__', 1.0],
                            ['X', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['Y', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['Z', 'd__; p__; c__; o__; f__; g__; s__', 1.0]],
                           columns=['Feature ID', 'Taxon', 'Confidence'])
        exp.set_index('Feature ID', inplace=True)
        obs = _classify(tree, set('abcdefXYZ'))
        obs = obs.loc[exp.index]
        pdt.assert_frame_equal(obs, exp)

    def test_taxonomy_from_features(self):
        exp = pd.DataFrame([['a', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['b', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['c', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['d', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['e', 'i5; i4; c__; o__; f__; g__; s__', 1.0],
                            ['f', 'i5; i4; c__; o__; f__; g__; s__', 1.0],
                            ['X', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['Y', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['Z', 'd__; p__; c__; o__; f__; g__; s__', 1.0]],
                           columns=['Feature ID', 'Taxon', 'Confidence'])
        exp.set_index('Feature ID', inplace=True)

        with tempfile.NamedTemporaryFile(delete=False) as tree, \
             tempfile.NamedTemporaryFile(delete=False) as reads:
            tree.write(self.treedata.read().encode('ascii'))
            tree.close()
            reads.write(self.seqdata.encode('ascii'))
            reads.close()
            obs = taxonomy_from_features(tree.name, reads.name)
        obs = obs.loc[exp.index]
        pdt.assert_frame_equal(obs, exp)

    def test_taxonomy_from_table(self):
        exp = pd.DataFrame([['a', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['b', 'i5; i3; i1; o__; f__; g__; s__', 1.0],
                            ['c', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['d', 'i5; i3; i2; o__; f__; g__; s__', 1.0],
                            ['e', 'i5; i4; c__; o__; f__; g__; s__', 1.0],
                            ['f', 'i5; i4; c__; o__; f__; g__; s__', 1.0],
                            ['X', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['Y', 'd__; p__; c__; o__; f__; g__; s__', 1.0],
                            ['Z', 'd__; p__; c__; o__; f__; g__; s__', 1.0]],
                           columns=['Feature ID', 'Taxon', 'Confidence'])
        exp.set_index('Feature ID', inplace=True)

        table = biom.Table(np.arange(27).reshape(9, 3),
                           list('abcdefXYZ'),
                           ['s1', 's2', 's3'])
        with tempfile.NamedTemporaryFile(delete=False) as tree:
            tree.write(self.treedata.read().encode('ascii'))
            tree.close()
            obs = taxonomy_from_table(tree.name, table)

        obs = obs.loc[exp.index]
        pdt.assert_frame_equal(obs, exp)

    def test_filter_features(self):
        table = biom.Table(np.arange(27).reshape(9, 3),
                           list('abcdefXYZ'),
                           ['s1', 's2', 's3'])
        tree = io.StringIO("(a,c,d,X);")
        exp = table.filter(set(['a', 'c', 'd', 'X']), axis='observation')
        obs = filter_features(table, tree)
        self.assertEqual(obs, exp)

    def test_relabel(self):
        table = biom.Table(np.arange(27).reshape(9, 3),
                           list('abcdefXYZ'),
                           ['s1', 's2', 's3'])
        mapping = pd.DataFrame([['a', 'a1', 'a2'],
                                ['b', 'b1', 'b2'],
                                ['c', 'c1', 'c2'],
                                ['d', 'd1', 'd2'],
                                ['e', 'e1', 'e2'],
                                ['f', 'f1', 'f2'],
                                ['X', 'X1', 'X2'],
                                ['Y', 'Y1', 'Y2'],
                                ['Z', 'Z1', 'Z2']],
                               columns=['id', 'md5', 'sequence'])

        with self.assertRaises(ValueError):
            relabel(table, mapping, as_md5=True, as_sequence=True)

        with self.assertRaises(ValueError):
            relabel(table, mapping, as_md5=True, as_sequence=True,
                    as_id=True)

        exp = ["%s1" % i for i in table.ids(axis='observation')]
        obs = relabel(table, mapping, as_md5=True)
        self.assertEqual(list(obs.ids(axis='observation')),
                         exp)

        exp = ["%s2" % i for i in table.ids(axis='observation')]
        obs = relabel(table, mapping, as_sequence=True)
        self.assertEqual(list(obs.ids(axis='observation')),
                         exp)


if __name__ == '__main__':
    unittest.main()
