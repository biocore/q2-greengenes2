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
import bp
import skbio

from q2_gg2._methods import (_load_tree_and_cache,
                             _fetch_taxonomy,
                             _fetch_unclassified,
                             _classify,
                             _clade_v4_asv_assessment,
                             _infer_feature_data_labels,
                             _sequence_v4_asv_assessment,
                             filter_features,
                             collapse_multifurcation,
                             taxonomy_from_features,
                             taxonomy_from_table,
                             relabel,
                             CladeAssessment,
                             ASVAssessment)


class GG2MethodTests(unittest.TestCase):
    def setUp(self):
        self.treedata = io.StringIO("((((a,b)i1,(c,d)i2)i3,(e,f)i4)i5);")
        self.tips = {'a', 'b', 'c', 'd', 'e', 'f'}
        self.seqdata = (">a\nA\n"
                        ">b\nA\n"
                        ">c\nA\n"
                        ">d\nA\n"
                        ">e\nA\n"
                        ">f\nA\n"
                        ">X\nA\n"
                        ">Y\nA\n"
                        ">Z\nA\n")

    def test_collapse_multifurcation(self):
        table = biom.Table(np.arange(24).reshape(6, 4),
                           ['10000000', '20000000', '30000000',
                            '40000000', '50000000', '60000000'],
                           list('abcd'))
        tree = io.StringIO("((((10000000,20000000)x,G1),(30000000,G2)y),((40000000,50000000,60000000)z,(G3,G4)));")  # noqa
        exp_tab = biom.Table(np.array([[4, 6, 8, 10],
                                       [8, 9, 10, 11],
                                       [48, 51, 54, 57]]),
                             ['x', 30000000, 'z'],
                             list('abcd'))
        exp_tree = skbio.TreeNode.read(["(((x,G1),(30000000,G2)y),(z,(G3,G4)));"])  # noqa
        obs_tab, obs_tree = collapse_multifurcation(table, tree)
        self.assertEqual(obs_tab, exp_tab)
        self.assertEqual(obs_tree.compare_rfd(exp_tree), 0.)

    def test_infer_feature_data_labels(self):
        taxa = pd.DataFrame([['MJ007-1-barcode27-umi40bins-ubs-7010', 's__foo',
                              0.1],
                             ['X80725', 's__foo', 0.1],
                             ['00000000', 's__foo', 0.1],
                             ['G012345678', 's__bar', 0.1],
                             ['11111111', 's__bar', 0.1]],
                            columns=['feature-id', 'Taxon', 'Confidence'])
        taxa.set_index('feature-id', inplace=True)
        exp = pd.DataFrame([['MJ007-1-barcode27-umi40bins-ubs-7010', 's__foo',
                             0.1, 'Operon'],
                            ['X80725', 's__foo', 0.1, 'LTP'],
                            ['00000000', 's__foo', 0.1, 'ASV'],
                            ['G012345678', 's__bar', 0.1, 'GTDB'],
                            ['11111111', 's__bar', 0.1, 'ASV']],
                           columns=['feature-id', 'Taxon', 'Confidence',
                                    'Type']).set_index('feature-id')
        obs = _infer_feature_data_labels(taxa)
        pdt.assert_series_equal(obs, exp['Type'])

    def test_clade_v4_asv_assessment(self):
        tree = bp.parse_newick("(((a,00000000,00000001,b)s__foo),(c,11111111)s__bar);")  # noqa
        taxa = pd.DataFrame([['a', 's__foo', 0.1],
                             ['b', 's__foo', 0.1],
                             ['00000000', 's__foo', 0.1],
                             ['00000001', 's__foo', 0.1],
                             ['c', 's__bar', 0.1],
                             ['11111111', 's__bar', 0.1]],
                            columns=['feature-id', 'Taxon', 'Confidence'])
        taxa.set_index('feature-id', inplace=True)
        taxa['Type'] = _infer_feature_data_labels(taxa)

        # note: this represents EXTRACTED isolate v4
        sequences = (">a\n"
                     "AACATGAA\n"
                     ">00000001\n"
                     "CATT\n"
                     ">b\n"
                     "CCCATGCC\n"
                     ">00000000\n"
                     "CATG\n"
                     ">c\n"
                     "AGGGCCAG\n"
                     ">11111111\n"
                     "GGCC\n")
        isolates_v4 = (">a\n"
                       "CATG\n"
                       ">b\n"
                       "CATG\n"
                       ">c\n"
                       "GGCC\n")
        seqs = {r.metadata['id']: str(r)
                for r in skbio.read(io.StringIO(sequences),
                                    format='fasta',
                                    constructor=skbio.DNA)}
        iso = {r.metadata['id']: str(r)
               for r in skbio.read(io.StringIO(isolates_v4),
                                   format='fasta',
                                   constructor=skbio.DNA)}
        obs = _clade_v4_asv_assessment(tree, taxa, iso, seqs, 's__foo',
                                       pd.DataFrame())
        exp = {'observed_v4_isolate': [['a', 'CATG'], ['b', 'CATG']],
               'observed_v4_fragment': [['00000000', 'CATG'],
                                        ['00000001', 'CATT']],
               'in_clade_matching_v4_fragment_and_full_length': [['a', 'CATG'],
                                                                 ['b', 'CATG']],  # noqa
               'out_clade_matching_v4_fragment_and_full_length': [],
               'out_clade_lineages': {},
               'unobserved_v4_isolate_fragments': [],
               'name': 's__foo',
               'lineage': 's__foo',
               'redbiom': {}
               }
        self.assertEqual(CladeAssessment(obs),
                         CladeAssessment({'clades': [exp, ]}))

        # example with v4 being insufficient to differentiate clades
        sequences = (">a\n"
                     "ATCATG\n"
                     ">b\n"
                     "CATGGG\n"
                     ">00000001\n"
                     "CAT\n"
                     ">00000000\n"
                     "ATG\n"
                     ">c\n"
                     "CATGA\n"
                     ">11111111\n"
                     "ATG\n")
        isolates_v4 = (">a\n"
                       "CATG\n"
                       ">b\n"
                       "CATG\n"
                       ">c\n"
                       "CATG\n")
        seqs = {r.metadata['id']: str(r)
                for r in skbio.read(io.StringIO(sequences),
                                    format='fasta',
                                    constructor=skbio.DNA)}
        iso = {r.metadata['id']: str(r)
               for r in skbio.read(io.StringIO(isolates_v4),
                                   format='fasta',
                                   constructor=skbio.DNA)}
        obs = _clade_v4_asv_assessment(tree, taxa, iso, seqs, 's__foo',
                                       pd.DataFrame())
        exp = {'observed_v4_isolate': [['a', 'CATG'], ['b', 'CATG']],
               'observed_v4_fragment': [['00000000', 'ATG'],
                                        ['00000001', 'CAT']],
               'in_clade_matching_v4_fragment_and_full_length': [['a', 'CATG'],
                                                                 ['b', 'CATG']],  # noqa
               'out_clade_matching_v4_fragment_and_full_length': [['c', 'CATG'], ],  # noqa
               'out_clade_lineages': {'CATG': ['s__bar', ]},
               'unobserved_v4_isolate_fragments': ['a', 'b'],
               'name': 's__foo',
               'lineage': 's__foo',
               'redbiom': {}
               }
        self.assertEqual(CladeAssessment(obs),
                         CladeAssessment({'clades': [exp, ]}))

    def test_sequence_v4_asv_assessment(self):
        tree = bp.parse_newick("((((a,(00000000,00000001)),b)s__foo),(c,(11111111))s__bar);")  # noqa
        taxa = pd.DataFrame([['a', 's__foo', 0.1],
                             ['b', 's__foo', 0.1],
                             ['00000000', 's__foo', 0.1],
                             ['00000001', 's__foo', 0.1],
                             ['c', 's__bar', 0.1],
                             ['11111111', 's__bar', 0.1]],
                            columns=['feature-id', 'Taxon', 'Confidence'])
        taxa.set_index('feature-id', inplace=True)
        taxa['Type'] = _infer_feature_data_labels(taxa)
        asv = 'CATG'

        # note: this represents EXTRACTED isolate v4
        sequences = (">a\n"
                     "AACATGAA\n"
                     ">00000001\n"
                     "CATT\n"
                     ">b\n"
                     "CCCATGCC\n"
                     ">00000000\n"
                     "CATG\n"
                     ">c\n"
                     "AGGGCCAG\n"
                     ">11111111\n"
                     "GGCC\n")
        isolates_v4 = (">a\n"
                       "CATG\n"
                       ">b\n"
                       "CATG\n"
                       ">c\n"
                       "GGCC\n")

        seqs = {r.metadata['id']: str(r)
                for r in skbio.read(io.StringIO(sequences),
                                    format='fasta',
                                    constructor=skbio.DNA)}
        iso = {r.metadata['id']: str(r)
               for r in skbio.read(io.StringIO(isolates_v4),
                                   format='fasta',
                                   constructor=skbio.DNA)}
        exp = {'asv': 'CATG',
               'id': '00000000',
               'lineage': 's__foo',
               'md5': '2619b1f331be094f6beb73877d078a7a',
               'observed_in_full_length': [['a', 's__foo'], ['b', 's__foo']],
               'full_length_in_enclosing_clade': ['a', ],
               'multifurcation_members': 1,  # ['00000001', ]}
               'redbiom': {}}
        obs = _sequence_v4_asv_assessment(tree, seqs, iso, taxa, asv,
                                          pd.DataFrame())
        self.assertEqual(obs, ASVAssessment(exp))

        tree = bp.parse_newick("((((a,00000000),b)s__foo),(c,(11111111))s__bar);");  # noqa
        exp = {'asv': 'CATG',
               'id': '00000000',
               'lineage': 's__foo',
               'md5': '2619b1f331be094f6beb73877d078a7a',
               'observed_in_full_length': [['a', 's__foo'], ['b', 's__foo']],
               'full_length_in_enclosing_clade': ['a', ],
               'multifurcation_members': 0,
               'redbiom': {}}  # []}
        obs = _sequence_v4_asv_assessment(tree, seqs, iso, taxa, asv,
                                          pd.DataFrame())
        self.assertEqual(obs, ASVAssessment(exp))

    def test_load_tree_and_cache(self):
        obs = _load_tree_and_cache(self.treedata, self.tips)
        for tip in obs.tips():
            ancestors = [n.name for n in tip.ancestors()][:-1]
            self.assertEqual(ancestors, tip.ancestor_cache)

    def test_fetch_taxonomy(self):
        tree = _load_tree_and_cache(self.treedata, self.tips)
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
        tree = _load_tree_and_cache(self.treedata, self.tips)
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
        mapping = pd.DataFrame([['a', 'a1', 'a2' * 500],
                                ['b', 'b1', 'b2'],
                                ['c', 'c1', 'c2'],
                                ['d', 'd1', 'd2'],
                                ['e', 'e1', 'e2'],
                                ['f', 'f1', 'f2'],
                                ['X', 'X1', 'X2'],
                                ['Y', 'Y1', 'Y2'],
                                ['Z', 'Z1', 'Z2' * 500]],
                               columns=['id', 'md5', 'sequence'])

        with self.assertRaises(ValueError):
            relabel(table, mapping, as_md5=True, as_asv=True)

        with self.assertRaises(ValueError):
            relabel(table, mapping, as_md5=True, as_asv=True,
                    as_id=True)

        exp = ["%s1" % i for i in table.ids(axis='observation')]
        obs = relabel(table, mapping, as_md5=True)
        self.assertEqual(list(obs.ids(axis='observation')),
                         exp)

        # for ASV, only use ASV if the sequence appears short
        exp = ["%s2" % i for i in table.ids(axis='observation')]
        exp[0] = 'a'
        exp[-1] = 'Z'
        obs = relabel(table, mapping, as_asv=True)
        self.assertEqual(list(obs.ids(axis='observation')),
                         exp)


if __name__ == '__main__':
    unittest.main()
