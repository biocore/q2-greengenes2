# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.tree import NewickFormat
from q2_types.feature_data import DNAFASTAFormat

import skbio
import biom
import pandas as pd
import bp


DF_COLUMNS = ['Feature ID', 'Taxon', 'Confidence']
PAD = [f'{r}__' for r in 'dpcofgs']


def _load_tree_and_cache(data):
    treedata = data.read()
    tree = bp.to_skbio_treenode(bp.parse_newick(treedata))
    tree.ancestor_cache = []
    for node in tree.preorder(include_self=False):
        node.ancestor_cache = [node.parent.name] + node.parent.ancestor_cache
    return tree


def _fetch_taxonomy(tree, features):
    results = []
    for feature in features:
        lineage = PAD[:]
        tip = tree.find(feature)

        for i, name in enumerate(tip.ancestor_cache[::-1]):
            lineage[i] = name
        results.append((feature, '; '.join(lineage), 1.0))

    classification = pd.DataFrame(results, columns=DF_COLUMNS)
    return classification


def _fetch_unclassified(features):
    flat = '; '.join(PAD)
    results = [(feature, flat, 1.0) for feature in features]
    return pd.DataFrame(results, columns=DF_COLUMNS)


def _classify(tree, features):
    names = {n.name for n in tree.tips()}
    classifiable = names & features
    unclassified = features - names

    df = pd.concat([_fetch_taxonomy(tree, classifiable),
                    _fetch_unclassified(unclassified)])

    return df.set_index('Feature ID')


def taxonomy_from_features(reference_taxonomy: NewickFormat,
                           reads: DNAFASTAFormat) -> pd.DataFrame:
    tree = _load_tree_and_cache(open(str(reference_taxonomy)))
    features = {r.metadata['id'] for r in skbio.read(str(reads),
                                                     format='fasta',
                                                     constructor=skbio.DNA)}
    return _classify(tree, features)


def taxonomy_from_table(reference_taxonomy: NewickFormat,
                        table: biom.Table) -> pd.DataFrame:
    tree = _load_tree_and_cache(open(str(reference_taxonomy)))
    features = set(table.ids(axis='observation'))

    return _classify(tree, features)


def filter_features(table: biom.Table, reference: NewickFormat) -> biom.Table:
    treedata = reference.read()
    tree = bp.parse_newick(treedata)

    names = {tree.name(i) for i, v in enumerate(tree.B) if v}

    return table.filter(set(table.ids(axis='observation')) & names,
                        axis='observation',
                        inplace=False).remove_empty()


def relabel(feature_table: biom.Table,
            reference_label_map: pd.DataFrame,
            as_md5: bool = False,
            as_sequence: bool = False,
            as_id: bool = False) -> biom.Table:
    if int(as_md5) + int(as_sequence) + int(as_id) > 1:
        raise ValueError("Only a single conversion type can be specified")

    if as_md5:
        key = 'md5'
    elif as_sequence:
        key = 'sequence'
    elif as_id:
        key = 'id'

    ids = set(feature_table.ids(axis='observation'))

    currently_as_ids = len(ids & set(reference_label_map['id']))
    currently_as_sequence = len(ids & set(reference_label_map['sequence']))
    currently_as_md5 = len(ids & set(reference_label_map['md5']))

    if currently_as_ids > max(currently_as_sequence, currently_as_md5):
        current = 'id'
    elif currently_as_md5 > max(currently_as_sequence, currently_as_ids):
        current = 'md5'
    else:
        current = 'sequence'

    src_dst = reference_label_map[[current, key]]
    src_dst_map = {getattr(r, current): getattr(r, key)
                   for r in src_dst.itertuples()}

    return feature_table.update_ids(src_dst_map, inplace=False,
                                    axis='observation')
