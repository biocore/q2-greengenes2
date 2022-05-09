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


def _load_tree_and_cache(data, features):
    treedata = data.read()
    tree = bp.parse_newick(treedata)
    names = {tree.name(i) for i, v in enumerate(tree.B) if v}
    tree = tree.shear(names & features)
    tree = bp.to_skbio_treenode(tree)

    tree.ancestor_cache = []
    for node in tree.preorder(include_self=False):
        if node.parent.is_root():
            node.ancestor_cache = []
        else:
            parent_cache = node.parent.ancestor_cache
            node.ancestor_cache = [node.parent.name] + parent_cache
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
    features = {r.metadata['id'] for r in skbio.read(str(reads),
                                                     format='fasta',
                                                     constructor=skbio.DNA)}
    tree = _load_tree_and_cache(open(str(reference_taxonomy)), features)
    return _classify(tree, features)


def taxonomy_from_table(reference_taxonomy: NewickFormat,
                        table: biom.Table) -> pd.DataFrame:
    features = set(table.ids(axis='observation'))
    tree = _load_tree_and_cache(open(str(reference_taxonomy)), features)

    return _classify(tree, features)


def filter_features(feature_table: biom.Table,
                    reference: NewickFormat) -> biom.Table:
    try:
        treedata = reference.read()
    except:
        treedata = open(str(reference)).read()
    tree = bp.parse_newick(treedata)

    names = {tree.name(i) for i, v in enumerate(tree.B) if v}
    overlap = set(feature_table.ids(axis='observation')) & names
    return feature_table.filter(overlap,
                                axis='observation',
                                inplace=False).remove_empty()


def relabel(feature_table: biom.Table,
            reference_label_map: pd.DataFrame,
            as_md5: bool = False,
            as_asv: bool = False,
            as_id: bool = False) -> biom.Table:
    if int(as_md5) + int(as_asv) + int(as_id) > 1:
        raise ValueError("Only a single conversion type can be specified")

    if as_md5:
        key = 'md5'
    elif as_asv:
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

    if as_asv:
        # consider ASVs sub 500nt
        is_asv = reference_label_map['sequence'].apply(lambda x: len(x) < 500)
        subset = reference_label_map[is_asv][[current, 'sequence']]
        src_dst_map = subset.set_index(current)['sequence'].to_dict()

        if current == 'id':
            ids = reference_label_map[~is_asv]['id']
            to_id = {i: i for i in ids}
        else:
            outofsubset = reference_label_map[~is_asv][[current, 'id']]
            to_id = outofsubset.set_index(current)['id'].to_dict()
        src_dst_map.update(to_id)
    else:
        src_dst = reference_label_map[[current, key]]
        src_dst_map = {getattr(r, current): getattr(r, key)
                       for r in src_dst.itertuples()}

    return feature_table.update_ids(src_dst_map, inplace=False,
                                    axis='observation')
