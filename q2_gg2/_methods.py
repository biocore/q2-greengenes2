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
import q2templates
import os
import pkg_resources
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from qiime2.plugins import evident
import qiime2


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


def _split_distance_matrix(distance_matrix,
                           metadata, strict=False):
    """Split a distance matrix into its paired components

    Parameters
    ----------
    distance_matrix : skbio.DistanceMatrix
        The distance matrix to split into paired components
    metadata : pd.DataFrame
        The metadata to use to inform the split. We are assuming the DataFrame
        contains a column called 'paired_sample', which describes pairing
        relationships, and a column called 'preparation' which describes
        whether the sample is '16S' or 'WGS'.
    strict : bool
        If true, limit the paired components to only samples which exist in
        both preparations, e.g., a sample's 'paired_sample' also exists in the
        distance matrix

    Notes
    -----
    The sample metadata may be a superset of the samples in the distance matrix

    Raises
    ------
    KeyError
        If 'paired_sample' or 'preparation' do not exist as columns in the
        DataFrame
    ValueError
        If more than 2 or less than 2 unique 'preparation's, with samples in
        the distance matrix, exists

    Returns
    -------
    skbio.DistanceMatrix, skbio.DistanceMatrix
        A tuple where the first distance matrix are the 16S data, and the
        second distance matrix are the WGS data.
    """
    # RAISE ERRORS
    if 'paired_sample' not in metadata.columns:
        raise KeyError('Metadata must include paired_sample')
    elif 'preparation' not in metadata.columns:
        raise KeyError('Metadata must include preparation')
    elif any(metadata['preparation'].isnull()):
        raise ValueError('Missing preparations')
    elif len(set(metadata['preparation'])) != 2:
        raise ValueError('Must have exactly 2 unique prepartions')
    elif any(i not in list(metadata.index) for i in distance_matrix.ids):
        raise ValueError('Sample data missing IDs in distance matrix')

    # Separate 16s and WGS
    _16s = metadata[metadata['preparation'] == '16S']
    _sg = metadata[metadata['preparation'] == 'WGS']

    # filter ID's that are in distance matrix
    _16s = _16s[_16s.index.isin(distance_matrix.ids)]
    _sg = _sg[_sg.index.isin(distance_matrix.ids)]

    # check strict
    if strict is True:
        distance_matrix_16s = distance_matrix.filter(
                              _16s[_16s.index.isin(
                                _sg['paired_sample'])].index)
        distance_matrix_sg = distance_matrix.filter(
                             _sg[_sg.index.isin(
                                _16s['paired_sample'])].index)
    else:
        distance_matrix_16s = distance_matrix.filter(_16s.index)
        distance_matrix_sg = distance_matrix.filter(_sg.index)

# return tuple
    return (distance_matrix_16s, distance_matrix_sg)


def _compute_effect_size(distance_matrix_16s, distance_matrix_wgs,
                         metadata, columns=None, _max_level_by_category=5):

    """ 
    Computes effect sizes of 16S and WGS distance
    matrices from categorical columns in metadata

    Parameters
    ----------
    distance_matrix_16s : skbio.DistanceMatrix
        Distance matrix of 16S genes
    distance_matrix_wgs : skbio.DistanceMatrix
        Distance matrix of WGS genes
    metadata : Metadata
        The metadata with categorical columns to compute
        the effect sizes.
    columns : list of strings
        If None, all categorical columns will be used to
        compute effect sizes. Else, only specified columns
        will be used to compute effect sizes.
    _max_level_by_category : int
        The maximum number of levels allowed per categorical
        column. Default is 5.

    Raises
    ------
    KeyError
        If there is less than one categorical column, or if
        specified columns are not present in the metadata.
    ValueError
        If either distance matrix dimensions are less than 2

    Returns
    -------
    pd.Dataframe
        Dataframe that includes the effect size and metric
        of both the 16s and WGS data along with the
        column that the effect size is computed from.
    """

    # set optional columns argument
    if columns is None:
        metadata.filter_columns(column_type='categorical')
        columns = list(metadata.to_dataframe().columns.values)
        columns.remove('paired_sample')
        columns.remove('preparation')

    # RAISE ERRORS IF
    # check if there exists at least one categorical column to
    # compute effect size
    if len(columns) < 1:
        raise KeyError('Must include at least one categorical column')
    # check if the split distance matrix dimensions are greater than 1
    elif distance_matrix_16s.shape[0] < 2 or distance_matrix_wgs.shape[0] < 2:
        raise ValueError('Distance matrix dimensions must be larger than 1x1')
    # check if columns argument are columns in the metadata
    elif any(i not in list(metadata.to_dataframe().columns.values)
             for i in columns):
        raise KeyError('Columns are not defined in metadata')

    # separate Metadata between 16s and WGS
    metadata_16s = metadata.filter_ids(distance_matrix_16s.ids)
    metadata_wgs = metadata.filter_ids(distance_matrix_wgs.ids)

    # only include relevant columns in metadata
    # before passing it into effect size calculations
    metadata_16s = metadata_16s.to_dataframe()[columns]
    metadata_wgs = metadata_wgs.to_dataframe()[columns]
    metadata_16s = qiime2.Metadata(metadata_16s)
    metadata_wgs = qiime2.Metadata(metadata_wgs)

    # convert skbio.DistanceMatrix to DistanceMatrix Artifact
    _16s = qiime2.Artifact.import_data('DistanceMatrix', distance_matrix_16s)
    wgs = qiime2.Artifact.import_data('DistanceMatrix', distance_matrix_wgs)

    # calculate effect sizes
    _16s_effect_size, = evident.methods.multivariate_effect_size_by_category(
                            data=_16s,
                            sample_metadata=metadata_16s,
                            group_columns=columns,
                            max_levels_per_category=_max_level_by_category)
    _wgs_effect_size, = evident.methods.multivariate_effect_size_by_category(
                            data=wgs,
                            sample_metadata=metadata_wgs,
                            group_columns=columns,
                            max_levels_per_category=_max_level_by_category)

    # convert effect size to pandas dataframe
    effect_size_16 = _16s_effect_size.view(pd.DataFrame)
    effect_size_wgs = _wgs_effect_size.view(pd.DataFrame)
    effect_size_16.rename(columns={'effect_size': 'effect_size_16s',
                                   'metric': 'metric_16s'}, inplace=True)
    effect_size_wgs.rename(columns={'effect_size': 'effect_size_wgs',
                                    'metric': 'metric_wgs'}, inplace=True)
    effect_sizes = pd.merge(effect_size_16, effect_size_wgs,
                            how='outer', on='column')
    return(effect_sizes)


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
    except AttributeError:
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


def compute_effect_size(output_dir: str,
                        distance_matrix: skbio.DistanceMatrix,
                        metadata: qiime2.Metadata,
                        strict: bool = False,
                        columns: list = None,
                        _max_level_by_category: int = 5) -> None:
    distance_matrix_16s, distance_matrix_wgs = _split_distance_matrix(
                                                     distance_matrix,
                                                     metadata.to_dataframe(),
                                                     strict)
    effect_sizes = _compute_effect_size(
                        distance_matrix_16s,
                        distance_matrix_wgs,
                        metadata, columns,
                        _max_level_by_category)

    # compute correlation coefficient and p-value
    if (len(effect_sizes['column']) > 1):
        pearson_result = stats.pearsonr(effect_sizes['effect_size_16s'],
                                        effect_sizes['effect_size_wgs'])
    else:
        pearson_result = (0, 0)

    # create scatterplot
    m, b = np.polyfit(effect_sizes['effect_size_16s'],
                      effect_sizes['effect_size_wgs'], 1)
    plt.scatter(effect_sizes['effect_size_16s'],
                effect_sizes['effect_size_wgs'])
    plt.plot(effect_sizes['effect_size_16s'], m *
             effect_sizes['effect_size_16s'] + b)
    plt.xlim([min(effect_sizes['effect_size_16s']) * .9,
              max(effect_sizes['effect_size_16s']) * 1.1])
    plt.ylim([min(effect_sizes['effect_size_16s']) * .9,
              max(effect_sizes['effect_size_wgs']) * 1.1])
    plt.title('Scatterplot of Effect Sizes')
    plt.xlabel('16s Effect Sizes')
    plt.ylabel('WGS Effect Sizes')
    plt.text(0.1, 0.02, 'r = %0.2f, p = %0.2e' % pearson_result,
             transform=plt.gcf().transFigure, color='red')
    plt.savefig(os.path.join(output_dir, 'scatter.png'))
    plt.close()

    # download effect size table
    effect_sizes.to_csv(os.path.join(output_dir, 'table.tsv'), sep='\t')

    TEMPLATES = pkg_resources.resource_filename(
        'q2_gg2', 'compute_effect_size_assets')
    index = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={})
