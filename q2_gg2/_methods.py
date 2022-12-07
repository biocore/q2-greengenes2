# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.tree import NewickFormat
from q2_types.feature_data import DNAFASTAFormat

import os
import json
import skbio
import biom
import pandas as pd
import bp
import re
from pathlib import Path
import hashlib

import redbiom.fetch as rbf
import redbiom.util as rbu


CONTEXTS = {
    90: 'Deblur_2021.09-Illumina-16S-V4-90nt-dd6875',
    100: 'Deblur_2021.09-Illumina-16S-V4-100nt-50b3a2',
    125: 'Deblur_2021.09-Illumina-16S-V4-125nt-92f954',
    150: 'Deblur_2021.09-Illumina-16S-V4-150nt-ac8c0b',
    200: 'Deblur_2021.09-Illumina-16S-V4-200nt-0b8b48',
    250: 'Deblur_2021.09-Illumina-16S-V4-250nt-8b2bff'
}


DF_COLUMNS = ['Feature ID', 'Taxon', 'Confidence']
PAD = [f'{r}__' for r in 'dpcofgs']


def _infer_feature_data_labels(df):
    gtdb = re.compile(r'G\d{9}')
    operon = re.compile(r'M.\d{3}-\d-barcode\d+-umi\d+bins-ubs-\d')
    asv = re.compile(r'\d{8}')
    silva = re.compile(r'[A-Z][A-Z0-9]?\d+\.\d+\.\d+')
    # LTP -> all else

    def labeler(id_):
        if gtdb.match(id_):
            return 'GTDB'
        elif operon.match(id_):
            return 'Operon'
        elif asv.match(id_):
            return 'ASV'
        elif silva.match(id_):
            return 'SILVA'
        else:
            return 'LTP'

    return pd.Series(df.index.map(labeler), index=df.index, name='Type')


class CladeAssessment(dict):
    SORTABLE = ('observed_v4_isolate', 'observed_v4_fragment',
                'in_clade_matching_v4_fragment_and_full_length',
                'out_clade_matching_v4_fragment_and_full_length',
                'unobserved_v4_isolate_fragments')
    UNSORTABLE = ('out_clade_lineages', 'redbiom')

    def __init__(self, *args, **kwargs):
        self['type'] = 'species-clade'
        super().__init__(*args, **kwargs)

    def __eq__(self, other):
        if len(self['clades']) != len(other['clades']):
            return False

        for s_clade, o_clade in zip(self['clades'], other['clades']):
            for k in self.SORTABLE:
                if sorted(s_clade[k]) != sorted(o_clade[k]):
                    return False

            for k in self.UNSORTABLE:
                if s_clade[k] != o_clade[k]:
                    return False

        return True


class ASVAssessment(dict):
    SORTABLE = ('observed_in_full_length', )
    UNSORTABLE = ('asv', 'id', 'lineage', 'md5', 'multifurcation_members',
                  'redbiom')

    def __init__(self, *args, **kwargs):
        self['type'] = 'asv-detail'
        super().__init__(*args, **kwargs)

    def __eq__(self, other):
        for k in self.SORTABLE:
            if sorted(self[k]) != sorted(other[k]):
                return False
        for k in self.UNSORTABLE:
            if self[k] != other[k]:
                return False
        return True


def clade_lookup(taxonomy_as_tree: NewickFormat,
                 version: str, output_filename: str) -> CladeAssessment:
    tree = str(taxonomy_as_tree)
    tree = bp.to_skbio_treenode(bp.parse_newick(open(str(tree)).read()))

    lookup_data = {}
    for n in tree.non_tips():
        name = n.name
        if n.name is None or len(n.name) == 3:  # eg s__
            continue

        if n.name not in lookup_data:
            lookup_data[n.name] = {'type': 'lookup', 'parent': set(),
                                   'children': [], 'name': n.name}

        parent = n.parent.name
        children = [c.name for c in n.children
                    if not c.is_tip() and len(c.name) > 3]

        # There is a thin edge case on tax2tree decoration where names
        # coming through secondary taxonomy *might* get multiple parents.
        # Due to a yet-to-be-determined caveat, secondary taxonomy decoration
        # is not constructing assured unique names when a name is polyphyletic.

        if parent is not None:
            lookup_data[n.name]['parent'].add(parent)

        lookup_data[n.name]['children'].extend(children)

    for k, v in lookup_data.items():
        v['parent'] = list(v['parent'])

    with open(output_filename, 'w') as fp:
        for k, v in sorted(lookup_data.items()):
            fp.write('\t'.join((k, version, json.dumps(v))))
            fp.write('\n')

    fake = {'observed_v4_isolate': '',
            'observed_v4_fragment': '',
            'in_clade_matching_v4_fragment_and_full_length': '',
            'out_clade_matching_v4_fragment_and_full_length': '',
            'unobserved_v4_isolate_fragments': '',
            'out_clade_lineages': '',
            'clades': []}
    return CladeAssessment(fake)


def _fast_parse_fasta(path):
    f = open(path)
    id_, seq = iter(f), iter(f)
    for i, s in zip(id_, seq):
        # drop > and \n
        yield i[1:-1], s[:-1]


def _stage_redbiom_metadata():
    restrict_to = ["empo_1", "empo_2", "empo_3", "sample_type", "latitude",
                   "longitude"]
    seen = set()
    md = []
    for v in CONTEXTS.values():
        # grab non-overlapping sets of samples as a sample may exist in
        # multiple contexts
        samples = rbf.samples_in_context(v, unambiguous=False)
        samples -= seen
        md.append(rbf.sample_metadata(samples, restrict_to=restrict_to)[0])
        seen.update(samples)

    md = pd.concat(md)
    md['latitude'] = pd.to_numeric(md['latitude'], errors='coerce')
    md['longitude'] = pd.to_numeric(md['longitude'], errors='coerce')

    return md


def _sequence_v4_redbiom_summary(asvs, rbmd):
    if rbmd.empty:
        return {}

    if len(asvs) == 0:
        return {}

    # partition asvs by their length
    by_length = {}
    for asv in asvs:
        asv_len = len(asv)
        if asv_len not in by_length:
            by_length[asv_len] = []
        by_length[asv_len].append(asv)

    # gather sample detail per context
    samples = set()
    for asv_len, asv_batch in by_length.items():
        ctx = CONTEXTS[asv_len]
        cur_samples = rbu.ids_from(asv_batch, exact=False, axis='feature',
                                   contexts=[ctx, ], min_count=1)
        _, _, _, sample_batch = rbu.partition_samples_by_tags(cur_samples)
        samples.update(sample_batch)

    # this is i think necessary as rbmd may lack a sample if it did not have
    # the desired sample information from staging
    samples = samples & set(rbmd.index)

    md_subset = rbmd.loc[samples]
    summary = {}
    for c in md_subset.columns:
        if c in ('latitude', 'longitude', '#SampleID'):
            continue
        summary[c] = md_subset[c].value_counts().to_dict()

    latlong = rbmd.loc[samples, ["latitude", "longitude"]].value_counts()
    latlongcount = []
    for (lat, long_), count in latlong.items():
        latlongcount.append([lat, long_, count])
    summary['latitude_longitude'] = latlongcount

    return {'context': CONTEXTS[min(by_length.keys())],
            'summary': summary}


def bulk_sequence_v4_asv_assessment(phylogeny: NewickFormat, taxa: pd.DataFrame,  # noqa
                                    full_length_v4: DNAFASTAFormat,
                                    sequences: DNAFASTAFormat,
                                    version: str,
                                    group: int, output_filename: str) -> ASVAssessment:  # noqa
    out = open(output_filename + '.%d' % group, 'w')

    # cache the things
    tree = bp.parse_newick(open(str(phylogeny)).read())
    full_length_v4 = {i: s for i, s in _fast_parse_fasta(str(full_length_v4))}
    sequences = {i: s for i, s in _fast_parse_fasta(str(sequences))
                 if i not in full_length_v4}
    taxa['Type'] = _infer_feature_data_labels(taxa)
    rbmd = _stage_redbiom_metadata()

    detail = {}
    for asv_id in taxa[taxa['Type'] == 'ASV'].index:
        asv_hash = hashlib.md5(asv_id.encode('ascii')).hexdigest()
        asv_group_int = int('0x%s' % asv_hash[:2], 0)
        if asv_group_int != group:
            continue

        asv = sequences[asv_id]
        start = time.time()
        detail = _sequence_v4_asv_assessment(tree, sequences, full_length_v4,
                                             taxa, asv, rbmd)

        asv_md5 = detail['md5']
        out.write('%s\t%s\t%s\n' % (asv_md5, version, json.dumps(detail)))
    out.close()
    return ASVAssessment(detail)


def sequence_v4_asv_assessment(phylogeny: NewickFormat, taxa: pd.DataFrame,
                               full_length_v4: DNAFASTAFormat,
                               sequences: DNAFASTAFormat,
                               asv: str) -> ASVAssessment:
    # cache the things
    tree = bp.parse_newick(open(str(phylogeny)).read())
    full_length_v4 = {i: s for i, s in _fast_parse_fasta(str(full_length_v4))}
    sequences = {i: s for i, s in _fast_parse_fasta(str(sequences))
                 if i not in full_length_v4}
    taxa['Type'] = _infer_feature_data_labels(taxa)
    start = time.time()
    rbmd = _stage_redbiom_metadata()
    return _sequence_v4_asv_assessment(tree, sequences, full_length_v4, taxa,
                                       asv, rbmd)


def _sequence_v4_asv_assessment(tree, sequences, full_length_v4, taxa, asv,
                                rbmd):
    # assumes tree is by ID
    # assumes taxonomy is by ID
    # assumes sequence data are by ID
    import time
    start = time.time()
    asv_id = None
    asv_hash = None
    for k, v in sequences.items():
        if v == asv:
            asv_id = k
            asv_hash = hashlib.md5(asv.encode('ascii')).hexdigest()
            break

    if asv_id is None:
        raise ValueError("ASV %s not found in sequences" % asv)

    start = time.time()
    node_idx = None
    for idx, v in enumerate(tree.B):
        if v:
            name = tree.name(idx)
            if name == asv_id:
                node_idx = idx
                break

    if node_idx is None:
        raise ValueError("ASV ID %s not found in phylogeny" % asv_id)

    lineage = taxa.loc[asv_id]['Taxon']

    # we expect the structure to be:
    # (full_length,(asv1,asv2,asv3,...))
    # so we need the parent of the current asv to obtain the multifurcation
    # and the parent of that to obtain full length fragments within the
    # clade

    # obtain other ASVs within the same multifurcation
    start = time.time()
    parent_start = tree.parent(node_idx)
    parent_close = tree.close(parent_start)
    multifurcation_members = []
    is_multifurcation = True
    for idx in range(parent_start, parent_close):
        # nsibling and psibling could be used but we'd have to iterate in
        # both directions so meh
        if tree.B[idx] and not tree.B[idx + 1]:
            name = tree.name(idx)
            if name != asv_id:
                if taxa.loc[name]['Type'] == 'ASV':
                    multifurcation_members.append(name)
                else:
                    is_multifurcation = False

    # determine what full length if any this ASV exists in
    start = time.time()
    observed_in = []
    for id_, seq in full_length_v4.items():
        if asv in seq:
            seq_lineage = taxa.loc[id_, 'Taxon']
            observed_in.append([id_, seq_lineage])


    start = time.time()
    rbdetail = _sequence_v4_redbiom_summary([asv, ], rbmd)

    return ASVAssessment({'asv': asv,
                          'id': asv_id,
                          'lineage': lineage,
                          'md5': asv_hash,
                          'observed_in_full_length': list(observed_in),
                          #'full_length_in_enclosing_clade': \
                          #    full_length_in_enclosing,
                          'multifurcation_members': len(multifurcation_members),
                          'redbiom': rbdetail})  # noqa


def bulk_clade_v4_asv_assessment(phylogeny: NewickFormat, taxa: pd.DataFrame,
                                 full_length_v4: DNAFASTAFormat,
                                 sequences: DNAFASTAFormat,
                                 version: str,
                                 group: int, output_filename: str) -> CladeAssessment:  # noqa
    out = open(output_filename + '.%d' % group, 'w')

    # cache the things
    tree = bp.parse_newick(open(str(phylogeny)).read())
    full_length_v4 = {i: s for i, s in _fast_parse_fasta(str(full_length_v4))}
    sequences = {i: s for i, s in _fast_parse_fasta(str(sequences))
                 if i not in full_length_v4}
    taxa['Type'] = _infer_feature_data_labels(taxa)

    rbmd = _stage_redbiom_metadata()

    detail = {}

    # for each species, if the species is in our processing group, summarize
    taxa['Species'] = taxa['Taxon'].apply(lambda x: x.split('; ')[-1])
    for species in taxa['Species'].unique():
        if species == 's__':
            continue
        species_group = hashlib.md5(species.encode('ascii')).hexdigest()

        # range is 0 -> 255
        species_group_int = int('0x%s' % species_group[:2], 0)
        if species_group_int != group:
            continue

        detail = _clade_v4_asv_assessment(tree, taxa, full_length_v4,
                                          sequences, species, rbmd)

        out.write('%s\t%s\t%s\n' % (species, version, json.dumps(detail)))

    out.close()

    return CladeAssessment(detail)


def clade_v4_asv_assessment(phylogeny: NewickFormat, taxa: pd.DataFrame,
                            full_length_v4: DNAFASTAFormat,
                            sequences: DNAFASTAFormat,
                            clade: str) -> CladeAssessment:
    tree = bp.parse_newick(open(str(phylogeny)).read())
    full_length_v4 = {i: s for i, s in _fast_parse_fasta(str(full_length_v4))}
    sequences = {i: s for i, s in _fast_parse_fasta(str(sequences))
                 if i not in full_length_v4}
    taxa['Type'] = _infer_feature_data_labels(taxa)

    rbmd = _stage_redbiom_metadata()
    return CladeAssessment(_clade_v4_asv_assessment(tree, taxa, full_length_v4,
                                                    sequences, clade, rbmd))


def _clade_v4_asv_assessment(tree, taxa, full_length_v4, sequences, clade,
                             rbmd):
    taxa['Species'] = taxa['Taxon'].apply(lambda x: x.split('; ')[-1])
    species_to_lineage = {r.Species: r.Taxon for r in taxa.itertuples()}

    # search for the clade in the tree
    clade_start = None
    details = []
    for idx, v in enumerate(tree.B):
        if v:
            name = tree.name(idx)

            # account for "g__foo; s__bar"
            if name and clade in name:
                clade_start = idx
                clade_end = tree.close(idx)
                detail = _clade_v4_assessment_tree_coordinates(tree, taxa,
                                                               full_length_v4,
                                                               sequences,
                                                               clade_start,
                                                               clade_end,
                                                               rbmd)
                detail['name'] = clade
                detail['lineage'] = species_to_lineage[clade]
                details.append(detail)

    if clade_start is None:
        return ValueError("%s not found in tree" % clade)

    return CladeAssessment({'clades': details})


def _clade_v4_assessment_tree_coordinates(tree, taxa, full_length_v4,
                                          sequences, clade_start, clade_end,
                                          rbmd):
    # obtain the set of identifiers at the tips of the clade
    clade_identifiers = {tree.name(idx)
                         for idx, v in enumerate(tree.B[clade_start:clade_end],
                                                 clade_start)
                         if v and not tree.B[idx + 1]}  # 10 indicates tip

    # filter the taxonomic information to any full_length in the clade
    clade_taxa = taxa.loc[clade_identifiers]
    isolates = clade_taxa[clade_taxa['Type'].isin(['GTDB', 'LTP'])]
    asvs = clade_taxa[clade_taxa['Type'].isin(['ASV'])]

    # filter to out-of-clade full length records
    outclade_taxa = taxa[(~taxa.index.isin(clade_identifiers)) &
                         (taxa['Type'].isin(['GTDB', 'LTP', 'Operon']))]

    if not len(isolates):
        # this should never happen as taxonomic information comes from the
        # isolates. Note it is permissible for there to not be any ASVs
        raise ValueError("No isolates observed within clade")

    clade_features = {full_length_v4[i] for i in clade_taxa.index
                      if i in full_length_v4}
    clade_features |= {sequences[i] for i in clade_taxa.index
                       if i in sequences}
    observed_v4_isolates = [[k, full_length_v4[k]] for k in isolates.index
                             if k in full_length_v4]
    observed_v4_fragments = [[k, sequences[k]] for k in asvs.index]
    observed_v4_fragments_lookup = {s for _, s in observed_v4_fragments}

    asv_rb_detail = _sequence_v4_redbiom_summary(
        list(observed_v4_fragments_lookup),
        rbmd
    )

    # test each fragment for exact match within clade isolates
    in_clade_v4_match = set()
    for isolate, isolate_frag in observed_v4_isolates:
        for asv, asv_frag in observed_v4_fragments:
            if asv_frag in isolate_frag:
                in_clade_v4_match.add((isolate, isolate_frag))
    in_clade_v4_match = [list(v) for v in in_clade_v4_match]

    # test each fragment for exact match outside of clade isolates
    out_clade_v4_match = set()
    for id_ in outclade_taxa.index:  # these are full length
        if id_ in full_length_v4:
            outclade_feature = full_length_v4[id_]
        elif id_ in sequences:
            outclade_feature = sequences[id_]
        else:
            # not all full length hybridize with EMP primers
            continue

        for inclade_feature in clade_features:  # these are any sequence
            if inclade_feature in outclade_feature:
                out_clade_v4_match.add((id_, outclade_feature))
    out_clade_v4_match = [list(v) for v in out_clade_v4_match]

    # test if any inclade isolate fragments are unobserved by ASV
    unobserved_v4_isolate_fragments = []
    for id_, seq in observed_v4_isolates:
        # our highest specificity fragments, which we also index in bulk
        # are 250nt. we could do substring testing against all outgroup
        # fragments but that may be annoyingly expensive

        # this isn't the greatest check, but checking even smaller substrings
        # is also not great. TODO: make better or remove
        if seq[:250] not in observed_v4_fragments_lookup:
            unobserved_v4_isolate_fragments.append(id_)

    # map fragments to associated IDs
    inverse_reads = {}
    for id_ in outclade_taxa.index:
        if id_ in full_length_v4:
            read = full_length_v4[id_]
        elif id_ in sequences:
            read = sequences[id_]
        else:
            # not all full length hybridize with EMP primers
            continue

        if read not in inverse_reads:
            inverse_reads[read] = set()
        lineage = outclade_taxa.loc[id_, 'Taxon']
        inverse_reads[read].add(lineage)

    # pull out the observed lineages for a given fragment
    out_clade_lineages = {}
    for _, fragment in out_clade_v4_match:
        lineages = list(inverse_reads.get(fragment, []))
        out_clade_lineages[fragment] = lineages

    return {'observed_v4_isolate': observed_v4_isolates,
            'observed_v4_fragment': observed_v4_fragments,
            'in_clade_matching_v4_fragment_and_full_length': in_clade_v4_match,
            'out_clade_matching_v4_fragment_and_full_length': out_clade_v4_match,
            'unobserved_v4_isolate_fragments': unobserved_v4_isolate_fragments,
            'out_clade_lineages': out_clade_lineages,
            'redbiom': asv_rb_detail}


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
