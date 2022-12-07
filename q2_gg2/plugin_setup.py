# ----------------------------------------------------------------------------
# Copyright (c) 2022-, Greengenes2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib
from qiime2.plugin import Plugin, Bool, Str, Int
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Taxonomy, Sequence
from q2_types.tree import Phylogeny, Rooted
import q2_gg2
from ._type import ReferenceMap, CladeAssessment, ASVAssessment
from ._format import (ReferenceMapDirFmt, CladeAssessmentDirFmt,
                      ASVAssessmentDirFmt)


plugin = Plugin(
    name='greengenes2',
    version=q2_gg2.__version__,
    website='https://github.com/biocore/q2-gg2',
    package='q2_gg2',
    description="Methods for interaction with Greengenes2",
    short_description="Methods for interaction with Greengenes2"
)


plugin.register_semantic_types(
    ReferenceMap, CladeAssessment, ASVAssessment
)


plugin.register_formats(ReferenceMapDirFmt,
                        CladeAssessmentDirFmt,
                        ASVAssessmentDirFmt)


plugin.register_semantic_type_to_format(
    ReferenceMap,
    artifact_format=ReferenceMapDirFmt
)


plugin.register_semantic_type_to_format(
    CladeAssessment,
    artifact_format=CladeAssessmentDirFmt
)


plugin.register_semantic_type_to_format(
    ASVAssessment,
    artifact_format=ASVAssessmentDirFmt
)


plugin.methods.register_function(
    function=q2_gg2.clade_v4_asv_assessment,
    inputs={'phylogeny': Phylogeny[Rooted],
            'taxa': FeatureData[Taxonomy],
            'full_length_v4': FeatureData[Sequence],
            'sequences': FeatureData[Sequence]},
    parameters={'clade': Str},
    outputs=[('characterization', CladeAssessment)],
    input_descriptions={
        'phylogeny': 'The ID annotated phylogeny',
        'taxa': 'The taxa indexed by ID',
        'full_length_v4': 'The extracted V4 sequences from the backbone',
        'sequences': 'The full set of GG2 sequences'},
    parameter_descriptions={
        'clade': 'The clade to search for',
        },
    output_descriptions={
        'characterization': "The clade characterization detail"},
    name='Clade ASV level assessment',
    description=("Test whether the clade is uniquely represented by ASVs "
                 "and other summary information about the clade"),
    citations=[]
)


plugin.methods.register_function(
    function=q2_gg2.sequence_v4_asv_assessment,
    inputs={'phylogeny': Phylogeny[Rooted],
            'taxa': FeatureData[Taxonomy],
            'full_length_v4': FeatureData[Sequence],
            'sequences': FeatureData[Sequence]},
    parameters={'asv': Str},
    outputs=[('characterization', ASVAssessment)],
    input_descriptions={
        'phylogeny': 'The ID annotated phylogeny',
        'taxa': 'The taxa indexed by ID',
        'full_length_v4': 'The extracted V4 sequences from the backbone',
        'sequences': 'The full set of GG2 sequences'},
    parameter_descriptions={
        'asv': 'The ASV sequence to summarize',
        },
    output_descriptions={
        'characterization': "The ASV characterization detail"},
    name='ASV assessment',
    description=("Check what full length records the ASV associates with "
                 "and members within its placement multifurcation"),
    citations=[]
)


plugin.methods.register_function(
    function=q2_gg2.bulk_sequence_v4_asv_assessment,
    inputs={'phylogeny': Phylogeny[Rooted],
            'taxa': FeatureData[Taxonomy],
            'full_length_v4': FeatureData[Sequence],
            'sequences': FeatureData[Sequence]},
    parameters={'group': Int,
                'version': Str,
                'output_filename': Str},
    outputs=[('characterization', ASVAssessment)],
    input_descriptions={
        'phylogeny': 'The ID annotated phylogeny',
        'taxa': 'The taxa indexed by ID',
        'full_length_v4': 'The extracted V4 sequences from the backbone',
        'sequences': 'The full set of GG2 sequences'},
    parameter_descriptions={
        'group': 'The hash group to process',
        'version': 'The database version to use',
        'output_filename': 'The file to write assessments too'
        },
    output_descriptions={
        'characterization': "This output is undefined"},
    name='Bulk ASV assessment',
    description=("Check what full length records the ASV associates with "
                 "and members within its placement multifurcation"),
    citations=[]
)


plugin.methods.register_function(
    function=q2_gg2.bulk_clade_v4_asv_assessment,
    inputs={'phylogeny': Phylogeny[Rooted],
            'taxa': FeatureData[Taxonomy],
            'full_length_v4': FeatureData[Sequence],
            'sequences': FeatureData[Sequence]},
    parameters={'group': Int,
                'version': Str,
                'output_filename': Str},
    outputs=[('characterization', CladeAssessment)],
    input_descriptions={
        'phylogeny': 'The ID annotated phylogeny',
        'taxa': 'The taxa indexed by ID',
        'full_length_v4': 'The extracted V4 sequences from the backbone',
        'sequences': 'The full set of GG2 sequences'},
    parameter_descriptions={
        'group': 'The hash group to process',
        'version': 'The database version to use',
        'output_filename': 'The file to write assessments too'
        },
    output_descriptions={'characterization': 'This output is undefined'},
    name='Bulk Clade ASV level assessment',
    description=("Test whether the clade is uniquely represented by ASVs "
                 "and other summary information about the clade"),
    citations=[]
)


plugin.methods.register_function(
    function=q2_gg2.clade_lookup,
    inputs={'taxonomy_as_tree': Phylogeny[Rooted]},
    parameters={'version': Str,
                'output_filename': Str},
    outputs=[('characterization', CladeAssessment)],
    input_descriptions={
        'taxonomy_as_tree': 'The taxonomy represented as a newick string'},
    parameter_descriptions={
        'version': 'The database version to use',
        'output_filename': 'The file to write the lookup too'
        },
    output_descriptions={'characterization': 'This output is undefined'},
    description=("Construct a taxonomy lookup suitable for a sqlite database"),
    citations=[],
    name='Construct a clade lookup'
)


plugin.methods.register_function(
    function=q2_gg2.taxonomy_from_table,
    inputs={'reference_taxonomy': Phylogeny[Rooted],
            'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[('classification', FeatureData[Taxonomy])],
    input_descriptions={
        'reference_taxonomy': 'The reference taxonomy to derive from',
        'table': 'The feature table to classify'
    },
    parameter_descriptions={},
    output_descriptions={
        'classification': 'The resulting classifications',
    },
    name='Classify features against Greengenes2',
    description=("Pull lineage information for each feature off the"
                 "reference phylogey"),
    citations=[]
)


plugin.methods.register_function(
    function=q2_gg2.taxonomy_from_features,
    inputs={'reference_taxonomy': Phylogeny[Rooted],
            'reads': FeatureData[Sequence]},
    parameters={},
    outputs=[('classification', FeatureData[Taxonomy])],
    input_descriptions={
        'reference_taxonomy': 'The reference taxonomy to derive from',
        'reads': 'The feature data to classify'
    },
    parameter_descriptions={},
    output_descriptions={
        'classification': 'The resulting classifications',
    },
    name='Classify ASVs against Greengenes2',
    description=("Pull lineage information for each feature off the"
                 "reference phylogey"),
    citations=[]
)


plugin.methods.register_function(
    function=q2_gg2.filter_features,
    inputs={'feature_table': FeatureTable[Frequency],
            'reference': Phylogeny[Rooted]},
    parameters={},
    outputs=[('filtered_feature_table', FeatureTable[Frequency])],
    input_descriptions={
        'feature_table': "The feature table to filter",
        'reference': "The reference phylogeny or taxonomy to filter against"},
    parameter_descriptions={},
    output_descriptions={
        'filtered_feature_table': 'The resulting filtered feature table'},
    name='Filter features against reference',
    description=("Filter the features in a feature table against the feature "
                 "set available from the reference"),
    citations=[]
)


plugin.methods.register_function(
    function=q2_gg2.relabel,
    inputs={'feature_table': FeatureTable[Frequency],
            'reference_label_map': ReferenceMap},
    parameters={'as_md5': Bool,
                'as_asv': Bool,
                'as_id': Bool},
    outputs=[('relabeled_table', FeatureTable[Frequency])],
    input_descriptions={
        'feature_table': "The feature table to relabel",
        'reference_label_map': "The reference label mapping"},
    parameter_descriptions={
        'as_md5': 'Convert to md5 labels',
        'as_asv': ('Convert to ASV labels, assumes ASVs are < 500nt, will '
                   'use sequence IDs as for non ASVs'),
        'as_id': 'Convert to identifiers'
    },
    output_descriptions={
        'relabeled_table': 'The resulting relabeled feature table'},
    name='Relabel features in a feature table',
    description=("Relabels the features in a feature table with a provided"
                 "reference mapping"),
    citations=[]
)
importlib.import_module('q2_gg2._transformer')
