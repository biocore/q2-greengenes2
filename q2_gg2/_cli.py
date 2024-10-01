import click
import qiime2
import pandas as pd
from q2_types.tree import NewickFormat
from q2_types.feature_data import DNAFASTAFormat
from ._methods import (
    bulk_sequence_v4_asv_assessment as _bulk_sequence_v4_asv_assessment,
)
from ._methods import bulk_clade_v4_asv_assessment as _bulk_clade_v4_asv_assessment
from ._methods import clade_lookup as _clade_lookup


@click.group()
def cli():
    pass


@cli.command()
@click.option("--i-taxonomy-as-tree", type=click.Path(exists=True), required=True)
@click.option("--p-output-filename", type=click.Path(exists=False), required=True)
@click.option("--p-version", type=str, required=True)
def clade_lookup(i_taxonomy_as_tree, p_output_filename, p_version):
    taxonomy_as_tree = qiime2.Artifact.load(i_taxonomy_as_tree).view(NewickFormat)
    _clade_lookup(taxonomy_as_tree, p_version, p_output_filename)


@cli.command()
@click.option("--i-phylogeny", type=click.Path(exists=True), required=True)
@click.option("--i-taxa", type=click.Path(exists=True), required=True)
@click.option("--i-full-length-v4", type=click.Path(exists=True), required=True)
@click.option("--i-sequences", type=click.Path(exists=True), required=True)
@click.option("--p-group", type=int, required=True)
@click.option("--p-output-filename", type=click.Path(exists=False), required=True)
@click.option("--p-version", type=str, required=True)
def bulk_clade_v4_asv_assessment(
    i_phylogeny,
    i_taxa,
    i_full_length_v4,
    i_sequences,
    p_group,
    p_output_filename,
    p_version,
):
    phylogeny = qiime2.Artifact.load(i_phylogeny).view(NewickFormat)
    taxa = qiime2.Artifact.load(i_taxa).view(pd.DataFrame)
    full_length_v4 = qiime2.Artifact.load(i_full_length_v4).view(DNAFASTAFormat)
    sequences = qiime2.Artifact.load(i_sequences).view(DNAFASTAFormat)
    _bulk_clade_v4_asv_assessment(
        phylogeny,
        taxa,
        full_length_v4,
        sequences,
        p_version,
        p_group,
        p_output_filename,
    )


@cli.command()
@click.option("--i-phylogeny", type=click.Path(exists=True), required=True)
@click.option("--i-taxa", type=click.Path(exists=True), required=True)
@click.option("--i-full-length-v4", type=click.Path(exists=True), required=True)
@click.option("--i-sequences", type=click.Path(exists=True), required=True)
@click.option("--o-characterization", type=click.Path(exists=False), required=True)
@click.option("--p-group", type=int, required=True)
@click.option("--p-output-filename", type=click.Path(exists=False), required=True)
@click.option("--p-version", type=str, required=True)
def bulk_sequence_v4_asv_assessment(
    i_phylogeny,
    i_taxa,
    i_full_length_v4,
    i_sequences,
    o_characterization,
    p_group,
    p_output_filename,
    p_version,
):
    phylogeny = qiime2.Artifact.load(i_phylogeny).view(NewickFormat)
    taxa = qiime2.Artifact.load(i_taxa).view(pd.DataFrame)
    full_length_v4 = qiime2.Artifact.load(i_full_length_v4).view(DNAFASTAFormat)
    sequences = qiime2.Artifact.load(i_sequences).view(DNAFASTAFormat)
    _bulk_sequence_v4_asv_assessment(
        phylogeny,
        taxa,
        full_length_v4,
        sequences,
        p_version,
        p_group,
        p_output_filename,
    )
