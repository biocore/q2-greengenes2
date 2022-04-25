# Background

The Greengenes2 phylogeny is based on whole genome information from the [Web of Life](https://biocore.github.io/wol/), and revised with high quality full length 16S from the [Living Tree Project](https://imedea.uib-csic.es/mmg/ltp/) and full length 16S extracted from [bacterial operons](https://www.nature.com/articles/s41592-020-01041-y) using [uDance](https://github.com/balabanmetin/uDance). A seed taxonomy is derived using the mappings from the Web of Life to [GTDB](https://gtdb.ecogenomic.org/). This taxonomy is then augmented using information from the Living Tree Project when possible. The augmented taxonomy is decorated onto the backbone using [tax2tree](https://github.com/biocore/tax2tree).

Using this decorated backbone, all public and private 16S V4 ASVs from [Qiita](https://qiita.ucsd.edu/) pulled from [redbiom](https://github.com/biocore/redbiom/) representing hundreds of thousands of samples, as well as full length mitochondrial and chloroplast 16S (sourced from [SILVA](https://www.arb-silva.de/), are then placed using [DEPP](https://github.com/yueyujiang/DEPP). Fragments are resolved. The resulting tree contains > 15,000,000 tips. 

Fragment resolution can result in fragments being placed on the parent edge of a named node. This can occur if the node representing a clade, such as d__Archaea, does not represent sufficient diversity for the input fragments to place. As a result, prior to reading taxonomy off of the tree, each name from the backbone is evaluated for whether its edge to parent has a single or multifurcation of placements. If this occurs, the name is “promoted”. The idea being that fragments off a named edge to its parent are more like the named node than a sibling.

Following this name promotion, the full taxonomy is then read off the tree providing lineage information for each fragment and sequence represented in the tree. This taxonomy information can be utilized within QIIME 2 by cross referencing your input feature set against what’s present in the tree. By doing so, we can obtain taxonomy for both WGS data (if processed by Woltka) and 16S V4 ASVs. There is an important caveat though: right now, we can only classify based sequences already represented by the tree, so unrepresented V4 ASVs will be unassigned.  

# Install

```
$ git clone https://github.com/wasade/q2-greengenes2.git
$ source activate qiime2.2022.2
$ cd q2-greengenes2
$ conda install -c bioconda iow
$ pip install -e .
```

# Reference database artifacts

The reference database release contains the following artifacts. <version> refers to the version of the database, which follows a YYYY.MM format.  

Unless otherwise indicated, the feature IDs present in the artifacts use the WoL namespace for genomes, and md5 hashes for all else. A label mapping artifact is provided should transformations to the original ID, hash or raw sequence be necessary.

* `gg2-<version>-taxonomy.newick.qza`
	A Newick version of the taxonomy where internal nodes are taxa, and tips are feature IDs
* `gg2-<version>-taxonomy.tsv.qza`
	The taxonomy data in the classic tab delimited format
* `gg2-<version>-phylogeny.qza`
	A phylogeny representing all of the Greengenes2 features
* `gg2-<version>-label_map.qza`
	A tab delimited map of sequence ID, md5 hash of a sequence, and the sequence itself
* `gg2-<version>-sequences.fna.gz`
	A FASTA file of the records used in Greengenes2

# To classify

Feature tables which contain either 16S V4, WoL genome identifiers, or a combination thereof can be taxonomically characterized against Greengenes using the following command:

```
$ qiime greengenes2 taxonomy-from-table \
	--i-reference-taxonomy <the_greengenes_reference> \
	--i-table <your_feature_table> \
    --o-classification <the_resulting_classifications>
```

The QIIME 2 Greengenes2 plugin also supports the classic method of classification through `FeatureData[Sequence]` artifacts:

```
$ qiime greengenes2 taxonomy-from-features \
    --i-reference-taxonomy <the_greengenes_reference> \
    --i-reads <your_ASVs> \ 
    --o-classification <the_resulting_classifications>
```

# Diversity calculations

The Greengenes2 reference tree can be readily used for feature tables based on WoL-assessed short read data, 16S V4 ASVs, or the combination of those data types. The only essential requirement is that the features represented by the table are also present in the tree.

The Greengenes2 plugin implements a rapid method for performing this filtering. The same filtering is also possible using the `q2-fragment-insertion` plugin, however that plugin does not presently use faster phylogeny parsing logic. For filtering, either the phylogeny or taxonomy tree can be used:

```
$ qiime greengenes2 filter-features \
    --i-feature-table <your_feature_table> \
    --i-reference <a_greengenes_tree> \
    --o-filtered-table <the_filtered_table>
```

From here, other phylogenetic-aware methods such as UniFrac can be performed as normal.

# Relabeling

The Greengenes2 taxonomy and phylogeny can be expressed in three different namespaces:

* Record IDs, such as Genbank accessions
* md5 sequence hashes 
* The actual sequence

If your input `FeatureTable[Frequency]` was computed using default parameters from `q2-dada2` or `q2-deblur`, then it is likely the features are expressed as md5 hashes. However, if your table is coming from a `redbiom` query or Qiita, then your features are most likely expressed as sequences and/or Woltka record IDs. 

The q2-gg2 plugin provides a way to relabel your table. For example, the following will relabel as "md5":

```
$ qiime greengenes2 relabel \
    --i-feature-table <your_feature_table> \
    --i-reference-label-map gg2-<version>-label_map.qza \
    --p-as-md5 \
    --o-relabeled-table <the_relabeled_table>
```
