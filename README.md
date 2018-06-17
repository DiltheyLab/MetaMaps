MetaMaps
========================================================================

MetaMaps is tool specifically developed for the analysis of long-read (PacBio/ONT) metagenomic datasets.

It simultaenously carries out read assignment and sample composition estimation.

It is faster than classical exact alignment-based approaches, and its output is more information-rich than that of kmer-spectra-based methods. For example, each MetaMaps alignment comes with approximate alignment locations, and estimated alignment identity and a mapping quality.

The approximate mapping algorithm employed by MetaMaps is identical to that of [https://github.com/marbl/MashMap](MashMap).


## Installation
Follow [`INSTALL.txt`](INSTALL.txt) to compile and install MetaMaps.

Then download a database, e.g. miniSeq+H (microbial genomes and the human reference genome). Extract the downloaded database into the `databases/` directory.

## Usage

Analysis of a dataset with MetaMap consists of two steps: mapping and classification:

```
./metamaps mapDirectly --all -r databases/miniSeq+H/DB.fa -q input.fastq -o classification_results
./metamaps classify --mappings classification_results --DB databases/miniSeq+H
```

### Memory-efficient mapping

You can use the '--maxmemory' parameter to specify a target for maximum memory consumption (in gigabytes). Note that this feature is implemented heuristically; actual memory usage will fluctuate around and may exceed the target. It is recommended to use around 70% of the available memory as a target amount.

Example:

```
./metamaps mapDirectly --all -r databases/miniSeq+H/DB.fa -q input.fastq -o classification_results --maxmemory 20
./metamaps classify --mappings classification_results --DB databases/miniSeq+H
```

## Databases

The 'miniSeq+H' database is a good place to start. It contains >12000 microbial genomes and the human reference genome. We provide miniSeq+H as a download.

You can also download and construct your own reference databases. For example, this is how to construct the miniSeq+H database:

1. Download the genomes you want to include. The easiest way to do this is by copying the RefSeq/Genbank directory structure of the taxonomic branches you're interested in. This can be done with the `downloadRefSeq.pl` script, which is easily customizable (e.g., `--targetBranches archaea,bacteria,fungi` to download these three branches). Example:

    ```
    mkdir testDownload
    perl downloadRefSeq.pl --seqencesOutDirectory testDownload/refseq --taxonomyOutDirectory testDownload/taxonomy
    ```

2. We need to make sure that each contig ID is annotated with a correct and unique taxon ID and we want the whole database as one file. `annotateRefSeqSequencesWithUniqueTaxonIDs.pl` can help:

    ```
    perl annotateRefSeqSequencesWithUniqueTaxonIDs.pl --refSeqDirectory downloads/refseq --taxonomyInDirectory downloads/taxonomy --taxonomyOutDirectory downloads/taxonomy_uniqueIDs
    ```
    
3. We might also manually want to include additional genomes, for example the human reference genome. Obtain the genome in one file (e.g. `hg38.primary.fna`) and add taxon IDs:

    ```
    perl util/addTaxonIDToFasta.pl --inputFA hg38.primary.fna --outputFA hg38.primary.fna.with9606 --taxonID 9606
    ```
    
4. Finally, construct the MetaMap database:

    ```
    perl buildDB.pl --DB databases/myDB --FASTAs downloads/refseq,hg38.primary.fna.with9606 --taxonomy downloads/taxonomy_uniqueIDs
    ```
	
    The NCBI taxonomy changes on a regular basis, and you might not want to repeat the complete database construction process every time that happens. You can update the utilized taxonomy as part of `buildDB.pl`, by specifying the "old" taxonomy (used for `addTaxonIDToFasta.pl`), `--updateTaxonomy 1`, and the path to a download of the new taxonomy (e.g. [ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)). Example:

    ```
    perl buildDB.pl --DB databases/myDB --FASTAs downloads/refseq/ref.fa,hg38.primary.fna.with9606 --taxonomy downloads/new_taxonomy --oldTaxonomy downloads/taxonomy_uniqueIDs --updateTaxonomy 1
    ```




