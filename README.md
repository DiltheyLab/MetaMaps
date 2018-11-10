MetaMaps
========================================================================

MetaMaps is tool specifically developed for the analysis of long-read (PacBio/ONT) metagenomic datasets.

It simultaenously carries out read assignment and sample composition estimation.

It is faster than classical exact alignment-based approaches, and its output is more information-rich than that of kmer-spectra-based methods. For example, each MetaMaps alignment comes with an approximate alignment location, an estimated alignment identity and a mapping quality.

The approximate mapping algorithm employed by MetaMaps is based on [MashMap](https://github.com/marbl/MashMap). MetaMaps adds a mapping quality model and EM-based estimation of sample composition.

## News
(28 August 2018) We're adding more flexible tools to construct your own databases - see [here](https://github.com/DiltheyLab/MetaMaps/issues/5#issuecomment-414301437) for details. Feedback welcome!

## Installation
Follow [`INSTALL.txt`](INSTALL.txt) to compile and install MetaMaps.

Then download a database, e.g. [miniSeq+H](https://www.dropbox.com/s/g2jzj8sklnlp3j2/miniSeq%2BH.tar.gz?dl=0) (~8G compressed, microbial genomes and the human reference genome). Extract the downloaded database into the `databases/` directory.

## Usage

Analysis of a dataset with MetaMaps consists of two steps: mapping and classification:

```
./metamaps mapDirectly --all -r databases/miniSeq+H/DB.fa -q input.fastq -o classification_results
./metamaps classify --mappings classification_results --DB databases/miniSeq+H
```

### Memory-efficient mapping

You can use the '--maxmemory' parameter to specify a target for maximum memory consumption (in gigabytes). Note that this feature is implemented heuristically; actual memory usage will fluctuate and execeed the target. We recommend using around 70% of the available memory as a target amount (for example, 20G when you have a 32G machine).

Example:

```
./metamaps mapDirectly --all -r databases/miniSeq+H/DB.fa -q input.fastq -o classification_results --maxmemory 20
./metamaps classify --mappings classification_results --DB databases/miniSeq+H
```

## Output

MetaMaps outputs both an overall compositional assignment and per-read taxonomic assignments. Specifically, it will (for `-o classification_results`) produce the following files:

1. `classification_results.EM.WIMP`: Sample composition at different taxonomic levels (WIMP = "What's in my pot"). The level "definedGenomes" represents strain-level resolution (i.e., the defined genomes in the classification database). The EM algorithm is carried out at this level.
   
   Output columns: `Absolute` specifies the number of reads assigned (by their maximum likelihood mapping estimate) to the taxonomic entity; `EMFrequency` specifies the estimated frequency of the taxonomic entity prior to taking into account unmapped reads; `PotFrequency` specifies the estimated final frequency of the taxonomic entity (i.e. after correcting for unmapped reads).

2. `classification_results.EM.reads2Taxon`: One line per read, and each line consists of the read ID and the taxon ID of the genome that the read was assigned to. Taxon IDs beginning with an 'x' represent MetaMaps-internal taxon IDs that disambiguate between source genomes attached to the same 'species' node. These can be interpreted using the extended database taxonomy (sub-directory `taxonomy` in the directory of the utilized database).

3. `classification_results.EM.reads2Taxon.krona`: Like `classification_results.EM.reads2Taxon`, but in [Krona](https://github.com/marbl/Krona/wiki) format. Each line has an additional quality value, and only taxon IDs from the standard NCBI taxonomy are used.

4. `classification_results.EM.contigCoverage`: Read coverage for contigs that appear in the final set of maximum-likelihood mappings. Contigs are split into windows of 1000 base pairs. Each line represents one window and specifies the MetaMaps taxonID of the contig (`taxonID`), a species/strain label (`equalCoverageUnitLabel`), the ID of the contig in the classification database FASTA file (`contigID`), start and stop coordinates of the window (`start` and `stop`), the number of read bases overlapping the window (`nBases`), and the number of read bases overlapping the window divided by window length, i.e. coverage (`readCoverage`).

5. `classification_results.EM.lengthAndIdentitiesPerMappingUnit`: Read length and estimated identity for all reads, sorted by the contig they are mapped to. Each line represents one read and has the fields `AnalysisLevel`, which is always equal to `EqualCoverageUnit`; `ID`, which is the ID of the contig in the classification database FASTA file; `readI`, which is a simple numerical read ID; `Identity`, which is the estimated alignment identity; and `Length`, specifying the length of the read.

5. `classification_results.EM`: The final and complete set of approximate read mappings. Each line represents one mapping in a simple space-delimited format. Fields: read ID, read length, beginning of the mapping in the read, end of the mapping in the read, strand, contig ID, contig length, beginning of the mapping in the contig, end of the mapping in the contig, estimated alignment identity using a Poisson model, MinHash intersection size, MinHash union size, estimated alignment identity using a binomial approximation, mapping quality. The mapping qualities of all mappings for one read sum up to 1.

6. `classification_results.EM.evidenceUnknownSpecies`: Various statistics on read identities and zero-coverage regions of identified genomes. These are not particularly useful in their current state; visual inspection of coverage and identity patterns should take precedence.

You can download [example output files](https://github.com/DiltheyLab/MetaMaps/blob/master/MetaMaps_example_output.zip).

## Databases

The 'miniSeq+H' database is a good place to start. It contains >12000 microbial genomes and the human reference genome. We provide miniSeq+H as a download (see above for the link).

You can also download and construct your own reference databases. For example, this is how to construct the miniSeq+H database:

1. Download the genomes you want to include. The easiest way to do this is by copying the RefSeq/Genbank directory structure of the taxonomic branches you're interested in. This can be done with the `downloadRefSeq.pl` script, which is easily customizable (e.g., `--targetBranches archaea,bacteria,fungi` to download these three branches). Example:

    ```
    mkdir download
	perl downloadRefSeq.pl --seqencesOutDirectory download/refseq --taxonomyOutDirectory download/taxonomy
    ```

2. We need to make sure that each contig ID is annotated with a correct and unique taxon ID and we want the whole database as one file. `annotateRefSeqSequencesWithUniqueTaxonIDs.pl` can help:

    ```
    perl annotateRefSeqSequencesWithUniqueTaxonIDs.pl --refSeqDirectory download/refseq --taxonomyInDirectory download/taxonomy --taxonomyOutDirectory download/taxonomy_uniqueIDs
    ```
    
    By default this script will only process complete (finished) assemblies - if you want to modify this behaviour, uncomment the line `next unless($assembly_level eq 'Complete Genome');`.
    
3. We might also manually want to include additional genomes, for example the human reference genome. Obtain the genome in one file (e.g. `hg38.primary.fna`) and add taxon IDs:

    ```
    perl util/addTaxonIDToFasta.pl --inputFA hg38.primary.fna --outputFA hg38.primary.fna.with9606 --taxonID 9606
    ```
    
4. Finally, construct the MetaMaps databasen (here `myDB`):

    ```
    perl buildDB.pl --DB databases/myDB --FASTAs download/refseq,hg38.primary.fna.with9606 --taxonomy download/taxonomy_uniqueIDs
    ```
	
    The NCBI taxonomy changes on a regular basis, and you might not want to repeat the complete database construction process every time that happens. You can update the utilized taxonomy as part of `buildDB.pl`, by specifying the "old" taxonomy (used for `addTaxonIDToFasta.pl`), `--updateTaxonomy 1`, and the path to a download of the new taxonomy (e.g. [ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)). Example:

    ```
    perl buildDB.pl --DB databases/myDB --FASTAs download/refseq/ref.fa,hg38.primary.fna.with9606 --taxonomy download/new_taxonomy --oldTaxonomy download/taxonomy_uniqueIDs --updateTaxonomy 1
    ```




