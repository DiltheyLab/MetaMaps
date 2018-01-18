use strict;
use List::MoreUtils qw/all mesh any /;
use List::Util qw/sum min max/;
use Data::Dumper;
use Getopt::Long;   
use File::Find;
use Math::GSL::Randist qw/gsl_ran_binomial_pdf/;
use FindBin;
use lib "$FindBin::Bin/../perlLib";
use Cwd qw/getcwd abs_path/;
use File::Copy;
use Storable qw/dclone store retrieve/;

use taxTree;
use validation;
use Util;

$| = 1;

my $masterTaxonomy_dir = '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';
my $DB = 'databases/miniSeq+H';
# my $MetaMap_results = '/scratch/tmp/MetaMap/hmp_set7';
my $MetaMap_results = 'tmp/hmp7_2_miniSeq+H';
my $kraken_results_dir = '/scratch/tmp/hmp_set7_combined_kraken_results';
my $truth = 'tmp/truthHMP7_blasr';


my %results_readLevel = (
	'MetaMap-EM' => $MetaMap_results . '.EM.reads2Taxon',
	'MetaMap-U' => $MetaMap_results . '.U.reads2Taxon',
	'Kraken' => $kraken_results_dir . '/results_kraken.txt.reads2Taxon',
);

my %results_distribution = (
	'MetaMap-EM' => $MetaMap_results . '.EM.WIMP',
	'MetaMap-U' => $MetaMap_results . '.U.WIMP',
	'Kraken' => $kraken_results_dir . '/results_kraken.txt',
	'Bracken' => $kraken_results_dir . '/results_bracken.txt',
);

die Dumper("Missing files", [grep {not -e $_} values %results_readLevel]) unless(all {-e $_} values %results_readLevel);
die unless(all {-e $_} values %results_distribution);


my $master_taxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
my $master_taxonomy_merged = taxTree::readMerged($masterTaxonomy_dir);

my $taxonomy_usedForInference = taxTree::readTaxonomy($DB . '/taxonomy');
(my $extendedMaster, my $extendedMaster_merged) = validation::prepare_masterTaxonomy_withX($masterTaxonomy_dir, $taxonomy_usedForInference);

my $truth_reads_href = validation::readTruthFileReads($extendedMaster, $extendedMaster_merged, $truth . '.perRead');

my $readLengths_href = Util::getReadLengths('/scratch/tmp/hmp_set7_combined.fastq');
my $mappableTaxonomy;
my $truth_reads_mappable;
my %reduced_taxonID_master_2_contigs;
{
		
	my %reduced_taxonID_original_2_contigs;
	my %reduced_contigLength;
	Util::read_taxonIDs_and_contigs($DB, \%reduced_taxonID_original_2_contigs, \%reduced_contigLength);
	
	# read reduced taxonomy
	foreach my $taxonID_original (keys %reduced_taxonID_original_2_contigs)
	{
		my $taxonID_master = taxTree::findCurrentNodeID($extendedMaster, $extendedMaster_merged, $taxonID_original);
		$reduced_taxonID_master_2_contigs{$taxonID_master} = 1;
	}
	
	$mappableTaxonomy = dclone $extendedMaster;
	taxTree::removeUnmappableParts($mappableTaxonomy, \%reduced_taxonID_master_2_contigs);

	# translate truth into reduced representation
	$truth_reads_mappable = validation::translateReadsTruthToReducedTaxonomy($extendedMaster, $mappableTaxonomy, $truth_reads_href);
}

my $truth_reads_href_noUnknown = { map {$_ => $truth_reads_mappable->{$_}} grep {$truth_reads_mappable->{$_}} keys %$truth_reads_mappable };
die unless(all {exists $mappableTaxonomy->{$_}} values %$truth_reads_href_noUnknown);

my $truth_mappingDatabase_distribution = validation::truthReadsToTruthSummary($mappableTaxonomy, $truth_reads_href_noUnknown, \%reduced_taxonID_master_2_contigs);


my $n_reads_correct_byVariety = {};
my $n_reads_correct_byVariety_byLevel = {};
my $freq_byVariety_byLevel = {};
my $n_reads_correct_byVariety_byLevel_byLength = {};
	
foreach my $label (keys %results_readLevel)
{
	print "Analysing $label -- at level of individual reads!\n";
	
	$n_reads_correct_byVariety->{$label} = {} unless(defined $n_reads_correct_byVariety->{$label});
	$n_reads_correct_byVariety_byLevel->{$label} = {} unless(defined $n_reads_correct_byVariety_byLevel->{$label});
	$n_reads_correct_byVariety_byLevel_byLength->{$label} = {} unless(defined $n_reads_correct_byVariety_byLevel_byLength->{$label});
	

	my $inferred_reads = validation::readInferredFileReads($extendedMaster, $extendedMaster_merged, $results_readLevel{$label});

	my @readIDs_no_truth = grep {not exists $truth_reads_href->{$_}} keys %$inferred_reads;
	if(scalar(@readIDs_no_truth))
	{
		die "Error: have ", scalar(@readIDs_no_truth), " (of ", scalar(keys %$inferred_reads), ") reads without defined truth.\n";
	}
	
	validation::readLevelComparison($extendedMaster, $truth_reads_href, $truth_reads_mappable, $inferred_reads, $label, $n_reads_correct_byVariety->{$label}, $n_reads_correct_byVariety_byLevel->{$label}, $n_reads_correct_byVariety_byLevel_byLength->{$label}, \%reduced_taxonID_master_2_contigs, $readLengths_href);
}

my @evaluateAccuracyAtLevels = validation::getEvaluationLevels();

my %distributions_byLevel_byLabel;

my %union_taxonIDs_byLevel;

foreach my $level (@evaluateAccuracyAtLevels)
{
	die unless(defined $truth_mappingDatabase_distribution->{$level});
	foreach my $trueTaxonID (keys %{$truth_mappingDatabase_distribution->{$level}})
	{
		$union_taxonIDs_byLevel{$level}{$trueTaxonID}++;
		$distributions_byLevel_byLabel{$level}{'truth'}{$trueTaxonID} = $truth_mappingDatabase_distribution->{$level}{$trueTaxonID};
	}	
}

foreach my $label (keys %results_distribution)
{
	print "Analysing $label -- at level of distribution!\n";

	$freq_byVariety_byLevel->{$label} = {} unless(defined $freq_byVariety_byLevel->{$label});

	my $inferred_distribution = validation::readInferredDistribution($extendedMaster, $extendedMaster_merged, $results_distribution{$label});	

	
	foreach my $level (@evaluateAccuracyAtLevels)
	{		
		if(defined $inferred_distribution->{$level})
		{
			foreach my $inferredTaxonID (keys %{$inferred_distribution->{$level}})
			{
				$union_taxonIDs_byLevel{$level}{$inferredTaxonID}++;
				$distributions_byLevel_byLabel{$level}{$label}{$inferredTaxonID} = $inferred_distribution->{$level}{$inferredTaxonID}[1];
			}
		}
	}	
}

my $fn_output = '_HMP_distributions.txt';
open(F, '>', $fn_output) or die "Cannot open $fn_output";
print F join("\t", "Level", "Source", "taxonID", "taxonLabel", "F"), "\n";
foreach my $level (keys %union_taxonIDs_byLevel)
{
	foreach my $label (keys %{$distributions_byLevel_byLabel{$level}})
	{
		foreach my $taxonID (keys %{$union_taxonIDs_byLevel{$level}})
		{
			$distributions_byLevel_byLabel{$level}{$label}{$taxonID} = 0 if(not defined $distributions_byLevel_byLabel{$level}{$label}{$taxonID});
			$distributions_byLevel_byLabel{$level}{$label}{$taxonID} = 0 if(not defined $distributions_byLevel_byLabel{$level}{$label}{$taxonID});
		}
		foreach my $taxonID (keys %{$distributions_byLevel_byLabel{$level}{$label}})
		{
			my $taxonLabel;
			my $taxonID2 = $taxonID;
			if($taxonID eq 'Unclassified')
			{
				$taxonID2 = 0;
				$taxonLabel = 'Unclassified';
			}
			elsif($taxonID eq 'Undefined')
			{
				$taxonLabel = 'Undefined';
			}
			elsif($taxonID eq 'NotLabelledAtLevel')
			{
				$taxonLabel = 'NotLabelledAtLevel';
			}			
			else
			{
				$taxonLabel = taxTree::taxon_id_get_name($taxonID, $master_taxonomy);			
			}
			print F join("\t", $level, $label, $taxonID2, $taxonLabel, $distributions_byLevel_byLabel{$level}{$label}{$taxonID}), "\n";
		}
	}
}

			
		
		