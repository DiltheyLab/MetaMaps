use strict;
use List::MoreUtils qw/all mesh any /;
use List::Util qw/sum min max/;
use Data::Dumper;
use Getopt::Long;   
use File::Find;
use Math::GSL::Randist qw/gsl_ran_binomial_pdf/;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;
use File::Copy;
use Storable qw/dclone store retrieve/;

use taxTree;
use validation;
use Util;

$| = 1;

my $masterTaxonomy_dir = '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';
 
my $varietyName = 'main';
my $methodName = 'MetaMap';
my $mappings;
my $truth;
my $database;
my $projectTruthOntoDBTaxonomy = 1;
GetOptions (
	'mappings:s' => \$mappings, 
	'DB:s' => \$database, 
	'truth:s' => \$truth, 
	'projectTruthOntoDBTaxonomy:s' => \$projectTruthOntoDBTaxonomy, 
);

die "Please specify --mappings" unless($mappings);
die "Please specify --truth" unless($truth);
die "Please specify --DB" unless($database);

die "Expected file not found: ${truth}.perRead" unless(-e $truth . '.perRead');
die "Expected file not found: ${mappings}.WIMP" unless(-e $mappings . '.WIMP');
die "Expected file not found: ${mappings}.reads2Taxon" unless(-e $mappings . '.reads2Taxon');

my $master_taxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
my $master_taxonomy_merged = taxTree::readMerged($masterTaxonomy_dir);

my $taxonomy_usedForInference = taxTree::readTaxonomy($database . '/taxonomy');
(my $extendedMaster, my $extendedMaster_merged) = validation::prepare_masterTaxonomy_withX($masterTaxonomy_dir, $taxonomy_usedForInference);


my $truth_reads_href = validation::readTruthFileReads($extendedMaster, $extendedMaster_merged, $truth . '.perRead');

my $specificTaxonomy;
my $truth_mappingDatabase_reads;
if($projectTruthOntoDBTaxonomy)
{
	my %reduced_taxonID_original_2_contigs;
	my %reduced_contigLength;
	Util::read_taxonIDs_and_contigs($database, \%reduced_taxonID_original_2_contigs, \%reduced_contigLength);
	
	# read reduced taxonomy
	my %reduced_taxonID_master_2_contigs;
	foreach my $taxonID_original (keys %reduced_taxonID_original_2_contigs)
	{
		my $taxonID_master = taxTree::findCurrentNodeID($extendedMaster, $extendedMaster_merged, $taxonID_original);
		$reduced_taxonID_master_2_contigs{$taxonID_master} = 1;
	}
	
	$specificTaxonomy = dclone $extendedMaster;
	taxTree::removeUnmappableParts($specificTaxonomy, \%reduced_taxonID_master_2_contigs);

	# translate truth into reduced representation
	$truth_mappingDatabase_reads = validation::translateReadsTruthToReducedTaxonomy($extendedMaster, $specificTaxonomy, $truth_reads_href);
}
else
{
	# warn Dumper($extendedMaster);
	$specificTaxonomy = dclone $extendedMaster;
	$truth_mappingDatabase_reads = $truth_reads_href;
}

my $truth_reads_href_noUnknown = { map {$_ => $truth_mappingDatabase_reads->{$_}} grep {$truth_mappingDatabase_reads->{$_}} keys %$truth_mappingDatabase_reads };
die unless(all {exists $specificTaxonomy->{$_}} values %$truth_reads_href_noUnknown);

my $truth_mappingDatabase_distribution = validation::truthReadsToTruthSummary($specificTaxonomy, $truth_reads_href_noUnknown);

my $inferred_reads = validation::readInferredFileReads($extendedMaster, $extendedMaster_merged, $mappings . '.reads2Taxon');
my $inferred_distribution = validation::readInferredDistribution($extendedMaster, $extendedMaster_merged, $mappings . '.WIMP');	


my $n_reads_correct_byVariety = {};
my $n_reads_correct_byVariety_byLevel = {};
my $freq_byVariety_byLevel = {};
	

$n_reads_correct_byVariety->{$varietyName} = {} unless(defined $n_reads_correct_byVariety->{$varietyName});
$n_reads_correct_byVariety_byLevel->{$varietyName} = {} unless(defined $n_reads_correct_byVariety_byLevel->{$varietyName});
$freq_byVariety_byLevel->{$varietyName} = {} unless(defined $freq_byVariety_byLevel->{$varietyName});

my @readIDs_no_truth = grep {not exists $truth_reads_href->{$_}} keys %$inferred_reads;
if(scalar(@readIDs_no_truth))
{
	die "Error: have ", scalar(@readIDs_no_truth), " (of ", scalar(keys %$inferred_reads), ") reads without defined truth.\n";
}

validation::readLevelComparison($extendedMaster, $truth_reads_href, $truth_mappingDatabase_reads, $inferred_reads, $methodName, $n_reads_correct_byVariety->{$varietyName}, $n_reads_correct_byVariety_byLevel->{$varietyName});
validation::distributionLevelComparison($extendedMaster, $truth_mappingDatabase_distribution, $inferred_distribution, $methodName, $freq_byVariety_byLevel->{$varietyName});

				

