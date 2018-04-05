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

my @resultsSets = (
	# [
		# 'Nanopore',
		# 'tmp/hmp-nanopore_miniSeq+H',
		# '/scratch/tmp/hmp_nanopore_set7_combined_kraken_results',
		# 'tmp/truthHMP7_bwa_nanopore',
		# '/scratch/tmp/hmp-nanopore.fasta.fastq'
	# ],
	[
		'PacBio',
		'tmp/hmp7_2_miniSeq+H',
		'/scratch/tmp/hmp_set7_combined_kraken_results',
		'tmp/truthHMP7_bwa_pacbio',
		'/scratch/tmp/hmp_set7_combined.fastq'
	],
);

foreach my $resultsSet (@resultsSets)
{
	my $MetaMap_results = $resultsSet->[1];
	my $kraken_results_dir = $resultsSet->[2];
	my $truth = $resultsSet->[3];
	my $fastq = $resultsSet->[4];

	my %resultsFiles = (
		'MetaMap-EM' => [$MetaMap_results . '.EM.reads2Taxon',  $MetaMap_results . '.EM.WIMP'],
		'MetaMap-U' => [$MetaMap_results . '.U.reads2Taxon', $MetaMap_results . '.U.WIMP'],
		#'Kraken' => [$kraken_results_dir . '/results_kraken.txt.reads2Taxon', $kraken_results_dir . '/results_kraken.txt.ignoreUnclassified'],
		#'Bracken' => [undef, $kraken_results_dir . '/results_bracken.txt.ignoreUnclassified'],
		'Kraken' => [$kraken_results_dir . '/results_kraken.txt.reads2Taxon', $kraken_results_dir . '/results_kraken.txt'],
		'Bracken' => [undef, $kraken_results_dir . '/results_bracken.txt'],
	);

	my $missingFiles = 0;
	foreach my $method (keys %resultsFiles)
	{	
		foreach my $file (@{$resultsFiles{$method}})
		{
			next unless(defined $file);
			unless(-e $file)
			{
				warn "Missing file for $method: $file";
				$missingFiles++;
			}	
		}
	}		
	die if($missingFiles);
		
	my $master_taxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);

	my $taxonomy_usedForInference = taxTree::readTaxonomy($DB . '/taxonomy');
	(my $extendedMaster, my $extendedMaster_merged) = validation::prepare_masterTaxonomy_withX($masterTaxonomy_dir, $taxonomy_usedForInference);

	my $truth_reads_href = validation::readTruthFileReads($extendedMaster, $extendedMaster_merged, $truth . '.perRead');
	my $truth_reads_href_noUnknown = { map {$_ => $truth_reads_href->{$_}} grep {$truth_reads_href->{$_}} keys %$truth_reads_href };
	
	my $readLengths_href = Util::getReadLengths($fastq, 1);
	
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
		$truth_reads_mappable = validation::translateReadsTruthToReducedTaxonomy($extendedMaster, $mappableTaxonomy, $truth_reads_href_noUnknown);
	}
	
	my $reads_below_2000 = 0;
	my %changedTaxonIDMap;
	foreach my $readID (keys %$truth_reads_href_noUnknown)
	{
		my $original_taxon_ID = $truth_reads_href_noUnknown->{$readID};
		my $new_taxon_ID = $truth_reads_mappable->{$readID};
		if($original_taxon_ID ne $new_taxon_ID)
		{
			$changedTaxonIDMap{$original_taxon_ID}{$new_taxon_ID}++;
		}
		die unless(exists $readLengths_href->{$readID});
		if($readLengths_href->{$readID} < 2000)
		{
			$reads_below_2000++;
		}
	}	
	print "Set $resultsSet->[0]\n";
	print "\tOf ", scalar(keys %$truth_reads_href_noUnknown), " reads, $reads_below_2000 are < 2000 in length\n";
	print "\tChanged taxon IDs:\n";
	foreach my $oldTaxon (sort keys %changedTaxonIDMap)
	{
		foreach my $newTaxon (sort keys %{$changedTaxonIDMap{$oldTaxon}})
		{
			print "\t\t", $oldTaxon, '[', ($changedTaxonIDMap{$oldTaxon} ? 1 : 0), '] -> ', $newTaxon, '[', ($changedTaxonIDMap{$newTaxon} ? 1 : 0), ']: ', $changedTaxonIDMap{$oldTaxon}{$newTaxon}, " reads.\n";
		}
	}
	# die unless(all {exists $mappableTaxonomy->{$_}} values %$truth_reads_href_noUnknown);

	my @methodNames;
	my @inferred_reads;
	my @inferred_distributions; 
	foreach my $method (keys %resultsFiles)
	{
		push(@methodNames, $method);
		
		my $fn_inference_reads = $resultsFiles{$method}[0];
		my $fn_inference_distribution = $resultsFiles{$method}[1];
		
		if(defined $fn_inference_reads)
		{
			my $inferred_reads = validation::readInferredFileReads($extendedMaster, $extendedMaster_merged, $fn_inference_reads);
			my @keys_with_defined_truth = grep {exists $truth_reads_href_noUnknown->{$_}} keys %$inferred_reads;
			$inferred_reads = {map {$_ => $inferred_reads->{$_}} @keys_with_defined_truth};
			push(@inferred_reads, $inferred_reads);
		}
		else
		{
			push(@inferred_reads, undef);
		}
		
		if(defined $fn_inference_distribution)
		{
			my $inferred_distribution = validation::readInferredDistribution($extendedMaster, $extendedMaster_merged, $fn_inference_distribution);	
			push(@inferred_distributions, $inferred_distribution);
		}
		else
		{
			push(@inferred_distributions, undef);
		}		
	}
	
	
	my $allSimulations_data_href = validation::analyseAndAddOneExperiment(
		$extendedMaster,
		\%reduced_taxonID_master_2_contigs,
		$truth_reads_href_noUnknown,
		$readLengths_href, 
		\@methodNames,
		\@inferred_reads,
		\@inferred_distributions,
		'fullDB',
	);	

	$allSimulations_data_href->{realizedN} = 1;
	
	validation::produceValidationOutputFiles($allSimulations_data_href, $extendedMaster, 'HMPresults/' . $resultsSet->[0], 'HMPresults/', $resultsSet->[0]);

	
	# my $n_reads_correct_byVariety = {};
	# my $n_reads_correct_byVariety_byLevel = {};
	# my $freq_byVariety_byLevel = {};
	# my $n_reads_correct_byVariety_byLevel_byLength = {};
		
	# foreach my $label (keys %results_readLevel)
	# {
		# print "Analysing $label -- at level of individual reads!\n";
		
		# $n_reads_correct_byVariety->{$label} = {} unless(defined $n_reads_correct_byVariety->{$label});
		# $n_reads_correct_byVariety_byLevel->{$label} = {} unless(defined $n_reads_correct_byVariety_byLevel->{$label});
		# $n_reads_correct_byVariety_byLevel_byLength->{$label} = {} unless(defined $n_reads_correct_byVariety_byLevel_byLength->{$label});
		

		# my $inferred_reads = validation::readInferredFileReads($extendedMaster, $extendedMaster_merged, $results_readLevel{$label});

		# my @readIDs_no_truth = grep {not exists $truth_reads_href->{$_}} keys %$inferred_reads;
		# if(scalar(@readIDs_no_truth))
		# {
			# die "Error: have ", scalar(@readIDs_no_truth), " (of ", scalar(keys %$inferred_reads), ") reads without defined truth.\n";
		# }
		
		# validation::readLevelComparison($extendedMaster, $truth_reads_href, $truth_reads_mappable, $inferred_reads, $label, $n_reads_correct_byVariety->{$label}, $n_reads_correct_byVariety_byLevel->{$label}, $n_reads_correct_byVariety_byLevel_byLength->{$label}, \%reduced_taxonID_master_2_contigs, $readLengths_href);
	# }


		
		
		
	my @evaluateAccuracyAtLevels = validation::getEvaluationLevels();

	# print "A\n";
	my $truth_mappingDatabase_distribution = validation::truthReadsToTruthSummary($mappableTaxonomy, $truth_reads_mappable, \%reduced_taxonID_master_2_contigs);
	# print "B\n";
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

	my $freq_byVariety_byLevel = {};
	for(my $methodI = 0; $methodI <= $#methodNames; $methodI++)
	{
		my $methodName = $methodNames[$methodI];
		my $inferred_distribution = $inferred_distributions[$methodI];
		
		print "Analysing $methodName -- at level of distribution!\n";

		$freq_byVariety_byLevel->{$methodName} = {} unless(defined $freq_byVariety_byLevel->{$methodName});

		foreach my $level (@evaluateAccuracyAtLevels)
		{		
			if(defined $inferred_distribution->{$level})
			{
				foreach my $inferredTaxonID (keys %{$inferred_distribution->{$level}})  
				{
					$union_taxonIDs_byLevel{$level}{$inferredTaxonID}++;
					$distributions_byLevel_byLabel{$level}{$methodName}{$inferredTaxonID} = $inferred_distribution->{$level}{$inferredTaxonID}[1];
				}
			}
		}			
	}
	


	my $fn_output = '_HMP_distributions_' . $resultsSet->[0] . '.txt';
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
					die;
					$taxonLabel = 'Undefined';
				}
				elsif($taxonID eq 'NotLabelledAtLevel')
				{ 
					die;
					$taxonLabel = 'NotLabelledAtLevel';
				}			
				else
				{
					# print "C\n";
					$taxonLabel = taxTree::taxon_id_get_name($taxonID, $master_taxonomy);			
				}
				print F join("\t", $level, $label, $taxonID2, $taxonLabel, $distributions_byLevel_byLabel{$level}{$label}{$taxonID}), "\n";
			}
		}
	}
	
	# foreach my $label (keys %results_distribution)
	# {
		# print "Analysing $label -- at level of distribution!\n";

		# $freq_byVariety_byLevel->{$label} = {} unless(defined $freq_byVariety_byLevel->{$label});

		# my $inferred_distribution = validation::readInferredDistribution($extendedMaster, $extendedMaster_merged, $results_distribution{$label});	
 
		
		# foreach my $level (@evaluateAccuracyAtLevels)
		# {		
			# if(defined $inferred_distribution->{$level})
			# {
				# foreach my $inferredTaxonID (keys %{$inferred_distribution->{$level}})
				# {
					# $union_taxonIDs_byLevel{$level}{$inferredTaxonID}++;
					# $distributions_byLevel_byLabel{$level}{$label}{$inferredTaxonID} = $inferred_distribution->{$level}{$inferredTaxonID}[1];
				# }
			# }
		# }	
	# }	
}
				
			
		