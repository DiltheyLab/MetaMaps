#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/all any shuffle sum/;
use List::MoreUtils qw/mesh/;
use Getopt::Long;   
use File::Path qw(make_path remove_tree);
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;
use File::Copy;
use Storable qw/dclone store retrieve/;
use Math::Random;

# TODO
# REMOVE ONE SPECIES REMOVALL!!!
# (currently we're always removing the same base node)

# todo real random seed?
srand(12345);
$| = 1;

use taxTree;
use Util;
use SimulationsKraken;
use SimulationsMetaPalette;
use validation;
use simulation;

Util::get_metaMap_bin_and_enforce_mainDir();
die unless(-e 'estimateSelfSimilarity.pl');

#my $create_simulations = 2;
#my $perSimulation_totalKnownGenomes = 5;

my $simulation_read_length = 5000;

# my $taxonomyDir = '/data/projects/phillippy/projects/mashsim/NCBI/refseq/taxonomy/';
my $masterTaxonomy_dir = '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';

# these are populated globally if required
my $masterTaxonomy; 
my $masterTaxonomy_merged;

#my $ref = '../db/refseq/ref.fa';
#my $contigsTaxa = $ref . '.taxa';

#my $treeSelfSimiliarityReads = '/data/projects/phillippy/projects/mashsim/NCBI/refseq/taxonomy/selfSimilarity/results.reads.byNode';
#(my $treeSelfSimilarityReadsMany = $treeSelfSimiliarityReads) =~ s/results\.reads\.byNode/results.reads.many.byNode/;


my $PBsim_cmd = qq(/data/projects/phillippy/projects/mashsim/PBSIM-PacBio-Simulator/src/pbsim --model_qc /data/projects/phillippy/projects/mashsim/PBSIM-PacBio-Simulator/data/model_qc_clr --data-type CLR --depth DEPTH --prefix --length-mean $simulation_read_length --accuracy-mean 0.88 REF);

my $metamap_bin = './metamap';
my $metaPalette_installation_dir = qq(/data/projects/phillippy/software/MetaPalette/);
my $jellyfish_2_bin = qq(/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish);

my $kraken_binPrefix = SimulationsKraken::getKrakenBinPrefix();
my $Bracken_dir = SimulationsKraken::getBrackenDir();
my $krakenDBTemplate = SimulationsKraken::getKrakenDBTemplate();

die unless(-e $metamap_bin);

# my $taxonomy_full = taxTree::readTaxonomy($fullReferenceTaxonomy);
# my $taxonomy_kraken = taxTree::readTaxonomy($krakenDBTemplate . '/taxonomy');

# my %combined_nodes = map {$_ => 1} ((keys %$taxonomy_full), (keys %$taxonomy_kraken));
# my $node_in_both = 0;
# my $node_in_full = 0;
# my $node_in_kraken = 0;
# foreach my $node (keys %combined_nodes)
# {
	# if((exists $taxonomy_full->{$node}) and (exists $taxonomy_kraken->{$node}))
	# {
		# $node_in_both++;
	# }
	# elsif((not exists $taxonomy_full->{$node}) and (exists $taxonomy_kraken->{$node}))
	# {
		# $node_in_kraken++;
	# }
	# elsif((exists $taxonomy_full->{$node}) and (not exists $taxonomy_kraken->{$node}))
	# {	
		# $node_in_full++;
	# }
# }
# print "Kraken / reference taxonomy comparison:\n";
# print "\tNodes in both       : $node_in_both \n";
# print "\tNodes in ref only   : $node_in_full \n";
# print "\tNodes in Kraken only: $node_in_kraken \n";

# exit;

# my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

# test Kraken/Bracken file conversion

# create_compatible_file_from_kraken(
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/results_kraken_manual.txt',
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/testDB/taxonomy',
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/combined_reads.fastq.classified.report',
# );
# exit;

# create_compatible_file_from_kraken_bracken(
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/results_bracken_manual.txt',
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/testDB/taxonomy',
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/combined_reads.fastq.classified.report',
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/combined_reads.fastq.classified.bracken_S',
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/combined_reads.fastq.classified.bracken_G',
	# '/data/projects/phillippy/projects/mashsim/src/simulations/0/combined_reads.fastq.classified.bracken_F'
# );
# exit; 

# create_compatible_file_from_metapalette( 
	# '/data/projects/phillippy/software/MetaPalette/simReads/combined_reads.fastq.profile.compatible',
	# '/data/projects/phillippy/software/MetaPalette/simReads/combined_reads.fastq.profile'
# );
# exit;


system("export PATH=/data/projects/phillippy/software/jellyfish-1.1.11/bin:\$PATH") and die "Couldn't add Jellyfish to PATH";

#system(". /etc/profile.d/modules.sh") and die;
#system('MODULEPATH=$MODULEPATH:/data/projects/phillippy/software/modules/') and die;
system("eval 'module load libs/hdf5/1.8.12'") and die "Could load hdf5";
system("eval 'module load libs/scipy/0.18.1'") and die "Could load scipy";
	
my $action;
my $DB = 'databases/miniSeq_100';
my $jobI; 
my $really;
my $jobIMethod = 'all';
my $n_simulations = 10;
my $useVarietyI;
my $suffix;
my $n_species = 10;
my $coverageMode = "equal";
my $desiredTaxa;
GetOptions (
	'DB:s' => \$DB,
	'action:s' => \$action,
	'jobI:s' => \$jobI,
	'really:s' => \$really,
	'n_simulations:s' => \$n_simulations,
	'jobIMethod:s' => \$jobIMethod,
	'varietyI:s' => \$useVarietyI,
	'suffix:s' => \$suffix,
	'n_species:s' => \$n_species,
	'desiredTaxa:s' => \$desiredTaxa,
	'coverageMode:s' => \$coverageMode,
);

die unless(($coverageMode eq 'equal') or ($coverageMode eq 'logNormal') or ($coverageMode eq 'file'));

if(not defined $action)
{
	die "Please specify --action, e.g. prepare, prepareII, analyzeAll, inferenceJobI, analyzeJobI";
}

unless(defined $DB)
{
	die "Please specify --DB, e.g. databases/miniSeq";
}

my $globalOutputDir = $DB . '/simulations';
if($suffix)
{
	$globalOutputDir .= ('_' . $suffix);
}

unless(-e $globalOutputDir)
{
	mkdir($globalOutputDir) or die "Cannot mkdir $globalOutputDir";
}

if($action eq 'prepare')
{
	my $DB_fa = $DB . '/DB.fa';
	my $existingSelfSimilarities = $DB . '/selfSimilarities.txt';
	
	die "Simulations in directory $globalOutputDir already prepared" if getFlag($globalOutputDir, 'simulated');

	die unless(-e $DB_fa);
	die unless(-e $existingSelfSimilarities);
	
	my $MetaMap_taxonomy = {};
	my %contigID_2_taxonID;
	my %taxonID_2_contigIDs;	
	my %contig_2_length;
	fill_contigID_taxonID($DB, $MetaMap_taxonomy, \%contigID_2_taxonID, \%taxonID_2_contigIDs, \%contig_2_length);
	my %taxon_2_genomeLength = map {$_ => getGenomeLength($_, \%taxonID_2_contigIDs, \%contig_2_length)} keys %taxonID_2_contigIDs;
	
	my $fn_log = $globalOutputDir. '/log.txt';	
	my $fn_selfSimilarities = $globalOutputDir. '/selfSimilarities.txt';	
	my $fn_selfSimilarities_collect = $globalOutputDir. '/selfSimilarities_collect.txt';	
	my $fn_selfSimilarities_fromTemplate = $globalOutputDir. '/selfSimilarities_fromTemplate.txt';	
	my $fn_selfSimilarities_fromTemplate_asArray = $globalOutputDir. '/selfSimilarities_fromTemplate.txt.asArray';	

	my $fh_log;
	open($fh_log, '>', $fn_log) or die "Cannot open $fn_log";

	my $fh_selfSimilarities;
	open($fh_selfSimilarities, '>', $fn_selfSimilarities) or die "Cannot open $fn_selfSimilarities";

	my $fh_selfSimilarities_collect;
	open($fh_selfSimilarities_collect, '>', $fn_selfSimilarities_collect) or die "Cannot open $fn_selfSimilarities_collect";
		
	my $fh_selfSimilarities_prepareFromOthers;
	open($fh_selfSimilarities_prepareFromOthers, '>', $fn_selfSimilarities_fromTemplate) or die "Cannot open $fn_selfSimilarities_fromTemplate";
		
	my @simulations_to_execute;
	for(my $simulationI = 0; $simulationI < $n_simulations; $simulationI++)
	{
		print {$fh_log} "Outer simulation $simulationI\n";
		
		my $thisSimulation_outputDirectory = $globalOutputDir . '/' . $simulationI;
		(mkdir($thisSimulation_outputDirectory) or die "Cannot mkdir $thisSimulation_outputDirectory") unless(-e $thisSimulation_outputDirectory);
		
		my $oneSimulation_href;
		if($simulationI < $n_simulations)
		{
			die unless(($coverageMode eq 'equal') or ($coverageMode eq 'logNormal'));
			my $equalCoverage = ($coverageMode eq 'equal');
			$oneSimulation_href = returnOneSimulation_noUnknown($DB, \%taxonID_2_contigIDs, \%taxon_2_genomeLength, $n_species, $thisSimulation_outputDirectory, $equalCoverage, $MetaMap_taxonomy, $fh_log);
		}
		
		addInferenceRoundsWithReducedDBs($oneSimulation_href, $MetaMap_taxonomy, $fh_log);
		
		push(@simulations_to_execute, $oneSimulation_href);
	}

	for(my $simulationI = 0; $simulationI <= $#simulations_to_execute; $simulationI++)
	{
		print {$fh_log} "Executing outer simulation $simulationI\n";
		my $oneSimulation_href = $simulations_to_execute[$simulationI];
		executeSimulation($oneSimulation_href, \%taxonID_2_contigIDs, \%contig_2_length, $MetaMap_taxonomy, $fh_log, $fh_selfSimilarities, $fh_selfSimilarities_collect, $fh_selfSimilarities_prepareFromOthers);
	}
	
	my $simulations_qsub_file = $globalOutputDir. '/qsub.txt';	 
	open(QSUB, '>', $simulations_qsub_file) or die "Cannot open $simulations_qsub_file";
print QSUB qq(#!/bin/bash
#\$ -t 1-${n_simulations}
#\$ -q public.q
#\$ -l mem_free=172G
jobID=\$(expr \$SGE_TASK_ID - 1)
cd $FindBin::Bin
perl simulate.pl --action inferenceJobI --DB $DB --suffix $suffix --jobI \$jobID
);	
	close(QSUB);
	close($fh_log);
	close($fh_selfSimilarities);
	close($fh_selfSimilarities_collect);
	close($fh_selfSimilarities_prepareFromOthers);
	
	my $n_simulations_file = $globalOutputDir . '/n_simulations.txt';
	open(N, '>', $n_simulations_file) or die "Cannot open $n_simulations_file";
	print N $n_simulations, "\n";
	close(N); 
	
	createArrayJob($fn_selfSimilarities_fromTemplate, $fn_selfSimilarities_fromTemplate_asArray);
	
	# print "$n_simulations prepared on DB $DB -- if you're in an SGE environment, you can now\n\n\texecute the commands in $fn_selfSimilarities\n\tand later these in $fn_selfSimilarities_collect\n\n\tqsub $simulations_qsub_file\n\n";
	print "$n_simulations prepared on DB $DB -- if you're in an SGE environment, you can now\n\n\tqsub $fn_selfSimilarities_fromTemplate_asArray\n\tand when all jobs are done call me with --mode prepareII\n\n";
	
	setFlag($globalOutputDir, 'simulated', 1);
}
elsif($action eq 'prepareFromFile')
{
	my $DB_fa = $DB . '/DB.fa';
	my $existingSelfSimilarities = $DB . '/selfSimilarities.txt';
	
	die "Simulations in directory $globalOutputDir already prepared" if getFlag($globalOutputDir, 'simulated');

	die unless(-e $DB_fa);
	die unless(-e $existingSelfSimilarities);
	die "Please set --n_simulations to 1 if doing --mode prepareFromFile" unless($n_simulations == 1);
	
	my $MetaMap_taxonomy = {};
	my %contigID_2_taxonID;
	my %taxonID_2_contigIDs;	
	my %contig_2_length;
	fill_contigID_taxonID($DB, $MetaMap_taxonomy, \%contigID_2_taxonID, \%taxonID_2_contigIDs, \%contig_2_length);
	my %taxon_2_genomeLength = map {$_ => getGenomeLength($_, \%taxonID_2_contigIDs, \%contig_2_length)} keys %taxonID_2_contigIDs;

	my $fn_log = $globalOutputDir. '/log.txt';	
	my $fn_selfSimilarities = $globalOutputDir. '/selfSimilarities.txt';	
	my $fn_selfSimilarities_collect = $globalOutputDir. '/selfSimilarities_collect.txt';	
	my $fn_selfSimilarities_fromTemplate = $globalOutputDir. '/selfSimilarities_fromTemplate.txt';	
	my $fn_selfSimilarities_fromTemplate_asArray = $globalOutputDir. '/selfSimilarities_fromTemplate.txt.asArray';	

	my $fh_log;
	open($fh_log, '>', $fn_log) or die "Cannot open $fn_log";

	my $fh_selfSimilarities;
	open($fh_selfSimilarities, '>', $fn_selfSimilarities) or die "Cannot open $fn_selfSimilarities";

	my $fh_selfSimilarities_collect;
	open($fh_selfSimilarities_collect, '>', $fn_selfSimilarities_collect) or die "Cannot open $fn_selfSimilarities_collect";
		
	my $fh_selfSimilarities_prepareFromOthers;
	open($fh_selfSimilarities_prepareFromOthers, '>', $fn_selfSimilarities_fromTemplate) or die "Cannot open $fn_selfSimilarities_fromTemplate";

	my @simulations_to_execute;
	{
		print {$fh_log} "Outer simulation from file $desiredTaxa\n";
		
		my $thisSimulation_outputDirectory = $globalOutputDir . '/' . 0;
		(mkdir($thisSimulation_outputDirectory) or die "Cannot mkdir $thisSimulation_outputDirectory") unless(-e $thisSimulation_outputDirectory);
		
		my $oneSimulation_href = returnOneSimulation_fromFile($DB, \%taxonID_2_contigIDs, \%taxon_2_genomeLength, $thisSimulation_outputDirectory, $coverageMode, $MetaMap_taxonomy, $fh_log, $desiredTaxa);
		
		if($coverageMode ne 'file')
		{
			addInferenceRoundsWithReducedDBs_allSpecies($oneSimulation_href, $MetaMap_taxonomy, $fh_log);
		}
		
		push(@simulations_to_execute, $oneSimulation_href);
	}

	for(my $simulationI = 0; $simulationI <= $#simulations_to_execute; $simulationI++)
	{
		print {$fh_log} "Executing outer simulation $simulationI\n";
		my $oneSimulation_href = $simulations_to_execute[$simulationI];
		executeSimulation($oneSimulation_href, \%taxonID_2_contigIDs, \%contig_2_length, $MetaMap_taxonomy, $fh_log, $fh_selfSimilarities, $fh_selfSimilarities_collect, $fh_selfSimilarities_prepareFromOthers);
	}
	
	# my $simulations_qsub_file = $globalOutputDir. '/qsub.txt';	 
	# open(QSUB, '>', $simulations_qsub_file) or die "Cannot open $simulations_qsub_file";
# print QSUB qq(#!/bin/bash
# #\$ -t 1-${n_simulations}
# #\$ -q public.q
# #\$ -l mem_free=172G
# jobID=\$(expr \$SGE_TASK_ID - 1)
# cd $FindBin::Bin
# perl simulate.pl --action inferenceJobI --DB $DB --jobI \$jobID
# );	
	# close(QSUB);
	
	close($fh_log);
	close($fh_selfSimilarities);
	close($fh_selfSimilarities_collect);
	close($fh_selfSimilarities_prepareFromOthers);
	
	my $n_simulations_file = $globalOutputDir . '/n_simulations.txt';
	open(N, '>', $n_simulations_file) or die "Cannot open $n_simulations_file";
	print N $n_simulations, "\n";
	close(N); 
	
	my $n_arrayJobs = createArrayJob($fn_selfSimilarities_fromTemplate, $fn_selfSimilarities_fromTemplate_asArray);
	
	# print "$n_simulations prepared on DB $DB -- if you're in an SGE environment, you can now\n\n\texecute the commands in $fn_selfSimilarities\n\tand later these in $fn_selfSimilarities_collect\n\n\tqsub $simulations_qsub_file\n\n";
	if($n_arrayJobs)
	{
		print "$n_simulations prepared on DB $DB -- if you're in an SGE environment, you can now\n\n\tqsub $fn_selfSimilarities_fromTemplate_asArray\n\tand when all jobs are done call me with --mode prepareII\n\n";
	}
	else
	{
		print "$n_simulations prepared on DB $DB -- call me with --mode prepareII\n\n";
	}
	
	setFlag($globalOutputDir, 'simulated', 1);
}
elsif($action eq 'prepareII')
{
	die "Simulations in $globalOutputDir not prepared yet (--mode prepare)" unless getFlag($globalOutputDir, 'simulated');

	my $n_simulations_file = $globalOutputDir . '/n_simulations.txt';
	open(N, '<', $n_simulations_file) or die "Cannot open $n_simulations_file";
	my $realizedN = <N>;
	chomp($realizedN);
	close(N);
	die unless($realizedN =~ /^\d+$/);
	
	my %n_reads_correct_byVariety;
	my %n_reads_correct_byVariety_byLevel;
	my %freq_byVariety_byLevel;
	
	# for(my $jobI = 0; $jobI < $realizedN; $jobI++)
	my @commands;
	for(my $jobI = 0; $jobI < $realizedN; $jobI++)  
	{
		my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
		my $simulation_href = retrieve $simulation_href_fn;
					
		for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++)
		{			
			my $DB_target_dir = $simulation_href->{outputDirectory} . '/DB_' . $simulation_href->{inferenceDBs}[$varietyI][2];
			die unless($DB_target_dir eq $simulation_href->{dbDirs_metamap}[$varietyI]);
		
			my $fn_selfSimilarities = $DB_target_dir . '/selfSimilarities.txt';
			die "File $fn_selfSimilarities missing" unless(-e $fn_selfSimilarities);
			
			my $cmd = qq(perl simulate.pl --action inferenceJobI --DB $DB --jobI $jobI --suffix $suffix --varietyI $varietyI);
			
			push(@commands, $cmd);
		}	
	}	
	
	my $fn_jobs_nonArray = $globalOutputDir . '/inferenceJobs.noArray';
	open(F, '>', $fn_jobs_nonArray) or die;
	print F join("\n", @commands), "\n";
	close(F);
	
	my $fn_jobs_asArray = $globalOutputDir . '/inferenceJobs.asArray';

	createArrayJob($fn_jobs_nonArray, $fn_jobs_asArray, 172, 'public.q');

	print "\n", scalar(@commands), " individual jobs, now qsub $fn_jobs_asArray \n\n";
}
elsif($action eq 'inferenceJobI')
{
	die"Please specify --jobI" unless(defined $jobI);
	my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
	my $simulation_href = retrieve $simulation_href_fn;
	inferenceOneSimulation($simulation_href);
}
elsif($action eq 'analyzeJobI')
{
	die"Please specify --jobI" unless(defined $jobI);
	my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
	my $simulation_href = retrieve $simulation_href_fn;
	evaluateOneSimulation($simulation_href);
}
elsif($action eq 'analyzeAll')
{
	my $n_simulations_file = $globalOutputDir . '/n_simulations.txt';
	open(N, '<', $n_simulations_file) or die "Cannot open $n_simulations_file";
	my $realizedN = <N>;
	chomp($realizedN);
	close(N);
	die unless($realizedN =~ /^\d+$/);

	my @n_reads_correct_byVariety_bySimulation;
	my @n_reads_correct_byVariety_byLevel_bySimulation;;
	my @freq_byVariety_byLevel_bySimulation;;
	
	my %n_reads_correct_byVariety;
	my %n_reads_correct_byVariety_byLevel;
	my %freq_byVariety_byLevel;
	my @frequencyComparisons_bySimulation;
	
	my $fullTaxonomy_simulation = taxTree::readTaxonomy($DB . '/taxonomy');
	
	my @highLevel_stats_keptSeparate_bySimulation;
	my %callRate_and_accuracy_byReadCategory;
	# for(my $jobI = 0; $jobI < $realizedN; $jobI++)
	for(my $jobI = 0; $jobI < $realizedN; $jobI++)  
	{
		my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
		my $simulation_href = retrieve $simulation_href_fn;
		my $frequencyComparison = {};
		my %n_reads_correct_byVariety_local;
		my %n_reads_correct_byVariety_byLevel_local;		
		my %freq_byVariety_byLevel_local;		
		evaluateOneSimulation($simulation_href, \%n_reads_correct_byVariety_local, \%n_reads_correct_byVariety_byLevel_local, \%freq_byVariety_byLevel_local, $frequencyComparison);
		push(@frequencyComparisons_bySimulation, $frequencyComparison);

		# variety = fullDB/removeOne_genus ...
		# label = MetaMap / Kraken ...
		# category = read category ...
		foreach my $variety (keys %n_reads_correct_byVariety_local)
		{		
			foreach my $label (keys %{$n_reads_correct_byVariety_local{$variety}})
			{
				foreach my $category (keys %{$n_reads_correct_byVariety_local{$variety}{$label}})
				{
					foreach my $key (keys %{$n_reads_correct_byVariety_local{$variety}{$label}{$category}})
					{
						my $value = $n_reads_correct_byVariety_local{$variety}{$label}{$category}{$key};
						die unless(not ref($value));
						$n_reads_correct_byVariety{$variety}{$label}{$category}{$key} += $value;
					}
					
					my $d = $n_reads_correct_byVariety_local{$variety}{$label}{$category};
					die unless(exists $d->{N});
					die unless(exists $d->{missing});
					die unless(exists $d->{correct});
					my $N = $d->{N} + $d->{missing};
					die unless($N > 0);
					die unless($d->{N} > 0);
				
					my $CR = $d->{N} / $N; die unless(($CR >= 0) and ($CR <= 1));
					my $accuracy = $d->{correct} / $d->{N}; die unless(($accuracy >= 0) and ($accuracy <= 1));
					
					$highLevel_stats_keptSeparate_bySimulation[$jobI]{$variety}{$label}{$category}{mappingTarget}{CR} = $CR;
					$highLevel_stats_keptSeparate_bySimulation[$jobI]{$variety}{$label}{$category}{mappingTarget}{Accuracy} = $accuracy;
					push(@{$callRate_and_accuracy_byReadCategory{$category}{$label}{mappingTarget}}, [$CR, $accuracy]);
				}
			}
		}
		

		foreach my $variety (keys %n_reads_correct_byVariety_byLevel_local)
		{		
			foreach my $label (keys %{$n_reads_correct_byVariety_byLevel_local{$variety}})
			{		
				foreach my $category (keys %{$n_reads_correct_byVariety_byLevel_local{$variety}{$label}})
				{
					foreach my $level (keys %{$n_reads_correct_byVariety_byLevel_local{$variety}{$label}{$category}})
					{
						foreach my $key (keys %{$n_reads_correct_byVariety_byLevel_local{$variety}{$label}{$category}{$level}})
						{
							my $value = $n_reads_correct_byVariety_byLevel_local{$variety}{$label}{$category}{$level}{$key};
							die unless(not ref($value));
							$n_reads_correct_byVariety_byLevel{$variety}{$label}{$category}{$level}{$key} += $value;
						}
						

						my $d = $n_reads_correct_byVariety_byLevel_local{$variety}{$label}{$category}{$level};
						die unless(exists $d->{N});
						die unless(exists $d->{missing});
						die unless(exists $d->{correct});
						my $N = $d->{N} + $d->{missing};
						die unless($N > 0);
						die unless($d->{N} > 0);
						
						my $CR = $d->{N} / $N; die unless(($CR >= 0) and ($CR <= 1));
						my $accuracy = $d->{correct} / $d->{N}; die unless(($accuracy >= 0) and ($accuracy <= 1));
						
						$highLevel_stats_keptSeparate_bySimulation[$jobI]{$variety}{$label}{$category}{$level}{CR} = $CR;
						$highLevel_stats_keptSeparate_bySimulation[$jobI]{$variety}{$label}{$category}{$level}{Accuracy} = $accuracy;						
						push(@{$callRate_and_accuracy_byReadCategory{$category}{$label}{$level}}, [$CR, $accuracy]);						
					}
				}
			}
		}
		
		foreach my $variety (keys %freq_byVariety_byLevel_local)
		{			
			foreach my $label (keys %{$freq_byVariety_byLevel_local{$variety}})
			{		
				foreach my $level (keys %{$freq_byVariety_byLevel_local{$variety}{$label}})
				{
					foreach my $key (keys %{$freq_byVariety_byLevel_local{$variety}{$label}{$level}})
					{
						my $value = $freq_byVariety_byLevel_local{$variety}{$label}{$level}{$key};
						die unless(defined $value);
						die unless((not ref($value)) or (ref($value) eq 'ARRAY'));
						if(not ref($value))
						{
							$freq_byVariety_byLevel{$variety}{$label}{$level}{$key} += $value;
						}
						else
						{
							push(@{$freq_byVariety_byLevel{$variety}{$label}{$level}{$key}}, @$value);
						}
					}
				}
			}
		}				

		push(@n_reads_correct_byVariety_bySimulation, \%n_reads_correct_byVariety_local);
		push(@n_reads_correct_byVariety_byLevel_bySimulation, \%n_reads_correct_byVariety_byLevel_local);
		push(@freq_byVariety_byLevel_bySimulation, \%freq_byVariety_byLevel_local);		
	}
	

	
	my @varieties = qw/fullDB removeOne_self removeOne_species removeOne_genus/;
	@varieties = grep {exists $n_reads_correct_byVariety{$_}} @varieties;
	#die Dumper(\@varieties, \%n_reads_correct_byVariety) unless(scalar(@varieties) == scalar(keys %n_reads_correct_byVariety));
	die Dumper(\@varieties, [keys %n_reads_correct_byVariety], "Issue I") unless(all {exists $n_reads_correct_byVariety{$_}} @varieties);
	
	my @levels_ordered = validation::getEvaluationLevels();
	my %level_to_i;
	for(my $levelI = 0; $levelI <= $#levels_ordered; $levelI++)
	{
		$level_to_i{$levels_ordered[$levelI]} = $levelI;
	}
	
	{
		open(BARPLOTSREADCAT, '>', '_forPlot_barplots_readCategory') or die;
		print BARPLOTSREADCAT join("\t", qw/readCategory evaluationLevel method callRateAvg accuracyAvg callRate_raw accuracy_raw/), "\n";
		
		foreach my $readCategory (keys %callRate_and_accuracy_byReadCategory)
		{
			foreach my $label (keys %{$callRate_and_accuracy_byReadCategory{$readCategory}})
			{
				foreach my $level (keys %{$callRate_and_accuracy_byReadCategory{$readCategory}{$label}})
				{
					my $v = $callRate_and_accuracy_byReadCategory{$readCategory}{$label}{$level};
					my @callRates;
					my @accuracies;
					foreach my $e (@$v)
					{
						push(@callRates, $e->[0]);
						push(@accuracies, $e->[1]);
					}	
					die unless(scalar(@callRates));
					die unless(scalar(@accuracies));
					my $avg_callRate = Util::mean(@callRates);
					my $avg_accuracy = Util::mean(@accuracies);
					print BARPLOTSREADCAT join("\t", $readCategory, $level, $label, $avg_callRate, $avg_accuracy, join(';', @callRates), join(';', @accuracies)), "\n";
				}
			}
		}
		close(BARPLOTSREADCAT);
				
	}

	{
		my %_methods;
		my %_readStratification;
		my %_evaluationLevels;
		foreach my $variety (@varieties)
		{
			foreach my $methodName (keys %{$n_reads_correct_byVariety{$variety}})
			{
				$_methods{$methodName}++;
				foreach my $category (keys %{$n_reads_correct_byVariety{$variety}{$methodName}})
				{
					$_readStratification{$category}++;
					
					if(exists $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category})
					{
						foreach my $k (keys %{$n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category}})
						{
							$_evaluationLevels{$k}++;
						}
					}
				}
			}
		}
		die unless(all {$_methods{$_} == scalar(@varieties)} keys %_methods);
		my @methods = sort keys %_methods;
		my @readLevels = sort {
			if(($a =~ /novel_to/) and ($b =~ /novel_to/))
			{
				die unless($a =~ /novel_to_(.+)/);
				my $a_level = $1;
				die unless($b =~ /novel_to_(.+)/);
				my $b_level = $1;
				die "Undefined level $a_level" unless(defined $level_to_i{$a_level});				
				die "Undefined level $b_level" unless(defined $level_to_i{$b_level});				
				$level_to_i{$a_level} <=> $level_to_i{$b_level}
			}
			else
			{
				$a cmp $b
			}
		} keys %_readStratification;
		# my @evaluationLevels = sort keys %_evaluationLevels;
		my @evaluationLevels = qw/species genus family/;
		die Dumper("Missing evaluation levels", \@evaluationLevels, \%_evaluationLevels) unless(all {exists $_evaluationLevels{$_}} @evaluationLevels);
		
		open(BARPLOTSFULLDB, '>', '_forPlot_barplots_fullDB') or die;
		print BARPLOTSFULLDB join("\t", qw/readLevel variety method level callRate accuracy/), "\n";
			
		{
			open(READSABSOLUTELYCORRECT, '>', '_readsAbsolutelyCorrect') or die;
			my @header_fields_1_absolutelyCorrect = ('ReadLevel');
			my @header_fields_2_absolutelyCorrect = ('');
			my @header_fields_3_absolutelyCorrect = ('');

			foreach my $variety (@varieties)
			{
				my $hf2_before = $#header_fields_2_absolutelyCorrect;
				foreach my $method (@methods)
				{
					push(@header_fields_2_absolutelyCorrect, $method, '', '', '', '');		
					push(@header_fields_3_absolutelyCorrect, 'Ntotal', 'OKtotal', 'NmadeCall', 'OKmadeCall', 'noCall');
				}
				my $hf2_after = $#header_fields_2_absolutelyCorrect;
				my $requiredFields = $hf2_after - $hf2_before;
				die unless($requiredFields > 0);
				my @addToHeader1 = ($variety, (('') x ($requiredFields - 1)));
				die unless(scalar(@addToHeader1) == $requiredFields);
				push(@header_fields_1_absolutelyCorrect, @addToHeader1);
			}
			
			print READSABSOLUTELYCORRECT join("\t", @header_fields_1_absolutelyCorrect), "\n";
			print READSABSOLUTELYCORRECT join("\t", @header_fields_2_absolutelyCorrect), "\n";
			print READSABSOLUTELYCORRECT join("\t", @header_fields_3_absolutelyCorrect), "\n";
			
			foreach my $readLevel (@readLevels)
			{
				my @output_fields_absolutelyCorrect = ($readLevel);
					
				foreach my $variety (@varieties)
				{
					foreach my $methodName (@methods)
					{
						my $missing = 0;
						my $NmadeCall = 0;
						my $correct = 0;
						
						my $Ntotal = 0;
						my $percOK_madeCall = 0;
						my $percOK_madeCall_fullAccuracy = 0;
						my $percOK_total = 0;
						my $perc_missing = 0;
						my $callRate = 'NA';
						
						if(exists $n_reads_correct_byVariety{$variety}{$methodName}{$readLevel})
						{
							$missing =  $n_reads_correct_byVariety{$variety}{$methodName}{$readLevel}{missing};
							$NmadeCall =  $n_reads_correct_byVariety{$variety}{$methodName}{$readLevel}{N};
							$correct =  $n_reads_correct_byVariety{$variety}{$methodName}{$readLevel}{correct};
							
							$Ntotal =  $missing + $NmadeCall;
							
						}

						$percOK_madeCall = sprintf("%.2f", ($correct / $NmadeCall)) if($NmadeCall > 0);
						$percOK_madeCall_fullAccuracy = ($correct / $NmadeCall) if($NmadeCall > 0);
						$percOK_total = sprintf("%.2f", ($correct / $Ntotal)) if($Ntotal > 0);
						$perc_missing = sprintf("%.2f", ($missing / ($Ntotal))) if($Ntotal > 0);
						
						$callRate = $NmadeCall / $Ntotal if($Ntotal > 0);
						
						push(@output_fields_absolutelyCorrect, $Ntotal, $percOK_total, $NmadeCall, $percOK_madeCall, $perc_missing);
						
						print BARPLOTSFULLDB join("\t", $readLevel, $variety, $methodName, 'mappingTarget', $callRate, $percOK_madeCall_fullAccuracy), "\n";

					}
				}
			
				print READSABSOLUTELYCORRECT join("\t", @output_fields_absolutelyCorrect), "\n";
			}	
			

			close(READSABSOLUTELYCORRECT);
		}
		
		{
			open(READSCORRECTBYLEVEL, '>', '_readsCorrectByLevel') or die;
			my @header_fields_1_byLevelCorrect = ('ReadLevel', 'EvaluationLevel');
			my @header_fields_2_byLevelCorrect = ('', '');
			my @header_fields_3_byLevelCorrect = ('', '');	
			
			foreach my $variety (@varieties)
			{
				my $hf2_before = $#header_fields_2_byLevelCorrect;
				foreach my $method (@methods)
				{
					push(@header_fields_2_byLevelCorrect, $method, '', '', '', '');							
					push(@header_fields_3_byLevelCorrect, 'Ntotal', 'OKtotal', 'NmadeCall', 'OKmadeCall', 'noCall');
				}
				my $hf2_after = $#header_fields_2_byLevelCorrect;
				my $requiredFields = $hf2_after - $hf2_before;
				die unless($requiredFields > 0);
				my @addToHeader1 = ($variety, (('') x ($requiredFields - 1)));
				die unless(scalar(@addToHeader1) == $requiredFields);
				push(@header_fields_1_byLevelCorrect, @addToHeader1);
			}
			
			print READSCORRECTBYLEVEL join("\t", @header_fields_1_byLevelCorrect), "\n";
			print READSCORRECTBYLEVEL join("\t", @header_fields_2_byLevelCorrect), "\n";
			print READSCORRECTBYLEVEL join("\t", @header_fields_3_byLevelCorrect), "\n";

			foreach my $readLevel (@readLevels)
			{
				foreach my $evaluationLevel (@evaluationLevels)
				{
					my @output_fields_byLevelCorrect = ($readLevel, $evaluationLevel);
					foreach my $variety (@varieties)
					{
						foreach my $methodName (@methods)
						{				
							my $missing = 0;
							
							my $N_total = 0;
							my $N_total_truthDefined = 0;
							
							my $N_madeCall = 0;
							my $N_madeCall_truthDefined = 0;
							
							my $correct = 0;
							my $correct_truthDefined = 0;
							my $callRate = 'NA';
							
							my $percOK_total = 0;
							my $percOK_madeCall_fullAccuracy = 0;
							my $percOK_total_truthDefined = 0;							
							my $percOK_madeCall = 0;
							my $percOK_madeCall_truthDefined = 0;
							my $percMissing = 0;
							
							if(exists $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$readLevel}{$evaluationLevel})
							{
								$N_madeCall = $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$readLevel}{$evaluationLevel}{N};
								$N_madeCall_truthDefined = $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$readLevel}{$evaluationLevel}{N_truthDefined};
								
								$correct = $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$readLevel}{$evaluationLevel}{correct};
								$correct_truthDefined = $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$readLevel}{$evaluationLevel}{correct_truthDefined};
								
								$missing = $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$readLevel}{$evaluationLevel}{missing};
								
								$N_total = $N_madeCall + $missing;
								$N_total_truthDefined = $N_madeCall_truthDefined + $missing;
							}		

							print "N $N_total \n";
							
							
							# die Dumper("Weird - N is $N, but missing is $missing?", [$readLevel, $evaluationLevel, $variety, $methodName]) unless($missing <= $N);
							
							$callRate = $N_madeCall / $N_total if($N_total > 0);
							$percOK_total = sprintf("%.2f", ($correct / $N_total)) if($N_total > 0);
							$percOK_total_truthDefined = sprintf("%.2f", ($correct_truthDefined / $N_total_truthDefined)) if($N_total_truthDefined > 0);
							
							$percOK_madeCall = sprintf("%.2f", ($correct / $N_madeCall)) if($N_madeCall > 0);
							$percOK_madeCall_fullAccuracy = ($correct / $N_madeCall) if($N_madeCall > 0);
							$percOK_madeCall_truthDefined = sprintf("%.2f", ($correct_truthDefined / $N_madeCall_truthDefined)) if($N_madeCall_truthDefined > 0);
							
							$percMissing = sprintf("%.2f", ($missing / $N_total)) if($N_total > 0);
														
							push(@header_fields_3_byLevelCorrect, 'Ntotal', 'OKtotal', 'NmadeCall', 'OKmadeCall', 'noCall');

												
							push(@output_fields_byLevelCorrect,
								($N_total ne $N_total_truthDefined) ? join(' / ', $N_total, $N_total_truthDefined) : $N_total,
								($percOK_total ne $percOK_total_truthDefined) ? join(' / ', $percOK_total, $percOK_total_truthDefined) : $percOK_total, 
								($N_madeCall ne $N_madeCall_truthDefined) ? join(' / ', $N_madeCall, $N_madeCall_truthDefined) : $N_madeCall,
								($percOK_madeCall ne $percOK_madeCall_truthDefined) ? join(' / ', $percOK_madeCall, $percOK_madeCall_truthDefined) : $percOK_madeCall,
								$percMissing
							);
							
							print BARPLOTSFULLDB join("\t", $readLevel, $variety, $methodName, $evaluationLevel, $callRate, $percOK_madeCall_fullAccuracy), "\n";						
						}
					}
					print READSCORRECTBYLEVEL join("\t", @output_fields_byLevelCorrect), "\n";
				}
			}
			
			close(READSCORRECTBYLEVEL);
		}
	}
	
	close(BARPLOTSFULLDB);

	{
		my %_methods;
		my %_readStratification;
		my %_evaluationLevels;
		foreach my $variety (@varieties)
		{
			foreach my $methodName (keys %{$freq_byVariety_byLevel{$variety}})
			{
				$_methods{$methodName}++;
				foreach my $level (keys %{$freq_byVariety_byLevel{$variety}{$methodName}})
				{
					$_evaluationLevels{$level}++;
				}
			}
		}
		
		my @methods = sort keys %_methods;
		
		my @evaluationLevels = qw/species genus family/;
		# my @evaluationLevels = sort keys %_evaluationLevels;
		die Dumper("Missing evaluation levels", \@evaluationLevels, \%_evaluationLevels) unless(all {exists $_evaluationLevels{$_}} @evaluationLevels);
				
		open(FREQEVALUATION, '>', '_frequenciesCorrectByLevel') or die;
		my @header_fields_1_freqCorrect = ('EvaluationLevel');
		my @header_fields_2_freqCorrect = ('');
		my @header_fields_3_freqCorrect = ('');
		
		foreach my $variety (@varieties)
		{
			my $hf2_before = $#header_fields_2_freqCorrect;
			foreach my $method (@methods)
			{
				push(@header_fields_2_freqCorrect, $method, '');	
				push(@header_fields_3_freqCorrect, 'fCorrect', 'L1');					
			}
			
			my $hf2_after = $#header_fields_2_freqCorrect;
			my $requiredFields = $hf2_after - $hf2_before;
			die unless($requiredFields > 0);
			my @addToHeader1 = ($variety, (('') x ($requiredFields - 1)));
			die unless(scalar(@addToHeader1) == $requiredFields);
			push(@header_fields_1_freqCorrect, @addToHeader1);
		}
		
		print FREQEVALUATION join("\t", @header_fields_1_freqCorrect), "\n";
		print FREQEVALUATION join("\t", @header_fields_2_freqCorrect), "\n";
		print FREQEVALUATION join("\t", @header_fields_3_freqCorrect), "\n";
		 
		foreach my $evaluationLevel (@evaluationLevels)
		{		
			my @output_fields_freqCorrect = ($evaluationLevel);	
		
			foreach my $variety (@varieties)
			{				
				foreach my $methodName (@methods)
				{				 
					my $freqOK = 'NA';
					my $M_AVGRE = 'NA';
					my $M_RRMSE = 'NA';
					my $M_L1 = 'NA';
					my $M_L2 = 'NA';
					
					if(exists $freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel})
					{
						$freqOK = $freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{correct}/$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{total};
						$M_AVGRE = sum(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{AVGRE}}) / scalar(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{AVGRE}});
						$M_RRMSE = sum(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{RRMSE}}) / scalar(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{RRMSE}});
						$M_L1 = sum(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{L1}}) / scalar(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{L1}});
						$M_L2 = sum(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{L2}}) / scalar(@{$freq_byVariety_byLevel{$variety}{$methodName}{$evaluationLevel}{L2}});
					}
					
					push(@output_fields_freqCorrect, $freqOK, $M_L1);
				}				
			}	
			
			print FREQEVALUATION join("\t", @output_fields_freqCorrect), "\n";
		}
		
		close(FREQEVALUATION);
	
		open(XYPLOTS, '>', '_forPlot_frequencies_xy') or die;
		print XYPLOTS join("\t", qw/simulationI variety method level taxonID taxonLabel freqTarget freqIs/), "\n";
		for(my $simulationI = 0; $simulationI < $realizedN; $simulationI++)  
		{
			next unless(defined $frequencyComparisons_bySimulation[$simulationI]);
			foreach my $variety (keys %{$frequencyComparisons_bySimulation[$simulationI]})
			{
				foreach my $label (keys %{$frequencyComparisons_bySimulation[$simulationI]{$variety}})
				{
					foreach my $level (keys %{$frequencyComparisons_bySimulation[$simulationI]{$variety}{$label}})
					{
						foreach my $taxonID (keys %{$frequencyComparisons_bySimulation[$simulationI]{$variety}{$label}{$level}})
						{
							my $taxonID_label = (($taxonID eq 'Unclassified') or ($taxonID eq 'NotLabelledAtLevel')) ? $taxonID : taxTree::taxon_id_get_name($taxonID, $fullTaxonomy_simulation);
							print XYPLOTS join("\t", $simulationI, $variety, $label, $level, $taxonID, $taxonID_label, @{$frequencyComparisons_bySimulation[$simulationI]{$variety}{$label}{$level}{$taxonID}}), "\n";
						}
					}
				}
			}
		}
		close(XYPLOTS);
	}
	
	if(1 == 0)
	{
		foreach my $variety (@varieties)
		{
			print $variety, "\n";
			print "\tReads correctly assigned:\n";
			foreach my $methodName (keys %{$n_reads_correct_byVariety{$variety}})
			{
				print "\t\t", $methodName, "\n";
				foreach my $category (keys %{$n_reads_correct_byVariety{$variety}{$methodName}})
				{
					print "\t\t\t", $category, "\n";
					print "\t\t\t\t", "Missing:  ", $n_reads_correct_byVariety{$variety}{$methodName}{$category}{missing}, "\n";
					print "\t\t\t\t", "N      :  ", $n_reads_correct_byVariety{$variety}{$methodName}{$category}{N}, "\n";
					print "\t\t\t\t", "Correct:  ", $n_reads_correct_byVariety{$variety}{$methodName}{$category}{correct}, "\n";
				}
			}
			
			print "\tReads correctly assigned, by level:\n";
			foreach my $methodName (keys %{$n_reads_correct_byVariety_byLevel{$variety}})
			{
				print "\t\t", $methodName, "\n";
				foreach my $category (keys %{$n_reads_correct_byVariety_byLevel{$variety}{$methodName}})
				{
					print "\t\t\t", $category, "\n";			
					foreach my $level (keys %{$n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category}})
					{
						print "\t\t\t\t", $level, "\n";			
						print "\t\t\t\t\t", "Missing   :  ", $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category}{$level}{missing}, "\n";
						print "\t\t\t\t\t", "N         :  ", $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category}{$level}{N}, "\n";
						print "\t\t\t\t\t", "N_td      :  ", $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category}{$level}{N_truthDefined}, "\n";
						print "\t\t\t\t\t", "Correct   :  ", $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category}{$level}{correct}, "\n";
						print "\t\t\t\t\t", "Correct_tD:  ", $n_reads_correct_byVariety_byLevel{$variety}{$methodName}{$category}{$level}{correct_truthDefined}, "\n";
					}
				}
			}		
		}
	}

}
else
{
	die "Unknown --action";
}

sub createArrayJob
{
	my $in = shift;
	my $out = shift;
	my $maxMem = shift;
	my $queue = shift;
	
	if(not defined $maxMem)
	{
		$maxMem = 32;
	}		
	if(not defined $queue)
	{
		$queue = 'phillippy.q';
	}		
	
	my @cmds;
	open(IN, '<', $in) or die "Cannot open $in";
	while(<IN>)
	{
		my $l = $_;
		chomp($l);
		push(@cmds, $l);
	}
	close(IN);
	
	my $nJobs = scalar(@cmds);
	
	my @cmds_switch;
	for(my $cmdI = 0; $cmdI <= $#cmds; $cmdI++)
	{
		my $cmd = $cmds[$cmdI];
		my $cmd_switch = 
qq(if [ "\$jobID" == "$cmdI" ]
then $cmd
fi
);
		push(@cmds_switch, $cmd_switch);
	}

	open(OUT, '>', $out) or die "Cannot open $out";
print OUT qq(#!/bin/bash
#\$ -t 1-${nJobs}
#\$ -q $queue
#\$ -l mem_free=${maxMem}G
jobID=\$(expr \$SGE_TASK_ID - 1)
cd $FindBin::Bin
);
	print OUT join("\n", @cmds_switch);
	
	close(OUT);
	
	return scalar(@cmds);
}

sub inferenceOneSimulation
{
	my $simulation_href = shift;
	
	my $fullTaxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
	
	for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++)
	{		
		if(defined $useVarietyI)
		{
			next unless($varietyI eq $useVarietyI);
		}
		
		print "Inference $simulation_href->{inferenceDBs}[$varietyI][2] \n";
			
		my $DB_target_dir = $simulation_href->{outputDirectory} . '/DB_' . $simulation_href->{inferenceDBs}[$varietyI][2];
		die unless($DB_target_dir eq $simulation_href->{dbDirs_metamap}[$varietyI]);
		
		my $inference_target_dir = $simulation_href->{outputDirectory} . '/inference_' . $simulation_href->{inferenceDBs}[$varietyI][2];
		
		(mkdir($inference_target_dir) or die "Cannot mkdir $inference_target_dir") unless(-d $inference_target_dir);
		
		print "Doing inference in $DB_target_dir\n";
		 
		doMetaMap($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq});
		SimulationsKraken::doKraken($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $krakenDBTemplate, $kraken_binPrefix, $Bracken_dir);
		if($simulation_href->{inferenceDBs}[$varietyI][2] eq 'fullDB')
		{
			# SimulationsMetaPalette::doMetaPalette($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $metaPalette_installation_dir, $jellyfish_2_bin, $masterTaxonomy_dir);			
		}
	} 
	
	# foreach my $oneDBdir (@{$simulation_href->{dbDirs_metamap}})
	# {
	
	# foreach my $inference_variety (@{$simulation_href->{inferenceDBs}})
	# {
		
		# print "Doing inference in $oneDBdir\n";
		# # doMetaMap($simulation_href->{outputDirectory}, $oneDBdir, $simulation_href->{readsFastq});
		# SimulationsKraken::doKraken($simulation_href->{outputDirectory}, $oneDBdir, $simulation_href->{readsFastq}, $krakenDBTemplate, $kraken_binPrefix, $Bracken_dir);
		# #SimulationsMetaPalette::doMetaPalette($simulation_href->{outputDirectory}, $oneDBdir, $simulation_href->{readsFastq}, $metaPalette_installation_dir, $jellyfish_2_bin, $fullTaxonomy);			
	# }
	
}	

sub get_files_for_evaluation
{
	my $simulation_href = shift;
	return (
		'Bracken-Dist' => ['distribution', 'results_bracken.txt'],
		'Kraken-Dist' => ['distribution', 'results_kraken.txt'],
		'MetaMap-EM-Dist' => ['distribution', 'metamap.EM.WIMP'],
		'MetaMap-U-Dist' => ['distribution', 'metamap.U.WIMP'],
		'Kraken-Reads' => ['reads', 'results_kraken.txt.reads2Taxon'],
		'Metamap-U-Reads' => ['reads', 'metamap.U.reads2Taxon'],
		'Metamap-EM-Reads' => ['reads', 'metamap.EM.reads2Taxon'],
		# 'MetaPalette' => ['distribution', 'results_metapalette.txt', 1]
	);
	# return ('Metamap-U' => 'metamap.U.WIMP', 'Metamap-EM' => 'metamap.EM.WIMP');
}

sub evaluateOneSimulation 
{
	my $simulation_href = shift;
	my $n_reads_correct_byVariety = shift;
	my $n_reads_correct_byVariety_byLevel = shift;
	my $freq_byVariety_byLevel = shift;
	my $frequencyComparison_href = shift;
	
	my $taxonomyFromSimulation_dir = $simulation_href->{outputDirectory} . '/DB_fullDB/taxonomy';
	my $taxonomy_usedForSimulation = taxTree::readTaxonomy($taxonomyFromSimulation_dir);
		
	(my $extendedMaster, my $extendedMaster_merged) = validation::prepare_masterTaxonomy_withX($masterTaxonomy_dir, $taxonomy_usedForSimulation);

	my $truth_fn = $simulation_href->{outputDirectory} . '/truth_reads.txt';
	my $truth_raw_reads_href = validation::readTruthFileReads($extendedMaster, $extendedMaster_merged, $truth_fn);
	
	# my %truth_raw_taxonIDs;
	# $truth_raw_taxonIDs{$taxonID_master}++;
	
	my %expected_results_files = get_files_for_evaluation();
	for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++)
	{
		my $varietyName = $simulation_href->{inferenceDBs}[$varietyI][2];
		print "Analyse $varietyName\n";
		die if($varietyName =~ /_1/); # we still need to deal with this and remove the numerical indices for storing
		
		$n_reads_correct_byVariety->{$varietyName} = {} unless(defined $n_reads_correct_byVariety->{$varietyName});
		$n_reads_correct_byVariety_byLevel->{$varietyName} = {} unless(defined $n_reads_correct_byVariety_byLevel->{$varietyName});
		$freq_byVariety_byLevel->{$varietyName} = {} unless(defined $freq_byVariety_byLevel->{$varietyName});
		$frequencyComparison_href->{$varietyName} = {} unless(defined $frequencyComparison_href->{$varietyName});
		
		my $simulation_results_dir = $simulation_href->{outputDirectory} . '/inference_' . $varietyName;
		my $DBdir = $simulation_href->{outputDirectory} . '/DB_' . $varietyName;
		
		# details for reduced DB
		my %reduced_taxonID_original_2_contigs;
		my %reduced_contigLength;
		Util::read_taxonIDs_and_contigs($DBdir, \%reduced_taxonID_original_2_contigs, \%reduced_contigLength);
		
		# read reduced taxonomy
		my %reduced_taxonID_master_2_contigs;
		foreach my $taxonID_original (keys %reduced_taxonID_original_2_contigs)
		{
			my $taxonID_master = taxTree::findCurrentNodeID($extendedMaster, $extendedMaster_merged, $taxonID_original);
			# store which taxon IDs *are* mappable
			$reduced_taxonID_master_2_contigs{$taxonID_master} = $reduced_taxonID_original_2_contigs{$taxonID_original};
		}
		
		# reduce master taxonomy
		my $specificTaxonomy = dclone $extendedMaster;
		taxTree::removeUnmappableParts($specificTaxonomy, \%reduced_taxonID_master_2_contigs);
	
		# translate truth into reduced representation
		my $truth_mappingDatabase_reads = validation::translateReadsTruthToReducedTaxonomy($extendedMaster, $specificTaxonomy, $truth_raw_reads_href);
		
		# get distribution
		my $truth_mappingDatabase_distribution = validation::truthReadsToTruthSummary($specificTaxonomy, $truth_mappingDatabase_reads);

		# die Dumper($varietyName, $truth_mappingDatabase_distribution);
		
		# my %truth_allReads;
		# my %taxonID_translation;
		# foreach my $taxonID (keys %truth_raw_taxonIDs)
		# {
			# $taxonID_translation{$taxonID} = getRank_phylogeny($extendedMaster, $specificTaxonomy, $taxonID);
			# my $count = $truth_raw_taxonIDs{$taxonID};
			# foreach my $rank (keys %{$taxonID_translation{$taxonID}})
			# {
				# my $v = $taxonID_translation{$taxonID}{$rank};
				# $truth_allReads{$rank}{$v} += $count;
			# }
		# }
		# foreach my $rank (keys %truth_allReads)
		# {
			# foreach my $taxonID (keys %{$truth_allReads{$rank}})
			# {
				# $truth_allReads{$rank}{$taxonID} /= scalar(keys %truth_raw_reads);
			# }
		# }
		
		foreach my $methodName (keys %expected_results_files)
		{
			my $methodDetails = $expected_results_files{$methodName};
			my $evaluationType = $methodDetails->[0];
			my $f = $simulation_results_dir . '/' . $methodDetails->[1];
			my $optional = $methodDetails->[2];
			unless(-e $f)
			{
				if($optional)
				{
					warn "Expected file $f for $methodName not present.";
					next;
				}
				else
				{
					die "Expected file $f for $methodName not present.";					
				}
			}
			die unless(($evaluationType eq 'reads') or ($evaluationType eq 'distribution'));
			
			if($evaluationType eq 'reads')
			{
				my $inferred_reads = validation::readInferredFileReads($specificTaxonomy, $extendedMaster_merged, $f);
				unless(all {exists $truth_raw_reads_href->{$_}} keys %$inferred_reads)
				{
					my @missing_readIDs = grep {not exists $truth_raw_reads_href->{$_}} keys %$inferred_reads;
					die Dumper("Missing some reads in truth file $truth_fn (inference file $f)", @missing_readIDs[0  .. 5]);
				}
				validation::readLevelComparison($extendedMaster, $truth_raw_reads_href, $truth_mappingDatabase_reads, $inferred_reads, $methodName, $n_reads_correct_byVariety->{$varietyName}, $n_reads_correct_byVariety_byLevel->{$varietyName});
			}
			else
			{
				print "Read $f\n";
				my $inferred_distribution = validation::readInferredDistribution($extendedMaster, $extendedMaster_merged, $f);
				print "Analyse $f\n";
				validation::distributionLevelComparison($extendedMaster, $truth_mappingDatabase_distribution, $inferred_distribution, $methodName, $freq_byVariety_byLevel->{$varietyName}, $frequencyComparison_href->{$varietyName});
				print "Done $f\n\n";
			}
			
			# unclassified = deliberately no caller
			# undefined = there would be an assignment, but there is no node
		
			# foreach my $level (@evaluateAccuracyAtLevels)
			# {
				# next unless(defined $inference{$level});
				# die unless(defined $truth_allReads{$level});
				# my $totalFreq = 0;
				# my $totalFreqCorrect = 0;
				# foreach my $inferredTaxonID (keys %{$inference{$level}})
				# {
					# my $isFreq = $inference{$level}{$inferredTaxonID}[1];
					# my $shouldBeFreq = 0;
					# if(exists $truth_allReads{$level}{$inferredTaxonID})
					# {
						# $shouldBeFreq = $truth_allReads{$level}{$inferredTaxonID};
					# }
					
					# $totalFreq += $isFreq;
					# if($isFreq <= $shouldBeFreq)					
					# {
						# $totalFreqCorrect += $isFreq;
					# }
					# else
					# {
						# $totalFreqCorrect += $shouldBeFreq;
					# }
				# }
				# die Dumper("Weird total freq", $name, $totalFreq, $f, $level) unless(abs(1 - $totalFreq) <= 1e-3);
				
				# print join("\t", $varietyName, $name, $level, $totalFreqCorrect), "\n";
			# }
		}
			
	}
}

sub addInferenceRoundsWithReducedDBs_allSpecies
{
	my $simulation_href = shift;
	my $taxonomy_href = shift;
	my $fh_log = shift;
	
	print {$fh_log}  "Remove, in turn, each of the contained taxa.\n";
	
	my %n_remove_key_counter;
	foreach my $taxonID2Remove (keys %{$simulation_href->{targetTaxons}})
	{
		print {$fh_log} "\tSelected origin for removal: ", $taxonID2Remove, " (", taxTree::taxon_id_get_name(taxTree::node_get_rank_value($taxonomy_href, $taxonID2Remove, 'superkingdom'), $taxonomy_href), " -- ", taxTree::species_or_strain($taxonomy_href, $taxonID2Remove, 1), ") -- ", taxTree::taxon_id_get_name($taxonID2Remove, $taxonomy_href), "\n";
		
		my $ancestors_href = taxTree::get_ancestors_by_rank($taxonomy_href, $taxonID2Remove);
		$ancestors_href->{self} = $taxonID2Remove;
		
		foreach my $removeRank (qw/self species genus/)
		{
			next unless(exists $ancestors_href->{$removeRank});
			my $removeNode = $ancestors_href->{$removeRank};
			my $remove_key = 'removeOne_' . $removeRank;
			my $already_removed = 0;
			if(exists $n_remove_key_counter{$remove_key})
			{
				$already_removed = $n_remove_key_counter{$remove_key};
			}
			$n_remove_key_counter{$remove_key}++;
			$remove_key .= '_' . $already_removed;
			push(@{$simulation_href->{inferenceDBs}}, ['remove', [$removeNode], $remove_key]);
			print {$fh_log} "\t\tRemoval at level $removeRank: $removeNode -- ", taxTree::taxon_id_get_name($removeNode, $taxonomy_href), "\n";
		}
	}
}


sub addInferenceRoundsWithReducedDBs
{
	my $simulation_href = shift;
	my $taxonomy_href = shift;
	my $fh_log = shift;
	
	my @targetTaxons_shuffled = shuffle keys %{$simulation_href->{targetTaxons}};
	die unless(scalar(@targetTaxons_shuffled) > 1);

	my $removal_origin = $targetTaxons_shuffled[0];
	
	print {$fh_log} "\tSelected origin for removal: ", $removal_origin, " (", taxTree::taxon_id_get_name(taxTree::node_get_rank_value($taxonomy_href, $removal_origin, 'superkingdom'), $taxonomy_href), " -- ", taxTree::species_or_strain($taxonomy_href, $removal_origin, 1), ") -- ", taxTree::taxon_id_get_name($removal_origin, $taxonomy_href), "\n";
	
	my $ancestors_href = taxTree::get_ancestors_by_rank($taxonomy_href, $removal_origin);
	$ancestors_href->{self} = $removal_origin;
	
	foreach my $removeRank (qw/self species genus/)
	{
		next unless(exists $ancestors_href->{$removeRank});
		my $removeNode = $ancestors_href->{$removeRank};
		push(@{$simulation_href->{inferenceDBs}}, ['remove', [$removeNode], 'removeOne_' . $removeRank]);
		print {$fh_log} "\t\tRemoval at level $removeRank: $removeNode -- ", taxTree::taxon_id_get_name($removeNode, $taxonomy_href), "\n";
	}
}

sub returnOneSimulation_noUnknown
{
	my $DB = shift;
	my $taxon_2_contig_href = shift;
	my $taxon_2_genomelength_href = shift;
	my $n_species = shift;
	my $outputDirectory = shift;
	my $equalCoverage = shift;
	my $taxonomy_href = shift;
	my $fh_log = shift;
	
	die unless(defined $taxonomy_href);
	
	die unless($n_species <= scalar(keys %$taxon_2_contig_href));
	
	my @availableTaxonIDs = shuffle(keys %$taxon_2_contig_href);
	
	@availableTaxonIDs = grep {die unless(defined $taxon_2_genomelength_href->{$_}); ($taxon_2_genomelength_href->{$_} >= 100000)} @availableTaxonIDs;
	
	my @selectedTaxa = @availableTaxonIDs[0 .. ($n_species-1)];
	
	print {$fh_log} "\tTaxa: ", join(', ', @selectedTaxa), "\n";
	
	my %targetTaxonIDs;
	if($equalCoverage)
	{
		my $fr_taxon = 1/scalar(@selectedTaxa);
		%targetTaxonIDs = map {$_ => $fr_taxon} @selectedTaxa;
	}
	else
	{
		my @frequencies;
		my $frequencies_sum = 0;
		for(@selectedTaxa)
		{
			my $normal = random_normal();
			my $logNormal = exp($normal);
			push(@frequencies, $logNormal);
			$frequencies_sum += $logNormal;
		}
		foreach my $fr (@frequencies)
		{
			$fr /= $frequencies_sum;
		}
		%targetTaxonIDs = (mesh @selectedTaxa, @frequencies);		
	}
	
	my $simulation_href = {
		coverageFactor => 20,
		outputDirectory => $outputDirectory,
		DB_simulation => $DB,
		targetTaxons => \%targetTaxonIDs,
		inferenceDBs => [['$SAME', undef, 'fullDB']],
	};
	
	return $simulation_href;
}

sub returnOneSimulation_fromFile
{
	my $DB = shift;
	my $taxon_2_contig_href = shift;
	my $taxon_2_genomelength_href = shift;
	my $outputDirectory = shift;
	my $coverageMode = shift;
	my $taxonomy_href = shift;
	my $fh_log = shift;
	my $fn_desired = shift;
	
	die unless(defined $taxonomy_href);
	die unless(defined $fn_desired);
	
	die unless(($coverageMode eq 'equal') or ($coverageMode eq 'logNormal') or ($coverageMode eq 'file'));

	my @selectedTaxa;
	my %selectedTaxa_frequencies;
	open(DESIRED, '<', $desiredTaxa) or die "Cannot open --desiredTaxa $desiredTaxa";
	while(<DESIRED>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless((scalar(@line_fields) == 1) or (scalar(@line_fields) == 2) or (scalar(@line_fields) == 3));
		if($coverageMode eq 'file')
		{
			die "File $desiredTaxa has only one field, but needs two (taxon ID and genome proportion), --coverageMode file in effect" unless(scalar(@line_fields) >= 2);
		}
		else
		{
			die "File $desiredTaxa has two fields, but needs one (taxon ID and genome proportion) unless --coverageMode file" unless(scalar(@line_fields) == 1);
		}
		my $useTaxonID;
		if(exists $taxon_2_contig_href->{$line_fields[0]})
		{
			$useTaxonID = $line_fields[0];
			if($line_fields[2])
			{
				my @contigs = @{$taxon_2_contig_href->{$line_fields[0]}};
				die "Weird - file $desiredTaxa specifies that we want NC $line_fields[2], but this doesn't seem to exist -- available NCs: " . join(', ', @contigs) unless(scalar(grep {index($_, $line_fields[2]) != -1} @contigs));
			}
		}
		else
		{
			die unless(defined $taxonomy_href->{$line_fields[0]}); 
			my @descendants = grep {exists $taxon_2_contig_href->{$_}} taxTree::descendants($taxonomy_href, $line_fields[0]);
			if(scalar(@descendants))
			{
				my @valid_descendants;
				foreach my $tID (@descendants)
				{
					my @contigs = @{$taxon_2_contig_href->{$tID}};
					if(scalar(grep {index($_, $line_fields[2]) != -1} @contigs))
					{
						push(@valid_descendants, $tID);
					}
				}
				if(scalar(@valid_descendants) == 1)
				{
					$useTaxonID = $valid_descendants[0];				
					warn "Taxon ID $line_fields[0] has exactly one genome that's compatible with $line_fields[2] -- use $useTaxonID";
					
				}
				elsif(scalar(@valid_descendants) > 1)
				{
					$useTaxonID = $valid_descendants[0];								
					warn "Taxon ID $line_fields[0] has more than one genome that's compatible with $line_fields[2] -- use $useTaxonID";
				}
				else
				{
					warn "Taxon ID $line_fields[0] has NO genome that's compatible with $line_fields[2]";
				}					
			}
			else
			{
				warn "Unknown taxon ID $line_fields[0] of $desiredTaxa -- skip";			
			}
		}
		
		# print join("\t", $line_fields[0], scalar(@{$taxon_2_contig_href->{$line_fields[0]}})), "\n";
		if(defined $useTaxonID)
		{	
			if(defined $line_fields[2])
			{
				my @contigs = @{$taxon_2_contig_href->{$useTaxonID}};
				die unless(scalar(grep {index($_, $line_fields[2]) != -1} @contigs));				
			}
			if(scalar(@line_fields) >= 2)
			{ 
				$selectedTaxa_frequencies{$useTaxonID} = $line_fields[1];
			}
			push(@selectedTaxa, $useTaxonID);
		}
	}
	close(DESIRED);
	my %_selectedTaxa = map {$_ => 1} @selectedTaxa;
	die "Duplicate selected taxa for simulation" unless(scalar(keys %_selectedTaxa) == scalar(@selectedTaxa));
		
	print {$fh_log} "\tTaxa from file: $fn_desired; ", join(', ', @selectedTaxa), "\n";
	
	my %targetTaxonIDs;
	if($coverageMode eq 'equal')
	{
		my $fr_taxon = 1/scalar(@selectedTaxa);
		%targetTaxonIDs = map {$_ => $fr_taxon} @selectedTaxa;
	} 
	elsif($coverageMode eq 'logNormal')
	{
		my @frequencies;
		my $frequencies_sum = 0;
		for(@selectedTaxa)
		{
			my $normal = random_normal();
			my $logNormal = exp($normal);
			push(@frequencies, $logNormal);
			$frequencies_sum += $logNormal;
		}
		foreach my $fr (@frequencies)
		{
			$fr /= $frequencies_sum;
		}
		%targetTaxonIDs = (mesh @selectedTaxa, @frequencies);		
	}
	elsif($coverageMode eq 'file')
	{
		my $S_fr = 0;
		for(values %selectedTaxa_frequencies)
		{
			$S_fr += $_;
		}
		die unless($S_fr);
		for(@selectedTaxa)
		{
			die unless(defined $selectedTaxa_frequencies{$_});
			$targetTaxonIDs{$_} = $selectedTaxa_frequencies{$_}/$S_fr;
		}	
	}
	
	my $simulation_href = {
		coverageFactor => 20,
		outputDirectory => $outputDirectory,
		DB_simulation => $DB,
		targetTaxons => \%targetTaxonIDs,
		inferenceDBs => [['$SAME', undef, 'fullDB']],
	};
	
	return $simulation_href;
}

sub executeSimulation
{
	my $simulation_href = shift;
	my $taxonID_2_contigIDs = shift;
	my $contig_2_length = shift;
	my $MetaMap_taxonomy = shift;
	my $fh_log = shift;
	my $fh_selfSimilarities = shift;
	my $fh_selfSimilarities_collect = shift;
	my $fh_selfSimilarities_prepareFromOthers = shift;
	
	setup_directory_and_simulate_reads($simulation_href, $MetaMap_taxonomy, $fh_log);
	prepareDBs($simulation_href, $fh_log, $fh_selfSimilarities, $fh_selfSimilarities_collect, $fh_selfSimilarities_prepareFromOthers);
	
	die unless($simulation_href->{dbDirs_metamap});
	die unless($simulation_href->{readsFastq});
	
	store $simulation_href, $simulation_href->{outputDirectory} . '/simulationStore';
	
	# inferenceOneSimulation($simulation_href);

	print "\n\nPrepared simulation!\n\n";
}



sub doMetaMap
{
	my $dir = shift;
	my $DB = shift;
	my $reads = shift;
	
	print "Carry out MetaMap mapping and classification\n";
	print "\tDB   : $DB\n";
	print "\tReads: $reads\n";
	
	my $metamap_output_dir = $dir . '/metamap';
	if(-e $metamap_output_dir)
	{
		if(scalar(glob("$metamap_output_dir/*")))
		{
			system("rm $metamap_output_dir/*") and die "Cannot remove $metamap_output_dir (I)";
		}
		system("rm -r $metamap_output_dir") and die "Cannot remove $metamap_output_dir (II)";
	}
	mkdir($metamap_output_dir) or die "Cannot mkdir $metamap_output_dir";
	
	my $file_mappings = $metamap_output_dir . '/metamap';
	my $file_res_mapping = $file_mappings . '/resources_mapping';
	my $file_res_classification = $file_mappings . '/resources_classification';
	
	die unless(-e $DB.'/DB.fa');
	my $cmd_map = qq(/usr/bin/time -v $metamap_bin mapDirectly --all -r $DB/DB.fa -q $reads -m 2000 --pi 80 -o $file_mappings &> file_res_mapping);
	$cmd_map = qq(/usr/bin/time -v $metamap_bin mapDirectly --all -r $DB/DB.fa -q $reads -m 2000 --pi 80 -o $file_mappings);
	system($cmd_map) and die "Cannot execute $cmd_map";
	die "No MetaMap mappings -- $file_mappings" unless(-e $file_mappings);
		
	my $minReads = 1000;
	if($DB =~ /miniSeq_100/)
	{
		warn "Assume that we're using a very small simulation DB, set --minReads to 1";
		$minReads = 1;
	}	
	my $cmd_classify = qq(/usr/bin/time -v $metamap_bin classify --DB $DB --mappings $file_mappings --minreads $minReads &> file_res_classification);
	$cmd_classify = qq(/usr/bin/time -v $metamap_bin classify --DB $DB --mappings $file_mappings --minreads $minReads); 
	system($cmd_classify) and die "Cannot execute $cmd_classify";
	

	my $resultsFile_EM = $file_mappings . '.EM.WIMP';
	my $resultsFile_EM_reads = $file_mappings . '.EM.reads2Taxon';
	my $resultsFile_U = $file_mappings . '.U.WIMP';
	my $resultsFile_U_reads = $file_mappings . '.U.reads2Taxon';
	
	die "No MetaMap EM classification -- $resultsFile_EM" unless(-e $resultsFile_EM);
	die "No MetaMap EM-reads classification -- $resultsFile_U" unless(-e $resultsFile_EM_reads);	
	die "No MetaMap EM-U classification -- $resultsFile_U" unless(-e $resultsFile_U);
	die "No MetaMap EM-U-reads classification -- $resultsFile_U" unless(-e $resultsFile_U_reads);
	
	copy($resultsFile_EM, $dir);
	copy($resultsFile_EM_reads, $dir);	
	copy($resultsFile_U, $dir);
	copy($resultsFile_U_reads, $dir);
}

sub prepareDBs
{
	my $simulation_href = shift;
	my $fh_log = shift;
	my $fh_selfSimilarities = shift;
	my $fh_selfSimilarities_collect = shift;
	my $fh_selfSimilarities_prepareFromOthers = shift;
	
	foreach my $inference_variety (@{$simulation_href->{inferenceDBs}})
	{
		my $DB_target_dir = $simulation_href->{outputDirectory} . '/DB_' . $inference_variety->[2];
		my $outputFile_truth_inferenceDB = $DB_target_dir . '/truth_inferenceDB.txt';
		
		print {$fh_log} "\t", "Preparing DB variety: $inference_variety->[0]\n";
		
		if($inference_variety->[0] eq '$SAME')
		{
			unless(-d $DB_target_dir)
			{
				mkdir($DB_target_dir) or die "Cannot mkdir $DB_target_dir";
			}
			
			Util::copyMetaMapDB($simulation_href->{DB_simulation}, $DB_target_dir);
			
			# todo perhaps resoncier
			# simulation::truthFileFromReadCounts($outputFile_truth_inferenceDB, \%reads_taxon, $taxonomy_simulation);
		}
		else
		{
			die unless($inference_variety->[0] eq 'remove');
			produceReducedDB($simulation_href->{DB_simulation}, $inference_variety->[1], $DB_target_dir, $fh_selfSimilarities, $fh_selfSimilarities_collect, $fh_selfSimilarities_prepareFromOthers);	
			
			# perhaps reconsider
			# truthFileForReducedInferenceDB();
		}
		
		push(@{$simulation_href->{dbDirs_metamap}}, $DB_target_dir);	
	}
}

sub setup_directory_and_simulate_reads
{
	my $simulation_href = shift;
	my $MetaMap_taxonomy = shift;
	my $fh_log = shift;

	my $verbose = 1;
	
	die unless(defined $simulation_href->{DB_simulation});
	die unless(defined $simulation_href->{coverageFactor});
	die unless(defined $simulation_href->{outputDirectory});
	
	my $simulation_FASTA = $simulation_href->{DB_simulation} . '/DB.fa';
	die unless(-e $simulation_FASTA);
	
	my $contigsTaxa_fn = $simulation_href->{DB_simulation} . '/taxonInfo.txt';
	die unless(-e $contigsTaxa_fn);
	
	my $taxonomy_simulation = {};
	my $contigID_2_taxonID_href = {};
	my $taxonID_2_contigIDs_href = {};	
	my $contig_2_length_href = {};
	fill_contigID_taxonID($simulation_href->{DB_simulation}, $taxonomy_simulation, $contigID_2_taxonID_href, $taxonID_2_contigIDs_href, $contig_2_length_href);
	
	my $outputFile_combinedReads = $simulation_href->{outputDirectory} . '/reads.fastq';
	my $outputFile_rawReads_truth = $simulation_href->{outputDirectory} . '/truth_reads.txt';
	my $outputFile_rawReadCounts = $simulation_href->{outputDirectory} . '/rawreadcounts.txt';
	my $outputFile_truth_readFrequencies_overCompleteTaxonomy = $simulation_href->{outputDirectory} . '/truth_readFrequencies_completeTaxonomy.txt';
	my $outputFile_truth_genomeFrequencies = $simulation_href->{outputDirectory} . '/truth_genomeFrequencies.txt';
	
	my %contigs_and_relative_coverage;
	
	my @targetTaxonIDs = sort keys %{$simulation_href->{targetTaxons}};
	die unless(all {$taxonomy_simulation->{$_}} @targetTaxonIDs );
	die unless(all {$taxonID_2_contigIDs_href->{$_}} @targetTaxonIDs );
	
	# check that target coverages are normalized
	{
		my %coverages_by_level;
		print {$fh_log} "\tTarget coverage proportions:\n";
		my $_s_taxons = 0;
		foreach my $taxonID (@targetTaxonIDs)
		{
			$_s_taxons += $simulation_href->{targetTaxons}{$taxonID};
			print {$fh_log} "\t\t", $taxonID, ": ", $simulation_href->{targetTaxons}{$taxonID}, " (", taxTree::taxon_id_get_name($taxonID, $MetaMap_taxonomy), ")\n";
			my $ancestors_href = taxTree::get_ancestors_with_specific_ranks($MetaMap_taxonomy, $taxonID, [taxTree::getRelevantRanks()]);
			foreach my $rank (keys %$ancestors_href)
			{
				$coverages_by_level{$rank}{$ancestors_href->{$rank}} += $simulation_href->{targetTaxons}{$taxonID};
			}
		}
		
		my $print_taxon_level = 'superkingdom';
		print {$fh_log} "\n\t\tBy $print_taxon_level:\n";
		foreach my $taxonID (sort keys %{$coverages_by_level{$print_taxon_level}})
		{	
			print {$fh_log}  "\t\t\t", $taxonID, " / ", taxTree::taxon_id_get_name($taxonID, $MetaMap_taxonomy), ": ", $coverages_by_level{$print_taxon_level}{$taxonID}, "\n";
		}
		die unless(abs(1 - $_s_taxons) <= 1e-3);
		print {$fh_log} "\n";
	}
	
	my $firstTaxon_genomeLength;
	my $firstTaxon_relativeCoverage;
	my %simulationRelativeCoverages;
	my %taxa_genome_lengths;
	foreach my $taxonID (@targetTaxonIDs)
	{
		my $thisTaxon_genomeLength = getGenomeLength($taxonID, $taxonID_2_contigIDs_href, $contig_2_length_href);
		$taxa_genome_lengths{$taxonID} = $thisTaxon_genomeLength;
		my $targetRelativeCoverage = $simulation_href->{targetTaxons}{$taxonID};
		if(not defined $firstTaxon_genomeLength)
		{
			$firstTaxon_genomeLength = $thisTaxon_genomeLength;
			$firstTaxon_relativeCoverage = $targetRelativeCoverage;
			$simulationRelativeCoverages{$taxonID} = 1;
		}
		else
		{
			my $targetCoverage_relative_first = $targetRelativeCoverage / $firstTaxon_relativeCoverage;
			my $genomeSize_differential = $thisTaxon_genomeLength / $firstTaxon_genomeLength;
			
			$simulationRelativeCoverages{$taxonID} = $targetCoverage_relative_first / $genomeSize_differential;
		}
	}

	if(1)
	{
		print {$fh_log} "\tGenome-size adjusted coverages:\n";
		foreach my $taxonID (sort keys %simulationRelativeCoverages)
		{
			print {$fh_log}  "\t\t", $taxonID, ": ", sprintf("%.3f", $simulationRelativeCoverages{$taxonID}), " (length ", sprintf("%.3f", $taxa_genome_lengths{$taxonID}/(1024**2)), "M)\n";
		}	
		print {$fh_log} "\n";
	}
	
	my $genome_files_href = extractTaxonIDsFromREF($simulation_FASTA, $simulation_href->{outputDirectory} . '/forSimulation', \%simulationRelativeCoverages, $taxonID_2_contigIDs_href);
	
	my $totalReads_allTaxa = 0;
	my %reads_taxon;
	my %readID_2_taxon;
	my %taxonID_2_bases;
	for(my $taxonI = 0; $taxonI <= $#targetTaxonIDs; $taxonI++)
	{
		my $taxonID = $targetTaxonIDs[$taxonI];
		
		my $fn_genome = $genome_files_href->{$taxonID};
		my $relativeCoverage = $simulationRelativeCoverages{$taxonID};
		die unless((defined $fn_genome) and (defined $relativeCoverage));
		$relativeCoverage *= $simulation_href->{coverageFactor};
		
		my $outputDir_reads = $fn_genome . '.reads';
		my $cmd_rm_outputDir_reads = "rm $outputDir_reads/*; rm -r $outputDir_reads";
		
		if(-e $outputDir_reads)
		{
			system($cmd_rm_outputDir_reads) and die "Couldn't execute: $cmd_rm_outputDir_reads";
		}
		
		mkdir($outputDir_reads) or die "Cannot mkdir $outputDir_reads";
		
		my $PBsim_thisSimulation = $PBsim_cmd;
		$PBsim_thisSimulation =~ s/REF/$fn_genome/;
		$PBsim_thisSimulation =~ s/DEPTH/$relativeCoverage/;
		$PBsim_thisSimulation =~ s/--prefix/--prefix $outputDir_reads\/ /;

		print "Executing simulation command:\n\t $PBsim_thisSimulation\n\n";
		system($PBsim_thisSimulation) and die "Command $PBsim_thisSimulation failed";

		my @files_reads = glob("$outputDir_reads/*.fastq");
		die "No reads files in $outputDir_reads ?" unless(scalar(@files_reads));

		my $doAppend = ($taxonI > 0);
		my $readIDs_aref = [];
		my $combinedReadLength = 0;
		my $reads_taxon = combineFASTQ(\@files_reads, $outputFile_combinedReads, $taxonI, $doAppend, $readIDs_aref, \$combinedReadLength);
		$reads_taxon{$taxonID} = $reads_taxon;
		$totalReads_allTaxa += $reads_taxon;
		foreach my $readID (@$readIDs_aref)
		{
			die if(defined $readID_2_taxon{$readID});
			$readID_2_taxon{$readID} = $taxonID;
		}
		system($cmd_rm_outputDir_reads) and die "Couldn't execute: $cmd_rm_outputDir_reads";
		
		unlink($fn_genome) or die "Cannot delete $fn_genome";
		
		$taxonID_2_bases{$taxonID} += $combinedReadLength;
	}
	
	my %realizedTaxonProportions;
	foreach my $taxonID (@targetTaxonIDs)
	{
		$realizedTaxonProportions{$taxonID} = $reads_taxon{$taxonID} / $totalReads_allTaxa;
	}
	
	if(1)
	{
		print {$fh_log} "\tRealized coverages:\n";
		foreach my $taxonID (@targetTaxonIDs)
		{
			print {$fh_log} "\t\t", $taxonID, ": ", $reads_taxon{$taxonID}, " reads, proportion ", sprintf("%.3f", $realizedTaxonProportions{$taxonID}), "; vs target ", sprintf("%.3f", $simulation_href->{targetTaxons}{$taxonID}), "\n";
		}	
		print {$fh_log} "\n";
	}
	
	open(READCOUNTS, '>', $outputFile_rawReadCounts) or die "Cannot open $outputFile_rawReadCounts";
	foreach my $taxonID (@targetTaxonIDs)
	{
		print READCOUNTS $taxonID, "\t", $reads_taxon{$taxonID}, "\n";
	}		
	close(READCOUNTS);
	
	open(READSTRUTH, '>', $outputFile_rawReads_truth) or die "Cannot open $outputFile_rawReads_truth";
	foreach my $readID (keys %readID_2_taxon)
	{
		print READSTRUTH join("\t", $readID, $readID_2_taxon{$readID}), "\n";
	}
	close(READSTRUTH);
	
	simulation::truthReadFrequenciesFromReadCounts($outputFile_truth_readFrequencies_overCompleteTaxonomy, \%reads_taxon, $taxonomy_simulation);
	simulation::truthGenomeFrequenciesFromReadCounts($outputFile_truth_genomeFrequencies, \%taxonID_2_bases, \%reads_taxon, \%taxa_genome_lengths, $taxonomy_simulation);
	
	$simulation_href->{readsFastq} = $outputFile_combinedReads;
}

sub produceReducedDB
{
	my $baseDB = shift;
	my $removeNodes_aref = shift;
	my $targetDir = shift;
	my $fh_selfSimilarities = shift;
	my $fh_selfSimilarities_collect = shift;
	my $fh_selfSimilarities_prepareFromOthers = shift;
	
	print "produceReducedDB(..)\n";
	print "\tBase  : ", $baseDB, "\n";
	print "\tRemove: ", join(', ', @{$removeNodes_aref}), "\n";
	print "\tTarget: ", $targetDir, "\n";
	
	unless(-d $targetDir)
	{
		mkdir($targetDir) or die "Cannot mkdir $targetDir";
	}
	
	my $taxonomy_base = taxTree::readTaxonomy($baseDB . '/taxonomy');
	die unless(all {exists $taxonomy_base->{$_}} @$removeNodes_aref);	
	
	my %removeNodes_and_descendants = map {$_ => 1} (map {($_, taxTree::descendants($taxonomy_base, $_))} @$removeNodes_aref);
	
	my $reducedTaxonomy = taxTree::cloneTaxonomy_removeNodes($taxonomy_base, \%removeNodes_and_descendants);
	die Dumper(\%removeNodes_and_descendants, $removeNodes_aref) unless(all {not exists $reducedTaxonomy->{$_}} @$removeNodes_aref);	
	
	# copy taxonomy
	Util::copyMetaMapTaxonomy($baseDB . '/taxonomy', $targetDir . '/taxonomy');
	
	# copy reduced taxon info
	open(TAXONIN, '<', $baseDB . '/taxonInfo.txt') or die;
	open(TAXONOUT, '>', $targetDir . '/taxonInfo.txt') or die;
	while(<TAXONIN>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/ /, $line);
		die unless(scalar(@f) == 2);
		if(exists $reducedTaxonomy->{$f[0]})
		{
			print TAXONOUT $line, "\n";
		}
	}
	close(TAXONIN);
	close(TAXONOUT);
	
	# copy taxon-windows file
	copy($baseDB . '/contigNstats_windowSize_1000.txt', $targetDir . '/contigNstats_windowSize_1000.txt');

	# create reduced FASTA
	my $baseDB_fa = $baseDB . '/DB.fa';
	my $reducedDB_fa = $targetDir . '/DB.fa';
	{
		open(DBIN, '<', $baseDB_fa) or die "Cannot open $baseDB_fa";
		open(DBOUT, '>', $reducedDB_fa) or die "Cannot open $reducedDB_fa";
		my $relevant = 0;
		while(<DBIN>)
		{
			if(substr($_, 0, 1) eq '>')
			{
				my $taxonID = Util::extractTaxonID($_, $baseDB_fa, $.);
				die unless(exists $taxonomy_base->{$taxonID});
				if(exists $reducedTaxonomy->{$taxonID})
				{
					$relevant = 1;
					print DBOUT $_;
				}
				else
				{
					$relevant = 0;
				}
			}
			else
			{
				if($relevant)
				{
					print DBOUT $_;
				}
			}
		}
		close(DBIN);
	}
	
	die unless(-e 'estimateSelfSimilarity.pl');
	

	my $cmd_self_similarity = qq(perl estimateSelfSimilarity.pl --DB $targetDir --templateDB $baseDB --mode prepareFromTemplate);
	my $cmd_self_similarity_qsub = qq(perl estimateSelfSimilarity.pl --DB $targetDir --mode prepareFromScratch --autoSubmit 1);
	my $cmd_self_similarity_collect = qq(perl estimateSelfSimilarity.pl --DB $targetDir --mode collect);
	
	# warn  "Check command:\n$cmd_self_similarity\n\n";
	# system($cmd_self_similarity) and die "Command $cmd_self_similarity failed";

	# print {$fh_selfSimilarities} $cmd_self_similarity_qsub, "\n";	
	# print {$fh_selfSimilarities_collect} $cmd_self_similarity_collect, "\n";	
	print {$fh_selfSimilarities_prepareFromOthers} $cmd_self_similarity, "\n";	
	
	print "\t\tCreated reduced DB ", $targetDir, " with ", scalar(keys %$reducedTaxonomy), " nodes instead of ", scalar(keys %$taxonomy_base), " nodes\n";

}

sub truthFileForReducedInferenceDB
{
	die;
}


sub combineFASTQ
{
	my $in_aref = shift;
	my $out_fn = shift;
	my $readID_prefix = shift;
	my $append = shift;
	my $readIDs_aref = shift;
	my $totalBases_sref = shift;
	
	my $totalReads = 0;
	if(defined $totalBases_sref)
	{
		$$totalBases_sref = 0;
	}
	if($append)
	{
		open(OUT, '>>', $out_fn) or die "Cannot open $out_fn";
	}
	else
	{
		open(OUT, '>', $out_fn) or die "Cannot open $out_fn";
	}
	foreach my $inF (@$in_aref)
	{
		open(F, '<', $inF) or die "Cannot open $inF";
		while(<F>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			die unless(substr($line, 0, 1) eq '@');
			$totalReads++;
			
			my $line_seq = <F>;
			chomp($line_seq);
			my $line_plus = <F>;
			my $line_qual = <F>;
			die unless(substr($line_plus, 0, 1) eq '+');
			
			substr($line, 0, 1) = ('@' . $readID_prefix);
			
			print OUT $line, "\n", $line_seq, "\n", $line_plus, $line_qual;
			
			if(defined($readIDs_aref))
			{
				push(@$readIDs_aref, substr($line, 1));
			}
			
			if(defined $totalBases_sref)
			{
				$$totalBases_sref += length($line_seq);
			}
		}
		close(F);
	}
	
	close(OUT);
	
	return $totalReads;
}

sub extractTaxonIDsFromREF
{
	my $fn_in = shift;
	my $fn_out_prefix = shift;
	my $extractTaxonIDs_href = shift;
	my $taxon_2_contig_href = shift;

	my %forReturn_files_by_taxonID;
	my %fh_per_taxonID;
	
	my %target_contig_2_taxonID;
	die unless(scalar(keys %$extractTaxonIDs_href));
	my @taxonIDs = sort keys %$extractTaxonIDs_href;
	for(my $i = 0; $i <= $#taxonIDs; $i++)
	{
		my $taxonID = $taxonIDs[$i];
		
		my $fnOut = $fn_out_prefix . '.' . $i;
		my $fh;
		open($fh, '>', $fnOut) or die "Cannot open $fnOut";
		
		$fh_per_taxonID{$taxonID} = $fh;
		$forReturn_files_by_taxonID{$taxonID} = $fnOut;
		
		die unless(defined $taxon_2_contig_href->{$taxonID});
		foreach my $contigID (@{$taxon_2_contig_href->{$taxonID}})
		{
			$target_contig_2_taxonID{$contigID} = $taxonID;
		}
	}	
	
	my $printedContigs = 0;
	my %sawContigs;
	open(IN, '<', $fn_in) or die "Cannot open $fn_in";
	my $relevant_taxonID = 0;
	while(<IN>)
	{
		my $line = $_;
		if(substr($line, 0, 1) eq '>')
		{
			my $contigID = substr($line, 1);
			chomp($contigID);

			if(exists $target_contig_2_taxonID{$contigID})
			{
				$printedContigs++;
				$relevant_taxonID = $target_contig_2_taxonID{$contigID};
				$sawContigs{$contigID} = 1;
				
				die unless(defined $fh_per_taxonID{$relevant_taxonID});
				print {$fh_per_taxonID{$relevant_taxonID}} $line;
			}
			else
			{
				$relevant_taxonID =  0;
			}
		}
		else
		{
			if($relevant_taxonID)
			{
				die unless(defined $fh_per_taxonID{$relevant_taxonID});
				print {$fh_per_taxonID{$relevant_taxonID}} $line;			
			}
		}
	}
	close(OUT);
	close(IN);
	
	print "Printed $printedContigs from $fn_in --> ${fn_out_prefix}.* \n";
	
	foreach my $contigID (keys %target_contig_2_taxonID)
	{
		die "Missed positive selection for contig $contigID in file $fn_in" unless($sawContigs{$contigID});
	}

	foreach my $taxonID (keys %fh_per_taxonID)
	{
		close($fh_per_taxonID{$taxonID});
	}
	
	return \%forReturn_files_by_taxonID;
}	


sub getGenomeLength
{
	my $taxonID = shift;
	my $taxon_2_contig = shift;
	my $contig_2_length = shift;
	
	my $gL = 0;
	die "Cannot determine genome length for taxon ID $taxonID" unless(defined $taxon_2_contig->{$taxonID});
	foreach my $contigID (@{$taxon_2_contig->{$taxonID}})
	{
		die unless(defined $contig_2_length->{$contigID});
		$gL += $contig_2_length->{$contigID};
	}
	
	return $gL;
}	

sub fill_contigID_taxonID
{
	my $DB = shift;
	
	my $taxonomy_href = shift;	
	my $contigID_2_taxonID_href = shift;	
	my $taxonID_2_contigIDs_href = shift;
	my $contig_2_length_href = shift;

	my $contigsTaxa = $DB . '/taxonInfo.txt';
	die unless(defined $contig_2_length_href);
	die unless(-e $contigsTaxa);
	
	
	my $taxonomy_href_inner = taxTree::readTaxonomy($DB . '/taxonomy');
	%{$taxonomy_href} = %{$taxonomy_href_inner};
		
	open(TAXA, '<', $contigsTaxa) or die "Cannot open $contigsTaxa";
	while(<TAXA>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/ /, $line, -1);
		die unless(scalar(@line_fields) == 2);
		
		my $taxonID = $line_fields[0];
		my $contigs = $line_fields[1];
		
		my @contigs = split(/;/, $contigs);
		
		foreach my $contigInfo (@contigs)
		{
			my @contigFields = split(/=/, $contigInfo);
			die unless(scalar(@contigFields)== 2);
			my $contigID = $contigFields[0];
			my $contigLength = $contigFields[1];
			
			die "Undefined taxon ID $taxonID" unless(defined $taxonomy_href->{$taxonID});

			die if(defined $contigID_2_taxonID_href->{$contigID});
			$contigID_2_taxonID_href->{$contigID} = $taxonID;

			push(@{$taxonID_2_contigIDs_href->{$taxonID}}, $contigID);
			
			die if(defined $contig_2_length_href->{$contigID});
			$contig_2_length_href->{$contigID} = $contigLength;
		}	
	}
	close(TAXA);	
}

sub flagFn
{
	my $directory = shift;
	my $flagName = shift;
	my $fn = $directory . '/flag_' . $flagName . '.txt';
	return $fn;
}

sub getFlag
{
	my $directory = shift;
	my $flagName = shift;
	my $fn = flagFn($directory, $flagName);
	if(-e $fn)
	{
		open(F, '<', $fn) or die;
		my $line = <F>;
		chomp($line);
		close(F);
		return $line;
	}
	else
	{
		return 0;
	}
}

sub setFlag
{
	my $directory = shift;
	my $flagName = shift;
	my $value = shift;
	
	my $fn = flagFn($directory, $flagName);
	open(F, '>', $fn) or die;
	print F $value;
	close(F);
}


__END__


if($action eq 'prepare')
{
	my @files_in_outputDir = glob($globalOutputDir . '/*');
	unless(scalar(@files_in_outputDir) == 0)
	{
		unless($really)
		{
			die "Directory $globalOutputDir not empty, assume already prepared, unless --really 1";
		}
		foreach my $obj (@files_in_outputDir)
		{
			remove_tree($obj);
		}
	}
	
	my $DB_fa = $DB . '/DB.fa';
	my $existingSelfSimilarities = $DB . '/selfSimilarities.txt';
	my $contigsTaxa = $DB . '/taxonInfo.txt';
	
	die unless(-e $DB_fa);
	die unless(-e $existingSelfSimilarities);
	die unless(-e $contigsTaxa);
	
	my $MetaMap_taxonomy = taxTree::readTaxonomy($DB . '/taxonomy');
	my %contigID_2_taxonID;
	my %taxonID_2_contigIDs;	
	my %contig_2_length;
	fill_contigID_taxonID($contigsTaxa, $MetaMap_taxonomy, \%contigID_2_taxonID, \%taxonID_2_contigIDs, \%contig_2_length);

	my @directly_mappable_taxonIDs_bacteria = grep {my $tI = taxTree::get_taxon_id_information($_, $taxonomy); ($tI->{'superkingdom'} eq 'Bacteria <prokaryotes>')} keys %taxonID_2_contigIDs;

	print "Have ", scalar(@directly_mappable_taxonIDs_bacteria), " directly mappable bacterial genomes.\n";

	my %mappingTargets_indirect;
	open(TREE, '<', $treeSelfSimilarityReadsMany) or die "Cannot open $treeSelfSimilarityReadsMany";
	while(<TREE>)
	{
		chomp;
		next unless($_);
		my $line = $_;
		my @fields = split(/\t/, $line, -1);
		die unless (scalar(@fields) > 4);
		my $treeID = $fields[0];
		die unless(defined $MetaMap_taxonomy->{$treeID});
		$mappingTargets_indirect{$treeID} = 1;
	}
	close(TREE);

	my @innernodes_for_possible_attachment = keys %mappingTargets_indirect;
	@innernodes_for_possible_attachment = grep {not $taxonID_2_contigIDs{$_}} @innernodes_for_possible_attachment;
	
	print "Have ", scalar(@innernodes_for_possible_attachment), " nodes for possible indirect attachment.\n";

	my @innernodes_bacteria_possible_attachment = grep {my $tI = taxTree::get_taxon_id_information($_, $taxonomy); ($tI->{'superkingdom'} eq 'Bacteria <prokaryotes>')} @innernodes_for_possible_attachment; 
	print "Have ", scalar(@innernodes_bacteria_possible_attachment), " bacterial nodes for possible attachment.\n";

	@innernodes_bacteria_possible_attachment = grep {(scalar(@{$taxonomy->{$_}{children}})) and (scalar(get_mappable_leaf_taxon_IDs($_, $taxonomy, \%taxonID_2_contigIDs)) > 2)} @innernodes_bacteria_possible_attachment;

	print "Have ", scalar(@innernodes_bacteria_possible_attachment), " bacterial nodes with more than 2 immediate ancestors and more than 2 mappable leafs.\n";

	@innernodes_bacteria_possible_attachment = shuffle @innernodes_bacteria_possible_attachment;

	print "\n";
	for(my $simulationI = 0; $simulationI < $create_simulations; $simulationI++)
	{
		my $outputDir = $globalOutputDir . '/' . $simulationI . '/';
		my $outputDir_reads = $outputDir . '/reads/';
		my $outputFile_combinedReads = $outputDir . '/combined_reads.fastq';
		my $outputFile_truth = $outputDir . '/truth.txt';
		my $outputFile_removedSpeciesAndContigs = $outputDir . '/toRemove.txt';
		
		unless(-d $outputDir)
		{
			mkdir($outputDir) or die "Cannot mkdir $outputDir";
		}
		unless(-d $outputDir_reads)
		{
			mkdir($outputDir_reads) or die "Cannot mkdir $outputDir_reads";
		}	
		
		my $nodeID_unknown_attachTo = $innernodes_bacteria_possible_attachment[$simulationI];
		my @children_possible_removal = @{$taxonomy->{$nodeID_unknown_attachTo}{children}};
		@children_possible_removal = shuffle @children_possible_removal;
		my $removeNodeIDFromReference = $children_possible_removal[0];
		
		my @remove_genome_IDs = get_mappable_leaf_taxon_IDs($removeNodeIDFromReference, $taxonomy, \%taxonID_2_contigIDs);
		@remove_genome_IDs = shuffle(@remove_genome_IDs); die unless(scalar(@remove_genome_IDs) >= 1);
		
		my %remove_genome_IDs = map {$_ => 1} @remove_genome_IDs;
		my @remove_contig_IDs = map {die "Don't have mappable contigs for $_" unless ($taxonID_2_contigIDs{$_}); @{$taxonID_2_contigIDs{$_}}} @remove_genome_IDs;
		die unless(scalar(@remove_contig_IDs) >= scalar(@remove_genome_IDs));
		my %remove_contig_IDs = map {$_ => 1} @remove_contig_IDs;
		
		my @simulate_reads_from_taxonIDs;
		
		print "Simulation $simulationI, remove node $removeNodeIDFromReference (attached to $nodeID_unknown_attachTo) from reference, affecting ", scalar(@remove_genome_IDs), " leaf nodes and ", scalar(@remove_contig_IDs), " contigs.\n";
		my $unknown_selected_genome = $remove_genome_IDs[0];
		push(@simulate_reads_from_taxonIDs, $unknown_selected_genome);
		
		print "\tFor node $removeNodeIDFromReference, we select ID $unknown_selected_genome for simulation-from-unknown.\n";
		
		my @knownGenomes_thisIteration = @directly_mappable_taxonIDs_bacteria;
		@knownGenomes_thisIteration = grep {not exists $remove_genome_IDs{$_}} @knownGenomes_thisIteration;
		@knownGenomes_thisIteration = shuffle @knownGenomes_thisIteration;
		
		print "\tRemaining bacterial genomes in DB after removal: ", scalar(@knownGenomes_thisIteration), " -- of which we select $perSimulation_totalKnownGenomes\n";
		my @selectedKnownGenomes = @knownGenomes_thisIteration[0 .. ($perSimulation_totalKnownGenomes-1)];
		print "\t\tSelection: ", join(', ', @selectedKnownGenomes), "\n";
		push(@simulate_reads_from_taxonIDs, @selectedKnownGenomes);
		
		my @simulate_reads_from_contigIDs = map {die "Don't have mappable contigs for $_" unless ($taxonID_2_contigIDs{$_}); @{$taxonID_2_contigIDs{$_}}} @simulate_reads_from_taxonIDs;
		print "\tSelect ", scalar(@simulate_reads_from_taxonIDs), " genomes for simulation --> ", scalar(@simulate_reads_from_contigIDs), " contigs.\n";

		# print truth files
		
		open(TRUTH, '>', $outputFile_truth) or die "Cannot open $outputFile_truth";
		foreach my $taxonID (@simulate_reads_from_taxonIDs)
		{
			my $unknownSpecies = (exists $remove_genome_IDs{$taxonID}) ? 1 : 0;
			print TRUTH join("\t", $taxonID, $unknownSpecies), "\n";
		}
		close(TRUTH);
		
		open(TOREMOVE, '>', $outputFile_removedSpeciesAndContigs) or die "Cannot open $outputFile_removedSpeciesAndContigs";
		foreach my $taxonID (@remove_genome_IDs)
		{
			foreach my $contigID (@{$taxonID_2_contigIDs{$taxonID}})
			{
				die unless($remove_contig_IDs{$contigID});
				print TOREMOVE join("\t", $taxonID, $contigID), "\n";				
			}
		}
		close(TOREMOVE);
		
		# generate simulated reads
		
		my $fasta_for_simulation = $outputDir . '/simulateReadsFrom.fa';
		selectContigsFromFASTA($ref, $fasta_for_simulation, {map {$_ => 1} @simulate_reads_from_contigIDs}, {});

		my $PBsim_thisSimulation = $PBsim_cmd;
		$PBsim_thisSimulation =~ s/REF/$fasta_for_simulation/;
		$PBsim_thisSimulation =~ s/--prefix/--prefix $outputDir_reads /;

		print "Executing simulation command:\n\t $PBsim_thisSimulation\n\n";
		system($PBsim_thisSimulation) and die "Command $PBsim_thisSimulation failed";

		my @files_reads = glob("$outputDir_reads/*.fastq");
		die "No reads files in $outputDir_reads ?" unless(scalar(@files_reads));
		
		my $cmd_combineReads = 'cat ' . join(' ', @files_reads) . ' > ' . $outputFile_combinedReads;
		system($cmd_combineReads) and die "Couldn't execute $cmd_combineReads";
		
		my $cmd_remove = "rm ${outputDir_reads}/*";
		system($cmd_remove) and die "Couldn't execute $cmd_remove";
	}
}
elsif($action eq 'doAll')
{
	my @dirs_in_outputDir = glob($globalOutputDir . '/*');
	@dirs_in_outputDir = grep {(-d $_) and ($_ =~ /\/\d+$/)} @dirs_in_outputDir;
	foreach my $dir (@dirs_in_outputDir)
	{
		die unless($dir =~ /\/(\d+)$/);
		my $jobI = $1;
		my $cmd = qq(perl createSimulations.pl --action doOneJob --jobI $jobI);
		print $cmd, "\n";
	}			
}	
elsif($action eq 'doOneJob')
{
	die "Please provide --jobI" unless(defined $jobI);
	print "do job $jobI\n";
	print "-jobIMethod: $jobIMethod\n\n";
	die "Please user either all, mashmap, kraken or metapalette as values for --jobIMethod" unless(($jobIMethod eq 'all') or ($jobIMethod eq 'mashmap') or ($jobIMethod eq 'kraken') or ($jobIMethod eq 'metapalette'));
	my %remove_contigs;
	
	my $jobDir = $globalOutputDir . '/' . $jobI . '/';
	my $fn_toRemove = $jobDir . '/toRemove.txt';
	open(TOREMOVE, '<', $fn_toRemove) or die "Cannot open $fn_toRemove";
	while(<TOREMOVE>)
	{
		chomp;
		next unless($_);
		my @f = split(/\t/, $_);
		die unless($#f == 1);
		$remove_contigs{$f[1]}++;
	}
	close(TOREMOVE);
	
	my $fasta_for_mapping = $jobDir . '/DB.fa';
	# selectContigsFromFASTA($ref, $fasta_for_mapping, {}, \%remove_contigs) unless(-e $fasta_for_mapping); # todo remove, this is dangerous! 
		
	onejob_mashmap($jobI) if(($jobIMethod eq 'all') or ($jobIMethod eq 'mashmap'));
	onejob_kraken($jobI) if(($jobIMethod eq 'all') or ($jobIMethod eq 'kraken'));
	onejob_metapalette($jobI) if(($jobIMethod eq 'all') or ($jobIMethod eq 'metapalette'));
}
else
{
	die "Please provid valid --action, e.g. prepare or doOneJob (with --jobI)";
}


sub onejob_kraken
{
	my $jobI = shift;
	my $jobDir = abs_path($globalOutputDir . '/' . $jobI . '/');
	
	my $kraken_dir = $jobDir . '/kraken/';
	my $kraken_dir_DB = $jobDir . '/kraken/DB';
	
	if(1 == 0) # todo remove
	{
	if(-e $kraken_dir)
	{
		system("rm -rf $kraken_dir") and die "Cannot delete $kraken_dir";
	}
	
	unless(-e $kraken_dir)
	{
		mkdir($kraken_dir) or die "Cannot open $kraken_dir";
	}
		
	my $fasta_for_mapping = $jobDir . '/DB.fa';	
	my $fasta_for_simulation = $jobDir . '/simulateReadsFrom.fa';
	my $simulatedReads = $jobDir . '/combined_reads.fastq';	
	die unless(all {-e $_} ($fasta_for_mapping, $fasta_for_simulation, $simulatedReads));
	chdir($kraken_dir) or die;  
	
	if(-e 'DB')
	{
		system('rm -rf DB') and die "Cannot rm";
	}
	
	my $cmd_copy_DB = qq(cp -r $krakenDBTemplate DB);
	system($cmd_copy_DB) and die "Cannot cp $krakenDBTemplate";
	die "DB missing" unless(-d 'DB');
	die unless(-e '../DB.fa');

	my $cmd_convert = qq(perl ${FindBin::Bin}/translateMashmapDBToKraken.pl --input ../DB.fa --taxonomyDir $taxonomyDir --krakenTemplate_taxonomy ${krakenDBTemplate}/taxonomy/);
	system($cmd_convert) and die "Could not execute command: $cmd_convert";
	die "Converted mashmap DB (mashmap -> kraken) missing!" unless(-e '../DB.fa.kraken');
	
	my $cmd_build_II = qq(/usr/bin/time -v ${kraken_binPrefix}-build --add-to-library ../DB.fa.kraken --db DB &> output_build_II.txt);
	system($cmd_build_II) and die "Could not execute command: $cmd_build_II";
	
	my $cmd_build_III = qq(export PATH=/data/projects/phillippy/software/jellyfish-1.1.11/bin:\$PATH; /usr/bin/time -v ${kraken_binPrefix}-build --build --threads 16 --db DB &> output_build_III.txt);
	system($cmd_build_III) and die "Could not execute command: $cmd_build_III";

	my $cmd_classify = qq(/usr/bin/time -v ${kraken_binPrefix} --preload --db DB ../combined_reads.fastq > reads_classified);
	system($cmd_classify) and die "Could not execute command: $cmd_classify";
	
	my $cmd_report = qq(/usr/bin/time -v ${kraken_binPrefix}-report --db DB reads_classified > reads_classified_report);
	system($cmd_report) and die "Could not execute command: $cmd_report";
		
	my $cmd_Bracken_selfSimilarity = qq(bash -c '/usr/bin/time -v ${kraken_binPrefix} --db DB --fasta-input --threads=10 <( find -L DB/library \\( -name "*.fna"  -o -name "*.fa" -o -name "*.fasta" \\) -exec cat {} + ) > database_kraken');
	system($cmd_Bracken_selfSimilarity) and die "Could not execute command: $cmd_Bracken_selfSimilarity";

	my $cmd_Bracken_countkMers = qq(/usr/bin/time -v perl ${Bracken_dir}/count-kmer-abundances.pl --db=DB --read-length=2000 --threads=10 database_kraken > database75mers.kraken_cnts);
	system($cmd_Bracken_countkMers) and die "Could not execute command: $cmd_Bracken_countkMers";
	
	
	my $cmd_Bracken_kMerDist = qq(/usr/bin/time -v python ${Bracken_dir}/generate_kmer_distribution.py -i database75mers.kraken_cnts -o database75mers.kraken_cnts.bracken);
	system($cmd_Bracken_kMerDist) and die "Could not execute command: $cmd_Bracken_kMerDist";
	
	foreach my $L (qw/S G F/)
	{
		my $cmd_Bracken_estAbundance = qq(/usr/bin/time -v python ${Bracken_dir}/est_abundance.py -i reads_classified_report -k database75mers.kraken_cnts.bracken -l $L -o reads_classified_report_bracken_${L});
		system($cmd_Bracken_estAbundance) and die "Could not execute command: $cmd_Bracken_estAbundance";
	}
	}
	
	chdir($kraken_dir) or die;  

	create_compatible_file_from_kraken(
		$jobDir . '/results_kraken.txt',
		'DB/taxonomy',
		'reads_classified_report',
	);
	
	create_compatible_file_from_kraken_bracken(
		$jobDir . '/results_bracken.txt',
		'DB/taxonomy',
		'reads_classified_report',
		'reads_classified_report_bracken_S',
		'reads_classified_report_bracken_G',
		'reads_classified_report_bracken_F');
		
	chdir($FindBin::Bin) or die;
}

sub onejob_mashmap
{
	my $jobI = shift;
	my $jobDir = $globalOutputDir . '/' . $jobI . '/';
	
	my $mashmap_dir = $jobDir . '/mashmap/';
	my $mashmap_dir_DB = $jobDir . '/mashmap/DB';
	
	unless(-e $mashmap_dir)
	{
		mkdir($mashmap_dir) or die "Cannot open $mashmap_dir";
	}
	
	unless(-e $mashmap_dir_DB)
	{
		mkdir($mashmap_dir_DB) or die "Cannot open $mashmap_dir_DB";
	}	
	
	my $fasta_for_mapping = $jobDir . '/DB.fa';	
	my $fasta_for_simulation = $jobDir . '/simulateReadsFrom.fa';
	my $simulatedReads = $jobDir . '/combined_reads.fastq';	
	die unless(all {-e $_} ($fasta_for_mapping, $fasta_for_simulation, $simulatedReads));
	chdir($FindBin::Bin) or die;  
	 
	
	my $cmd_DB_construction = qq(perl mashMapSpecies.pl --action buildDB --buildDB_multiFile $fasta_for_mapping --buildDB_outputFASTA ${mashmap_dir_DB}/ref.fa --buildDB_taxonomyDir $taxonomyDir);
	print "Now executing: $cmd_DB_construction\n";
	system($cmd_DB_construction) and die "Command $cmd_DB_construction failed.\n";
	
	my $cmd_analysis = qq(perl mashMapSpecies.pl --action allSteps --dbDir $mashmap_dir_DB --fastq $simulatedReads --minIdentity 80 --minReadLength 2000 --prefix ${mashmap_dir}/mashmap);
	print "Now executing: $cmd_analysis\n";
	system($cmd_analysis) and die "Command $cmd_analysis failed.\n";
	
	my $resultsFile_EM = $mashmap_dir . '/mashmap.EM.WIMP';
	my $resultsFile_U = $mashmap_dir . '/mashmap.U.WIMP';

	my $cmd_copy = qq(cp $resultsFile_EM ${jobDir}/results_mashmap_EM.txt; cp $resultsFile_U ${jobDir}/results_mashmap_U.txt);
	system($cmd_copy) and die "Command $cmd_copy was not successful";
	
	my $cmd_delete = qq(rm -rf $mashmap_dir);
	print "Now DO NOT DO $cmd_delete, but activate at some point";
}

sub selectContigsFromFASTA
{
	my $fn_in = shift;
	my $fn_out = shift;
	my $positiveSelection_href = shift;
	my $negativeSelection_href = shift;
	
	die unless((scalar(keys %$positiveSelection_href)) xor (scalar(keys %$negativeSelection_href)));
	
	my $printedContigs = 0;
	my %sawContigs;
	open(IN, '<', $fn_in) or die "Cannot open $fn_in";
	open(OUT, '>', $fn_out) or die "Cannot open $fn_out";
	my $relevant = 0;
	while(<IN>)
	{
		my $line = $_;
		if(substr($line, 0, 1) eq '>')
		{
			my $contigID = substr($line, 1);
			chomp($contigID);
			if(keys %$positiveSelection_href)
			{
				if(exists $positiveSelection_href->{$contigID})
				{
					$printedContigs++;
					$relevant = 1;
					$sawContigs{$contigID} = 1;
					print OUT $line;
				}
				else
				{
					$relevant =  0;
				}
			}
			else 
			{
				if(exists $negativeSelection_href->{$contigID})
				{
					$relevant = 0;
					$sawContigs{$contigID} = 1;
				}
				else
				{
					$printedContigs++;
					$relevant =  1;
					print OUT $line;
				}			
			}
		}
		else
		{
			print OUT $line if ($relevant);
		}
	}
	close(OUT);
	close(IN);
	
	print "Printed $printedContigs from $fn_in --> $fn_out \n";
	
	if(keys %$positiveSelection_href)
	{
		foreach my $contigID (keys %$positiveSelection_href)
		{
			die "Missed positive selection for contig $contigID in file $fn_in" unless($sawContigs{$contigID});
		}
	}
	else
	{
		foreach my $contigID (keys %$negativeSelection_href)
		{
			die "Missed negative selection against contig $contigID in file $fn_in" unless($sawContigs{$contigID});
		}	
	}			
}	



sub get_mappable_leaf_taxon_IDs
{
	my $nodeID = shift;
	my $taxonomy = shift;
	my $taxonID_2_contigIDs_href = shift;
	
	my @descendant_nodes = taxTree::descendants($taxonomy, $nodeID);
	push(@descendant_nodes, $nodeID);
	
	@descendant_nodes = grep {exists $taxonID_2_contigIDs_href->{$_}} @descendant_nodes;
	
	return @descendant_nodes;
}
