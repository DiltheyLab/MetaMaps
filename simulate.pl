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
use File::stat;
use Time::localtime;

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
my $masterTaxonomy_dir =  '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';

# these are populated globally if required
my $masterTaxonomy; 
my $masterTaxonomy_merged;

#my $ref = '../db/refseq/ref.fa';
#my $contigsTaxa = $ref . '.taxa';

#my $treeSelfSimiliarityReads = '/data/projects/phillippy/projects/mashsim/NCBI/refseq/taxonomy/selfSimilarity/results.reads.byNode';
#(my $treeSelfSimilarityReadsMany = $treeSelfSimiliarityReads) =~ s/results\.reads\.byNode/results.reads.many.byNode/;


my $PBsim_cmd = qq(/data/projects/phillippy/projects/mashsim/PBSIM-PacBio-Simulator/src/pbsim --model_qc /data/projects/phillippy/projects/mashsim/PBSIM-PacBio-Simulator/data/model_qc_clr --data-type CLR --depth DEPTH --prefix --length-mean $simulation_read_length --accuracy-mean 0.88 REF);

my $metamap_bin = './metamaps';
my $metaPalette_installation_dir = qq(/data/projects/phillippy/software/MetaPalette/);
my $jellyfish_2_bin = qq(/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish);

my $kraken_binPrefix = SimulationsKraken::getKrakenBinPrefix();
my $Bracken_dir = SimulationsKraken::getBrackenDir();
my $krakenDBTemplate = SimulationsKraken::getKrakenDBTemplate();
my $kraken2DBTemplate = SimulationsKraken::getKraken2DBTemplate();
my $kraken2_binPrefix = SimulationsKraken::getKraken2BinPrefix();
my $centrifugeBinDir = SimulationsKraken::getCentrifugeDir();

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
# this stuff is necessary for MetaPalette
#system("eval 'module load libs/hdf5/1.8.12'") and die "Could load hdf5";
#system("eval 'module load libs/scipy/0.18.1'") and die "Could load scipy";
	
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
my $targetTotalSimulationInGigabytes;
my $desiredTaxa;
my $coverageTargetsAreOrganismAbundances = 0;
my $ignoreFullDB;
my $skipKraken;
my $FASTA;
my $FASTA_taxon_id;
my $maxMemory;
my $simulationMinReadLength;
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
	'coverageTargetsAreOrganismAbundances:s' => \$coverageTargetsAreOrganismAbundances,
	'targetTotalSimulationInGigabytes:s' => \$targetTotalSimulationInGigabytes,
	'ignoreFullDB:s' => \$ignoreFullDB,
	'skipKraken:s' => \$skipKraken,
	'FASTA:s' => \$FASTA,
	'FASTA_taxon_id:s' => \$FASTA_taxon_id,
	'maxMemory:s' => \$maxMemory,
	'simulationMinReadLength:s' => \$simulationMinReadLength,
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

if($action eq 'testTaxonomy')
{
	my $fullTaxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
	my %ranks;
	foreach my $taxonID (keys %$fullTaxonomy)
	{
		next unless(scalar(@{$fullTaxonomy->{$taxonID}{children}}) == 0);
		my $rank = $fullTaxonomy->{$taxonID}{rank};
		$ranks{$rank}++;
	}
	print "Rank summary:\n";
	foreach my $rank (sort {$ranks{$b} <=> $ranks{$a}} keys %ranks)
	{
		print " - $rank: ", $ranks{$rank}, "\n";
	}

	foreach my $taxonID (keys %$fullTaxonomy)
	{
		die unless(defined $fullTaxonomy->{$taxonID});
		my @ancestors = taxTree::get_ancestors($fullTaxonomy, $taxonID);
		my $haveSpecies = 0;
		my $rank = $fullTaxonomy->{$taxonID}{rank};
		next unless(($rank eq 'subspecies') or ($rank eq 'varietas') or ($rank eq 'forma'));
		$haveSpecies = 1 if($rank eq 'species');		
		foreach my $ancestorID (@ancestors)
		{
			my $rank_ancestor = $fullTaxonomy->{$ancestorID}{rank};
			$haveSpecies = 1 if($rank_ancestor eq 'species');
		}
		unless($haveSpecies)
		{
			warn "No species info for taxon $taxonID";
		}
	}
}
elsif($action eq 'prepare')
{
	die "Not in use at the moment";
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
elsif($action eq 'repeatReadSimulations')
{
	die "Simulations in directory $globalOutputDir not prepared" unless getFlag($globalOutputDir, 'simulated');
	
	my $fn_log = $globalOutputDir. '/log_repeatReadSimulation.txt';	

	my $fh_log;
	open($fh_log, '>', $fn_log) or die "Cannot open $fn_log";
	
	
	my $MetaMap_taxonomy = {};
	my %contigID_2_taxonID;
	my %taxonID_2_contigIDs;	
	my %contig_2_length;
	fill_contigID_taxonID($DB, $MetaMap_taxonomy, \%contigID_2_taxonID, \%taxonID_2_contigIDs, \%contig_2_length);
	my %taxon_2_genomeLength = map {$_ => getGenomeLength($_, \%taxonID_2_contigIDs, \%contig_2_length)} keys %taxonID_2_contigIDs;
	
	my $n_simulations_file = $globalOutputDir . '/n_simulations.txt';
	open(N, '<', $n_simulations_file) or die "Cannot open $n_simulations_file";
	my $realizedN = <N>;
	chomp($realizedN);
	close(N);
	die unless($realizedN =~ /^\d+$/);
	for(my $simulationI = 0; $simulationI < $realizedN; $simulationI++)
	{
		my $simulation_href_fn = $globalOutputDir . '/' . $simulationI . '/simulationStore';
		my $simulation_href = retrieve $simulation_href_fn;
		
		setup_directory_and_simulate_reads($simulation_href, $MetaMap_taxonomy, $fh_log);
	}
	
	close($fh_log);
}
elsif($action eq 'prepareFromNonDB')
{
	my $DB_fa = $DB . '/DB.fa';
	my $existingSelfSimilarities = $DB . '/selfSimilarities.txt';
	
	die "Simulations in directory $globalOutputDir already prepared" if getFlag($globalOutputDir, 'simulated');
	unless($FASTA and (-e $FASTA))
	{
		die "Please provide valid --FASTA";
	}
	unless($FASTA_taxon_id)
	{
		die "Please provide valid --FASTA_taxon_id";
	}
	
	die unless(-e $DB_fa);
	die unless(-e $existingSelfSimilarities);
	$n_simulations = 1;
	
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
		print {$fh_log} "File-based simulation from file $FASTA\n";
		
		my $thisSimulation_outputDirectory = $globalOutputDir . '/' . 0;
		(mkdir($thisSimulation_outputDirectory) or die "Cannot mkdir $thisSimulation_outputDirectory") unless(-e $thisSimulation_outputDirectory);
		
			
		my $oneSimulation_href = returnOneSimulation_fromFASTA($DB, \%taxonID_2_contigIDs, \%taxon_2_genomeLength, $thisSimulation_outputDirectory, $coverageMode, $MetaMap_taxonomy, $fh_log, $FASTA, $FASTA_taxon_id);
		
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
	
	# check for required files
	my @warnings;
	for(my $jobI = 0; $jobI < $realizedN; $jobI++)  
	{
		my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
		my $simulation_href = retrieve $simulation_href_fn;
					
		for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++)
		{			
			my $DB_target_dir = $simulation_href->{outputDirectory} . '/DB_' . $simulation_href->{inferenceDBs}[$varietyI][2];
			die unless($DB_target_dir eq $simulation_href->{dbDirs_metamap}[$varietyI]);
		
			my $fn_selfSimilarities = $DB_target_dir . '/selfSimilarities.txt';
			push(@warnings, "File $fn_selfSimilarities missing") unless(-e $fn_selfSimilarities);
		}
	}
	die Dumper(@warnings) if(@warnings);
	
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
	unless($simulation_href->{outputDirectory} eq $globalOutputDir . '/' . $jobI)
	{
		die Dumper("I think you were too smart I", $simulation_href->{outputDirectory}, $globalOutputDir . '/' . $jobI);
	}
	foreach my $dir (@{$simulation_href->{dbDirs_metamap}}[0])
	{
		unless($dir eq $globalOutputDir . '/' . $jobI  . '/DB_fullDB')
		{
			die Dumper("I think you were too smart II", $simulation_href->{dbDirs_metamap}, $globalOutputDir . '/' . $jobI  . '/DB_fullDB');
		}	
	}
					
	inferenceOneSimulation($simulation_href, $skipKraken, $maxMemory);
}
elsif($action eq 'analyzeJobI')
{
	die"Please specify --jobI" unless(defined $jobI);
	my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
	my $simulation_href = retrieve $simulation_href_fn;			
	unless($simulation_href->{outputDirectory} eq $globalOutputDir . '/' . $jobI)
	{
		die Dumper("I think you were too smart", $simulation_href->{outputDirectory}, $globalOutputDir . '/' . $jobI);
	}
	unless($simulation_href->{dbDirs_metamap} eq $globalOutputDir . '/' . $jobI  . '/DB_fullDB')
	{
		die Dumper("I think you were too smart", $simulation_href->{dbDirs_metamap}, $globalOutputDir . '/' . $jobI  . '/DB_fullDB');
	}		
	evaluateOneSimulation($simulation_href);
}
elsif($action eq 'printVarietiesAndDirectories')
{
	my $n_simulations_file = $globalOutputDir . '/n_simulations.txt';
	open(N, '<', $n_simulations_file) or die "Cannot open $n_simulations_file";
	my $realizedN = <N>;
	chomp($realizedN);
	close(N);
	die unless($realizedN =~ /^\d+$/);
	
	my $haveAllFiles = 0;
	for(my $jobI = 0; $jobI < $realizedN; $jobI++)  
	{
		my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
		my $simulation_href = retrieve $simulation_href_fn;
		
		for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++)
		{
			my $varietyName = $simulation_href->{inferenceDBs}[$varietyI][2];
			print "Simulation $jobI - variety $varietyI - name $varietyName\n";
		}
	}
}
elsif($action eq 'analyzeAll')
{
	my $n_simulations_file = $globalOutputDir . '/n_simulations.txt';
	open(N, '<', $n_simulations_file) or die "Cannot open $n_simulations_file";
	my $realizedN = <N>;
	chomp($realizedN);
	close(N);
	die unless($realizedN =~ /^\d+$/);

	my $outputDir_allSimulations = $DB;
	
	# check that required files are present
	my $haveAllFiles = 1;
	for(my $jobI = 0; $jobI < $realizedN; $jobI++)  
	{
		my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
		my $simulation_href = retrieve $simulation_href_fn;
		
		#for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++) # todo!!!
		for(my $varietyI = 0; $varietyI <= 0; $varietyI++)
		{
			my $varietyName = $simulation_href->{inferenceDBs}[$varietyI][2];
			
			if($ignoreFullDB)
			{
				if($varietyName =~ /fullDB/)
				{
					#die "Skipping job $jobI $varietyName\n";
					next;
				}
			}
			# print "Check $varietyName\n";
			# die if($varietyName =~ /_1/); # we still need to deal with this and remove the numerical indices for storing

			my $simulation_results_dir = $simulation_href->{outputDirectory} . '/inference_' . $varietyName;
		
			my %expected_results_files = get_files_for_evaluation();

			foreach my $methodName (keys %expected_results_files)
			{
				my $methodDetails = $expected_results_files{$methodName};
				my $evaluationType = $methodDetails->[0];
				my $f = $simulation_results_dir . '/' . $methodDetails->[1];
				my $optional = $methodDetails->[2];
				if(-e $f)
				{ 
					if($methodName =~ /meta/i)
					{	
						my $fh ;
						open($fh, '<', $f) or die;
						my $timestamp = ctime(stat($fh)->mtime);
						# print "\t $f $timestamp\n";
					}				
				}
				else
				{

					if($optional)
					{
						warn "Expected file $f for $methodName not present -- ignore.";
					}
					else
					{
						warn "Expected file $f for $methodName not present.";	
						$haveAllFiles = 0;
					}
				}
			}
		}
	}
	
	unless($haveAllFiles)
	{
		die "Not all inference files present, abort." ;
	}
	
	my $allSimulations_data_href = validation::getEmptyGlobalResltsStore();
	$allSimulations_data_href->{realizedN} = $realizedN;
	
	# my @{$allSimulations_data_href->{n_reads_correct_byVariety_bySimulation}};
	# my @{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel_bySimulation}};
	# my @{$allSimulations_data_href->{freq_byVariety_byLevel_bySimulation}};
	# my @{$allSimulations_data_href->{directlyMappable_bySimulation}};
	
	# my %n_reads_correct_byVariety;
	# my %n_reads_correct_byVariety_byLevel;
	# my %freq_byVariety_byLevel;
	# my @frequencyComparisons_bySimulation;
	# my @frequencyComparisons_details_bySimulation;
	
	my $fullTaxonomy_simulation = taxTree::readTaxonomy($DB . '/taxonomy');
	
	# my @highLevel_stats_keptSeparate_bySimulation;
	# my %callRate_and_accuracy_byReadCategory;
	# my %callRate_and_accuracy_byReadCategory_byLength;
	# my %attachedTo_byReadCategory;
	
	
	
	# for(my $jobI = 0; $jobI < $realizedN; $jobI++)
	for(my $jobI = 0; $jobI < $realizedN; $jobI++)  
	{
		my %n_reads_correct_byVariety_byLevel_byLength;

		my $simulation_href_fn = $globalOutputDir . '/' . $jobI . '/simulationStore';
		my $simulation_href = retrieve $simulation_href_fn;
		my $frequencyComparison = {};
		my $frequencyComparison_details = {};
		my %n_reads_correct_byVariety_local;
		my %n_reads_correct_byVariety_byLevel_local;	
		my %n_reads_correct_byVariety_byLevel_byLength_local;		
		my %n_reads_unknownStats_byLevel;		
		my %freq_byVariety_byLevel_local;		
		my %freq_details_byVariety_byLevel_local;	
		my %directlyMappable;
		evaluateOneSimulation($simulation_href, \%n_reads_correct_byVariety_local, \%n_reads_correct_byVariety_byLevel_local, \%n_reads_correct_byVariety_byLevel_byLength_local, \%freq_byVariety_byLevel_local, $frequencyComparison, $frequencyComparison_details, \%directlyMappable, \%n_reads_unknownStats_byLevel);

		validation::addResultsToGlobalStore(
			$jobI,
			$allSimulations_data_href,
			$frequencyComparison,
			$frequencyComparison_details,
			\%n_reads_correct_byVariety_local,
			\%n_reads_correct_byVariety_byLevel_local,
			\%n_reads_correct_byVariety_byLevel_byLength_local,
			\%n_reads_unknownStats_byLevel,
			\%freq_byVariety_byLevel_local,
			\%freq_details_byVariety_byLevel_local,
			\%directlyMappable,
		);
		
	} ## end read data from one simulation

	## output summary stats

	validation::produceValidationOutputFiles($allSimulations_data_href, $fullTaxonomy_simulation, $globalOutputDir . '/', $outputDir_allSimulations, $suffix);
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
	my $skipKraken = shift;
	my $maxMemory = shift;

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
		 
		doMetaMap($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $maxMemory);
		unless($skipKraken)
		{
			SimulationsKraken::doKraken($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $krakenDBTemplate, $kraken_binPrefix, $Bracken_dir);
			SimulationsKraken::doKraken2($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $kraken2DBTemplate, $kraken2_binPrefix);			
			if($simulation_href->{inferenceDBs}[$varietyI][2] eq 'fullDB')
			{
				# SimulationsMetaPalette::doMetaPalette($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $metaPalette_installation_dir, $jellyfish_2_bin, $masterTaxonomy_dir);			
			}
		}
		
		SimulationsKraken::doCentrifuge($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $centrifugeBinDir); 
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
		# 'Bracken-Dist' => ['distribution', 'results_bracken.txt.ignoreUnclassified'],
		# 'Kraken-Dist' => ['distribution', 'results_kraken.txt.ignoreUnclassified'],
		'Bracken-Dist' => ['distribution', 'results_bracken.txt'],
		'Kraken-Dist' => ['distribution', 'results_kraken.txt'],		
		'MetaMap-EM-Dist' => ['distribution', 'metamap.EM.WIMP'],
		# 'MetaMap-U-Dist' => ['distribution', 'metamap.U.WIMP'],
		'Kraken-Reads' => ['reads', 'results_kraken.txt.reads2Taxon'],
		# 'Metamap-U-Reads' => ['reads', 'metamap.U.reads2Taxon'],
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
	my $n_reads_correct_byVariety_byLevel_byLength = shift;
	my $freq_byVariety_byLevel = shift;
	my $frequencyComparison_href = shift;
	my $unknown_and_frequencyContributions_href = shift;
	my $directlyMappable_href = shift;
	my $n_reads_unknownStats_byLevel_href = shift;
	
	die unless(defined $directlyMappable_href);
	die unless(defined $n_reads_unknownStats_byLevel_href);
	
	my $taxonomyFromSimulation_dir = $simulation_href->{outputDirectory} . '/DB_fullDB/taxonomy';
	my $taxonomy_usedForSimulation = taxTree::readTaxonomy($taxonomyFromSimulation_dir);
		
	(my $extendedMaster, my $extendedMaster_merged) = validation::prepare_masterTaxonomy_withX($masterTaxonomy_dir, $taxonomy_usedForSimulation);

	my $truth_fn = $simulation_href->{outputDirectory} . '/truth_reads.txt';
	my $truth_raw_reads_href = validation::readTruthFileReads($extendedMaster, $extendedMaster_merged, $truth_fn);
	my $readLengths_href = Util::getReadLengths($simulation_href->{outputDirectory} . '/reads.fastq');
	
	my %taxonIDs_in_direct_truth = map {$_ => 1} values %$truth_raw_reads_href;
	foreach my $taxonID (keys %taxonIDs_in_direct_truth)
	{
		my @descendants = taxTree::descendants($extendedMaster, $taxonID);
		foreach my $descendantID (@descendants)
		{
			die unless(exists $extendedMaster->{$descendantID});
			die if(exists $taxonIDs_in_direct_truth{$descendantID});
		}
	}
	# my %truth_raw_taxonIDs;
	# $truth_raw_taxonIDs{$taxonID_master}++;
	
	my %expected_results_files = get_files_for_evaluation(); 
	# for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++) # todo!!!!
	for(my $varietyI = 0; $varietyI <= 0; $varietyI++)
	{
		# last if($varietyI > 2);
		
		my $varietyName = $simulation_href->{inferenceDBs}[$varietyI][2];
		
		my @varietyNames_forStorage = ($varietyName);
		push(@varietyNames_forStorage, 'allCombined_' . $varietyI); 
		if($varietyName =~ /remove/)
		{
			push(@varietyNames_forStorage, 'incompleteCombined_' . $varietyI);
		}

		if($ignoreFullDB)
		{
			if($varietyName =~ /fullDB/)
			{
				#die "Skipping job $jobI $varietyName\n";
				next;
			}
		}		
		print "Analyse $varietyName\n";
		# if($varietyI >= 3)
		# {
			# warn "Fix this!!!!!!!!!!!!!!!!!!";
			# last;
		# }	
		# die if($varietyName =~ /_1/);
		# last if($varietyName =~ /_3/);
		
		foreach my $varietyName_forStorage (@varietyNames_forStorage)
		{
			$n_reads_correct_byVariety->{$varietyName_forStorage} = {} unless(defined $n_reads_correct_byVariety->{$varietyName_forStorage});
			$n_reads_correct_byVariety_byLevel->{$varietyName_forStorage} = {} unless(defined $n_reads_correct_byVariety_byLevel->{$varietyName_forStorage});
			$n_reads_correct_byVariety_byLevel_byLength->{$varietyName_forStorage} = {} unless(defined $n_reads_correct_byVariety_byLevel_byLength->{$varietyName_forStorage});
			$freq_byVariety_byLevel->{$varietyName_forStorage} = {} unless(defined $freq_byVariety_byLevel->{$varietyName_forStorage});
			$frequencyComparison_href->{$varietyName_forStorage} = {} unless(defined $frequencyComparison_href->{$varietyName_forStorage});
			$n_reads_unknownStats_byLevel_href->{$varietyName_forStorage} = {} unless(defined $n_reads_unknownStats_byLevel_href->{$varietyName_forStorage});
		}
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
		
		foreach my $tID (keys %reduced_taxonID_master_2_contigs)
		{
			$directlyMappable_href->{$varietyName}{$tID} = 1;
		}
		
		# reduce master taxonomy
		my $specificTaxonomy = dclone $extendedMaster;
		taxTree::removeUnmappableParts($specificTaxonomy, \%reduced_taxonID_master_2_contigs);
	
		# translate truth into reduced representation
		my $truth_mappingDatabase_reads = validation::translateReadsTruthToReducedTaxonomy($extendedMaster, $specificTaxonomy, $truth_raw_reads_href);

		# truth reads tree
		my $truth_reads_novelTree = validation::truthReadsTree($extendedMaster, $truth_mappingDatabase_reads, \%taxonIDs_in_direct_truth);
		
		# get distribution
		my $truth_mappingDatabase_distribution = validation::truthReadsToTruthSummary($specificTaxonomy, $truth_mappingDatabase_reads, \%reduced_taxonID_master_2_contigs);

		
		# print Dumper($truth_mappingDatabase_distribution);
		
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
					die Dumper("Missing some reads in truth file $truth_fn (inference file $f)", @missing_readIDs[0 .. 5]);
				}
				
				foreach my $varietyName_forStorage (@varietyNames_forStorage)
				{		
					print "$methodName $varietyName_forStorage reads\n";
					validation::readLevelComparison(
						$extendedMaster,
						$truth_raw_reads_href,
						$truth_mappingDatabase_reads,
						$inferred_reads,
						$methodName,
						$n_reads_correct_byVariety->{$varietyName_forStorage},
						$n_reads_correct_byVariety_byLevel->{$varietyName_forStorage},
						$n_reads_correct_byVariety_byLevel_byLength->{$varietyName_forStorage},
						\%reduced_taxonID_master_2_contigs,
						$readLengths_href,
						$n_reads_unknownStats_byLevel_href->{$varietyName_forStorage}
					);
					print "$methodName $varietyName_forStorage reads done\n";
				}
			}
			else 
			{
				print "Read $f\n";
				my $inferred_distribution = validation::readInferredDistribution($extendedMaster, $extendedMaster_merged, $f);
				print "Analyse $f\n";
				foreach my $varietyName_forStorage (@varietyNames_forStorage)
				{
					print "$methodName $varietyName_forStorage freqs\n";
					validation::distributionLevelComparison($extendedMaster, $truth_mappingDatabase_distribution, $inferred_distribution, $methodName, $freq_byVariety_byLevel->{$varietyName_forStorage}, $frequencyComparison_href->{$varietyName_forStorage});
					print "$methodName $varietyName_forStorage freqs done\n";					
				}
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
		
		$unknown_and_frequencyContributions_href->{$varietyName}{mappable} = {};
		foreach my $taxonID (keys %reduced_taxonID_original_2_contigs)
		{
			$unknown_and_frequencyContributions_href->{$varietyName}{mappable}{$taxonID} = 1;
		}
		$unknown_and_frequencyContributions_href->{$varietyName}{truthReadsNovelTree} = $truth_reads_novelTree;
		$unknown_and_frequencyContributions_href->{$varietyName}{nReads} = scalar(keys %$truth_raw_reads_href);		
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

sub returnOneSimulation_fromFASTA
{
	my $DB = shift;
	my $taxon_2_contig_href = shift;
	my $taxon_2_genomelength_href = shift;
	my $outputDirectory = shift;
	my $coverageMode = shift;
	my $taxonomy_href = shift;
	my $fh_log = shift;
	my $FASTA = shift;
	my $FASTA_taxon_id = shift;

	die unless(defined $taxonomy_href);
	
	die unless(($coverageMode eq 'equal'));
	die "Please provide FASTA argument" unless($FASTA and (-e $FASTA));
	die "Please provide taxon ID for file $FASTA" unless($FASTA_taxon_id);
	
	my @selectedTaxa = ($FASTA_taxon_id);
	my %selectedTaxa_frequencies = ($FASTA_taxon_id => 1);
	
	my %_selectedTaxa = map {$_ => 1} @selectedTaxa;
	die "Duplicate selected taxa for simulation" unless(scalar(keys %_selectedTaxa) == scalar(@selectedTaxa));
			
	my %targetTaxonIDs;
	if($coverageMode eq 'equal')
	{
		my $fr_taxon = 1/scalar(@selectedTaxa);
		%targetTaxonIDs = map {$_ => $fr_taxon} @selectedTaxa;
	} 
	else
	{
		die;
	}
	
	my $simulation_href = {
		coverageFactor => 20,
		outputDirectory => $outputDirectory,
		DB_simulation => $DB,
		targetTaxons => \%targetTaxonIDs,
		inferenceDBs => [['$SAME', undef, 'fullDB']],
		externalTaxonIDs => {
			$FASTA_taxon_id => {
				file => $FASTA
			}
		},
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
	my $maxMemory = shift;
	
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
	my $file_res_mapping = $metamap_output_dir . '/resources_mapping';
	my $file_res_classification = $metamap_output_dir . '/resources_classification';
	
	die unless(-e $DB.'/DB.fa');
	# my $cmd_map = qq(/usr/bin/time -v $metamap_bin mapDirectly --all --maxmemory 20 -r $DB/DB.fa -q $reads -m 2000 --pi 80 -o $file_mappings &> file_res_mapping);
	my $maxmemory_switch = '';
	if($maxMemory)
	{
		$maxmemory_switch = " --maxmemory $maxMemory "
	}
	my $cmd_map = qq(/usr/bin/time -v $metamap_bin mapDirectly --all ${maxmemory_switch}-r $DB/DB.fa -q $reads -m 2000 --pi 80 -o $file_mappings &> $file_res_mapping);

	print "Now executing: $cmd_map\n"; # todo remove?
	# $cmd_map = qq(/usr/bin/time -v $metamap_bin mapDirectly --all -r $DB/DB.fa -q $reads -m 2000 --pi 80 -o $file_mappings);
	system($cmd_map) and die "Cannot execute $cmd_map";
	die "No MetaMap mappings -- $file_mappings" unless(-e $file_mappings);
		
	my $minReads = 1000;
	if($DB =~ /miniSeq_100/)
	{
		warn "Assume that we're using a very small simulation DB, set --minReads to 1";
		$minReads = 1; 
	}	
	my $cmd_classify = qq(/usr/bin/time -v $metamap_bin classify --DB $DB --mappings $file_mappings --minreads $minReads &> $file_res_classification);
	$cmd_classify =    qq(/usr/bin/time -v $metamap_bin classify --DB $DB --mappings $file_mappings --minreads $minReads &> $file_res_classification); 
	system($cmd_classify) and die "Cannot execute $cmd_classify";
	

	my $resultsFile_EM = $file_mappings . '.EM.WIMP';
	my $resultsFile_EM_reads = $file_mappings . '.EM.reads2Taxon';
	my $resultsFile_U = $file_mappings . '.U.WIMP';
	my $resultsFile_U_reads = $file_mappings . '.U.reads2Taxon';
	
	die "No MetaMap EM classification -- $resultsFile_EM" unless(-e $resultsFile_EM);
	die "No MetaMap EM-reads classification -- $resultsFile_U" unless(-e $resultsFile_EM_reads);	
	#die "No MetaMap EM-U classification -- $resultsFile_U" unless(-e $resultsFile_U);
	#die "No MetaMap EM-U-reads classification -- $resultsFile_U" unless(-e $resultsFile_U_reads);
	
	copy($resultsFile_EM, $dir);
	copy($resultsFile_EM_reads, $dir);	
	#copy($resultsFile_U, $dir);
	#copy($resultsFile_U_reads, $dir);
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
	
	my $external_taxonIDs_href = (exists $simulation_href->{externalTaxonIDs}) ? $simulation_href->{externalTaxonIDs} : {};
	
	my @targetTaxonIDs = sort keys %{$simulation_href->{targetTaxons}};
	die unless(all {$taxonomy_simulation->{$_} or $external_taxonIDs_href->{$_}} @targetTaxonIDs );
	die unless(all {$taxonID_2_contigIDs_href->{$_} or $external_taxonIDs_href->{$_}} @targetTaxonIDs );
	
	my $haveExternalTaxonData = scalar(grep {$external_taxonIDs_href->{$_}} @targetTaxonIDs);
	
	my $fullMasterTaxonomy = $taxonomy_simulation;
	if($haveExternalTaxonData)
	{
		(my $extendedMaster, my $extendedMaster_merged) = validation::prepare_masterTaxonomy_withX($masterTaxonomy_dir, $taxonomy_simulation);
		$fullMasterTaxonomy = $extendedMaster;
	}
	die unless(all {$fullMasterTaxonomy->{$_}} @targetTaxonIDs );			
	
	# check that target coverages are normalized
	{
		my %coverages_by_level;
		print {$fh_log} "\tTarget coverage proportions:\n";
		my $_s_taxons = 0;
		foreach my $taxonID (@targetTaxonIDs)
		{
			$_s_taxons += $simulation_href->{targetTaxons}{$taxonID};
			print {$fh_log} "\t\t", $taxonID, ": ", $simulation_href->{targetTaxons}{$taxonID}, " (", taxTree::taxon_id_get_name($taxonID, $fullMasterTaxonomy), ")\n";
			{
				my $ancestors_href = taxTree::get_ancestors_with_specific_ranks($fullMasterTaxonomy, $taxonID, [taxTree::getRelevantRanks()]);
				foreach my $rank (keys %$ancestors_href)
				{
					$coverages_by_level{$rank}{$ancestors_href->{$rank}} += $simulation_href->{targetTaxons}{$taxonID};
				}
			}
		}
		
		{
			my $print_taxon_level = 'superkingdom';
			print {$fh_log} "\n\t\tBy $print_taxon_level:\n";
			foreach my $taxonID (sort keys %{$coverages_by_level{$print_taxon_level}})
			{	
				print {$fh_log}  "\t\t\t", $taxonID, " / ", taxTree::taxon_id_get_name($taxonID, $fullMasterTaxonomy), ": ", $coverages_by_level{$print_taxon_level}{$taxonID}, "\n";
			}
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
		my $thisTaxon_genomeLength = (exists $external_taxonIDs_href->{$taxonID}) ? getExternalGenomeLength($external_taxonIDs_href->{$taxonID}{file}) : getGenomeLength($taxonID, $taxonID_2_contigIDs_href, $contig_2_length_href);	
		$taxa_genome_lengths{$taxonID} = $thisTaxon_genomeLength;
	
		if(not $coverageTargetsAreOrganismAbundances)
		{
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
		else
		{
			$simulationRelativeCoverages{$taxonID} = $simulation_href->{targetTaxons}{$taxonID};
		}
	}

	if(1)
	{
		print {$fh_log} "\tFinal (potentially genome-size adjusted) coverages:\n";
		foreach my $taxonID (sort keys %simulationRelativeCoverages)
		{
			print {$fh_log}  "\t\t", $taxonID, ": ", sprintf("%.3f", $simulationRelativeCoverages{$taxonID}), " (length ", sprintf("%.3f", $taxa_genome_lengths{$taxonID}/(1024**2)), "M)\n";
		}	
		print {$fh_log} "\n";
	}
	
	# get estimated total megabases
	my $totalExpectedBases = 0;
	for(my $taxonI = 0; $taxonI <= $#targetTaxonIDs; $taxonI++)
	{
		my $taxonID = $targetTaxonIDs[$taxonI];
		my $thisTaxon_genomeLength = (exists $external_taxonIDs_href->{$taxonID}) ? getExternalGenomeLength($external_taxonIDs_href->{$taxonID}{file}) : getGenomeLength($taxonID, $taxonID_2_contigIDs_href, $contig_2_length_href);	
		my $relativeCoverage = $simulationRelativeCoverages{$taxonID};
		die unless(defined $relativeCoverage);
		$relativeCoverage *= $simulation_href->{coverageFactor};
		$totalExpectedBases += ($relativeCoverage * $thisTaxon_genomeLength);
	}
		
	my $local_coverageFactor = $simulation_href->{coverageFactor};
	
	if(defined $targetTotalSimulationInGigabytes)
	{
		my $expectGB = $totalExpectedBases / 1e9;
		print {$fh_log} "Coverage correction: --targetTotalSimulationInGigabytes in effect -- expect about ${expectGB} GB, and want $targetTotalSimulationInGigabytes GB.\n";
		my $factor = $expectGB / $targetTotalSimulationInGigabytes;
		my $new_local_coverageFactor = $local_coverageFactor * (1 / $factor);
		print {$fh_log} "Correct by factor $factor, i.e. set coverageFactor from $local_coverageFactor to $new_local_coverageFactor\n"; 
		$local_coverageFactor = $new_local_coverageFactor;
	}
		
	my $genome_files_href = extractTaxonIDsFromREF($simulation_FASTA, $simulation_href->{outputDirectory} . '/forSimulation', \%simulationRelativeCoverages, $taxonID_2_contigIDs_href, $external_taxonIDs_href);	

	my $totalReads_allTaxa = 0;
	my $totalBases_allTaxa = 0;
	my %reads_taxon; 
	my %readID_2_taxon;
	my %taxonID_2_bases;
	for(my $taxonI = 0; $taxonI <= $#targetTaxonIDs; $taxonI++)
	{
		my $taxonID = $targetTaxonIDs[$taxonI];
		
		my $fn_genome = $genome_files_href->{$taxonID};
		my $relativeCoverage = $simulationRelativeCoverages{$taxonID};
		die unless((defined $fn_genome) and (defined $relativeCoverage));
		$relativeCoverage *= $local_coverageFactor;
		
		my $outputDir_reads = $fn_genome . '.reads';
		my $cmd_rm_outputDir_reads = "rm $outputDir_reads/*; rm -r $outputDir_reads";
		
		if(-e $outputDir_reads)
		{
			system($cmd_rm_outputDir_reads) and die "Couldn't execute: $cmd_rm_outputDir_reads";
		}
		
		mkdir($outputDir_reads) or die "Cannot mkdir $outputDir_reads";
		
		my $PBsim_thisSimulation = $PBsim_cmd;
		if($simulationMinReadLength)
		{
			die unless($PBsim_thisSimulation =~ /--data-type CLR/);
			$PBsim_thisSimulation =~ s/--data-type CLR/--data-type CLR --length-min $simulationMinReadLength /;
		}
		
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
		my $reads_taxon = combineFASTQ(\@files_reads, $outputFile_combinedReads, $taxonI, $doAppend, $readIDs_aref, \$combinedReadLength, $simulationMinReadLength);
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
		$totalBases_allTaxa += $combinedReadLength;
	}
	
	my %realizedTaxonProportions_reads;
	my %realizedOrganismAbundances;
	my $sum_organismAbundances = 0;
	foreach my $taxonID (@targetTaxonIDs)
	{
		$realizedTaxonProportions_reads{$taxonID} = $reads_taxon{$taxonID} / $totalReads_allTaxa;
		$realizedOrganismAbundances{$taxonID} = $taxonID_2_bases{$taxonID} / $taxa_genome_lengths{$taxonID};
		$sum_organismAbundances += $realizedOrganismAbundances{$taxonID};
	}
	my %realizedRelativeOrganismAbundances;	
	foreach my $taxonID (@targetTaxonIDs)
	{
		$realizedRelativeOrganismAbundances{$taxonID} = $realizedOrganismAbundances{$taxonID} / $sum_organismAbundances;
	}
	
	if(1)
	{
		foreach my $taxonID (@targetTaxonIDs)
		{
			my $bases = 0;
		}
		print {$fh_log} "\tRealized read proportions:\n";
		foreach my $taxonID (@targetTaxonIDs)
		{
			print {$fh_log} "\t\t", $taxonID, ": ", $reads_taxon{$taxonID}, " reads, read proportion ", sprintf("%.3f", $realizedTaxonProportions_reads{$taxonID}), "; copies of organism $realizedOrganismAbundances{$taxonID} / relative: $realizedRelativeOrganismAbundances{$taxonID}; vs target ", sprintf("%.3f", $simulation_href->{targetTaxons}{$taxonID}), "\n";
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
	
	simulation::truthReadFrequenciesFromReadCounts($outputFile_truth_readFrequencies_overCompleteTaxonomy, \%reads_taxon, $fullMasterTaxonomy);
	simulation::truthGenomeFrequenciesFromReadCounts($outputFile_truth_genomeFrequencies, \%taxonID_2_bases, \%reads_taxon, \%taxa_genome_lengths, $fullMasterTaxonomy);
	
	$simulation_href->{readsFastq} = $outputFile_combinedReads;
}

sub getExternalGenomeLength
{
	my $file = shift;
	my $genome_href = Util::readFASTA($file);
	my $L = 0;
	foreach my $contig (values %$genome_href)
	{
		$L += length($contig);
	}
	return $L;
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
	my $simulationMinReadLength = shift;
	
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
			
			if((not $simulationMinReadLength) or (length($line_seq) >= $simulationMinReadLength))
			{
				print OUT $line, "\n", $line_seq, "\n", $line_plus, $line_qual;
			}
			
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
	my $external_taxonIDs_href = shift;
	
	my %forReturn_files_by_taxonID;
	my %fh_per_taxonID;
	
	my %target_contig_2_taxonID;
	die unless(scalar(keys %$extractTaxonIDs_href));
	my @taxonIDs = sort keys %$extractTaxonIDs_href;
	my $externalGenome = 0;
	for(my $i = 0; $i <= $#taxonIDs; $i++)
	{
		my $taxonID = $taxonIDs[$i];
		die unless((exists $external_taxonIDs_href->{$taxonID}) xor (exists $taxon_2_contig_href->{$taxonID}));
		
		my $fnOut = $fn_out_prefix . '.' . $i;
		my $fh;
		open($fh, '>', $fnOut) or die "Cannot open $fnOut";
		
		$fh_per_taxonID{$taxonID} = $fh;
		$forReturn_files_by_taxonID{$taxonID} = $fnOut;
		
		# die unless(defined $taxon_2_contig_href->{$taxonID});
		if(exists $taxon_2_contig_href->{$taxonID})
		{
			foreach my $contigID (@{$taxon_2_contig_href->{$taxonID}})
			{
				$target_contig_2_taxonID{$contigID} = $taxonID;
			}
		}
		else
		{
			$externalGenome++;
			die unless($external_taxonIDs_href->{$taxonID}{file});
			open(COMPLETEGENOME, '<', $external_taxonIDs_href->{$taxonID}{file}) or die "Cannot open $external_taxonIDs_href->{$taxonID}{file}";
			while(<COMPLETEGENOME>)
			{
				print {$fh_per_taxonID{$taxonID}} $_;
			}
			close(COMPLETEGENOME);
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
	
	print "Printed $printedContigs from $fn_in --> ${fn_out_prefix}.* (plus $externalGenome external genomes)\n";
	
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
