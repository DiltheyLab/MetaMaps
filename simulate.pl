#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use List::Util qw/all any shuffle/;
use List::MoreUtils qw/mesh/;
use Getopt::Long;   
use File::Path qw(make_path remove_tree);
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;
use File::Copy;
use Storable qw/dclone store retrieve/;

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

Util::get_metaMap_bin_and_enforce_mainDir();
die unless(-e 'estimateSelfSimilarity.pl');

#my $create_simulations = 2;
#my $perSimulation_totalKnownGenomes = 5;
my $globalOutputDir = 'new_simulations/';
unless(-e $globalOutputDir)
{
	mkdir($globalOutputDir) or die "Cannot mkdir $globalOutputDir";
}
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

my $metamap_bin = './mashmap';
my $kraken_binPrefix = qq(/data/projects/phillippy/software/kraken-0.10.5-beta/bin/kraken);
my $Bracken_dir = qq(/data/projects/phillippy/software/Bracken/);
my $metaPalette_installation_dir = qq(/data/projects/phillippy/software/MetaPalette/);
my $jellyfish_2_bin = qq(/data/projects/phillippy/software/jellyfish-2.2.6/bin/jellyfish);
my $krakenDBTemplate = '/data/projects/phillippy/projects/mashsim/src/krakenDBTemplate/'; # make sure this is current!

my @evaluateAccuracyAtLevels = qw/species genus family superkingdom/;
{
	my %_knowRank = map {$_ => 1} taxTree::getRelevantRanks();
	die unless(all {$_knowRank{$_}} @evaluateAccuracyAtLevels);
}
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
GetOptions (
	'action:s' => \$action,
	'jobI:s' => \$jobI,
	'really:s' => \$really,
	'jobIMethod:s' => \$jobIMethod,
);

if(not defined $action)
{
	die "Please specify --action, e.g. prepare, inferenceJobI, analyzeJobI";
}

if($action eq 'prepare')
{
	my $DB_fa = $DB . '/DB.fa';
	my $existingSelfSimilarities = $DB . '/selfSimilarities.txt';
	
	die unless(-e $DB_fa);
	die unless(-e $existingSelfSimilarities);
	
	my $MetaMap_taxonomy = {};
	my %contigID_2_taxonID;
	my %taxonID_2_contigIDs;	
	my %contig_2_length;
	fill_contigID_taxonID($DB, $MetaMap_taxonomy, \%contigID_2_taxonID, \%taxonID_2_contigIDs, \%contig_2_length);
	my %taxon_2_genomeLength = map {$_ => getGenomeLength($_, \%taxonID_2_contigIDs, \%contig_2_length)} keys %taxonID_2_contigIDs;
		
	my @simulations_to_execute;
	my $n_simulations = 1;
	for(my $simulationI = 0; $simulationI < $n_simulations; $simulationI++)
	{
		my $thisSimulation_outputDirectory = $globalOutputDir . '/' . $simulationI;
		(mkdir($thisSimulation_outputDirectory) or die "Cannot mkdir $thisSimulation_outputDirectory") unless(-e $thisSimulation_outputDirectory);
		
		my $oneSimulation_href;
		if($simulationI < $n_simulations)
		{
			$oneSimulation_href = returnOneSimulation_noUnknown($DB, \%taxonID_2_contigIDs, \%taxon_2_genomeLength, 4, $thisSimulation_outputDirectory, 1, $MetaMap_taxonomy);
		}
		
		addInferenceRoundsWithReducedDBs($oneSimulation_href, $MetaMap_taxonomy);
		
		push(@simulations_to_execute, $oneSimulation_href);
	}
		
	foreach my $oneSimulation_href (@simulations_to_execute)
	{
		executeSimulation($oneSimulation_href, \%taxonID_2_contigIDs, \%contig_2_length);
	}
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
else
{
	die "Unknown --action";
}

sub inferenceOneSimulation
{
	my $simulation_href = shift;
	
	my $fullTaxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
	
	for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++)
	{
		print "Inference $simulation_href->{inferenceDBs}[$varietyI][2] \n";
		
		my $DB_target_dir = $simulation_href->{outputDirectory} . '/DB_' . $simulation_href->{inferenceDBs}[$varietyI][2];
		die unless($DB_target_dir eq $simulation_href->{dbDirs_metamap}[$varietyI]);
		
		my $inference_target_dir = $simulation_href->{outputDirectory} . '/inference_' . $simulation_href->{inferenceDBs}[$varietyI][2];
		
		(mkdir($inference_target_dir) or die "Cannot mkdir $inference_target_dir") unless(-d $inference_target_dir);
		
		print "Doing inference in $DB_target_dir\n";
		# doMetaMap($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq});
		SimulationsKraken::doKraken($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $krakenDBTemplate, $kraken_binPrefix, $Bracken_dir);
		if($simulation_href->{inferenceDBs}[$varietyI][2] eq 'fullDB')
		{
			SimulationsMetaPalette::doMetaPalette($inference_target_dir, $DB_target_dir, $simulation_href->{readsFastq}, $metaPalette_installation_dir, $jellyfish_2_bin, $fullTaxonomy);			
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
	return ('Bracken' => 'results_bracken.txt', 'Kraken' => 'results_kraken.txt', 'Metamap-U' => 'metamap.U.WIMP', 'Metamap-EM' => 'metamap.EM.WIMP');
	# return ('Metamap-U' => 'metamap.U.WIMP', 'Metamap-EM' => 'metamap.EM.WIMP');
}

sub evaluateOneSimulation
{
	my $simulation_href = shift;
	
	my $masterTaxonomy_merged;

	unless(defined $masterTaxonomy)
	{
		$masterTaxonomy_merged = taxTree::readMerged($masterTaxonomy_dir);
		$masterTaxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
	}
	
	my $taxonomyFromSimulation_dir = $simulation_href->{outputDirectory} . '/DB_fullDB/taxonomy';
	my $taxonomy_usedForSimulation = taxTree::readTaxonomy($taxonomyFromSimulation_dir);
		
	my $extendedMaster = taxTree::cloneTaxonomy_integrateX($masterTaxonomy, $masterTaxonomy_merged, $taxonomy_usedForSimulation);
	
	my %truth_raw_reads;
	my %truth_raw_taxonIDs;
	my $truth_fn =  $simulation_href->{outputDirectory} . '/truth_reads.txt';
	open(T, '<', $truth_fn) or die "Cannot open $truth_fn";
	while(<T>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		my $readID = $f[0];
		my $taxonID_original = $f[1];
		
		# get the taxon ID in the master taxonomy
		my $taxonID_master = taxTree::findCurrentNodeID($extendedMaster, $masterTaxonomy_merged, $taxonID_original);
		
		die unless(exists $extendedMaster->{$taxonID_master});
		die if(defined $truth_raw_reads{$readID});
		
		$truth_raw_reads{$readID} = $taxonID_master;
		$truth_raw_taxonIDs{$taxonID_master}++;
	}
	close(T);
			
	my %expected_results_files = get_files_for_evaluation();
	for(my $varietyI = 0; $varietyI <= $#{$simulation_href->{dbDirs_metamap}}; $varietyI++)
	{
		my $varietyName = $simulation_href->{inferenceDBs}[$varietyI][2];
		print "Analyse $varietyName\n";
		
		my $simulation_results_dir = $simulation_href->{outputDirectory} . '/inference_' . $varietyName;
		my $DBdir = $simulation_href->{outputDirectory} . '/DB_' . $varietyName;
		
		# my $specificTaxonomy = taxTree::readTaxonomy($DBdir . '/taxonomy');
		
		# details for reduced DB
		my %reduced_taxonID_original_2_contigs;
		my %reduced_contigLength;
		Util::read_taxonIDs_and_contigs($DBdir, \%reduced_taxonID_original_2_contigs, \%reduced_contigLength);
		
		# read reduced taxonomy
		my %reduced_taxonID_master_2_contigs;
		foreach my $taxonID_original (keys %reduced_taxonID_original_2_contigs)
		{
			my $taxonID_master = taxTree::findCurrentNodeID($extendedMaster, $masterTaxonomy_merged, $taxonID_original);
			# store which taxon IDs *are* mappable
			$reduced_taxonID_master_2_contigs{$taxonID_master} = $reduced_taxonID_original_2_contigs{$taxonID_original};
		}
		
		# reduce master taxonomy
		my $specificTaxonomy = dclone $extendedMaster;
		taxTree::removeUnmappableParts($specificTaxonomy, \%reduced_taxonID_master_2_contigs);
	
		# translate truth
		my %truth_allReads;
		my %taxonID_translation;
		foreach my $taxonID (keys %truth_raw_taxonIDs)
		{
			$taxonID_translation{$taxonID} = getRank_phylogeny($extendedMaster, $specificTaxonomy, $taxonID);
			my $count = $truth_raw_taxonIDs{$taxonID};
			foreach my $rank (keys %{$taxonID_translation{$taxonID}})
			{
				my $v = $taxonID_translation{$taxonID}{$rank};
				$truth_allReads{$rank}{$v} += $count;
			}
		}
		foreach my $rank (keys %truth_allReads)
		{
			foreach my $taxonID (keys %{$truth_allReads{$rank}})
			{
				$truth_allReads{$rank}{$taxonID} /= scalar(keys %truth_raw_reads);
			}
		}
		
		foreach my $name (keys %expected_results_files)
		{
			my $f = $simulation_results_dir . '/' . $expected_results_files{$name};
			unless(-e $f)
			{
				warn "Expected file $f for $name not present.";
				next;
			}
			
			# unclassified = deliberately no caller
			# undefined = there would be an assignment, but there is no node
			my %inference;
			open(I, '<', $f) or die "Cannot open $f";
			my $header_line = <I>;
			chomp($header_line);
			my @header_fields = split(/\t/, $header_line);
			while(<I>)
			{
				my $line = $_;
				chomp($line);
				next unless($line);
				my @line_fields = split(/\t/, $line, -1);
				die unless($#line_fields == $#header_fields);
				my %line = (mesh @header_fields, @line_fields);
				
				my $taxonID_nonMaster = ($line{ID} // $line{taxonID});
				die unless(defined $taxonID_nonMaster);
				
				next if($line{Name} eq 'TooShort');
				next if($line{Name} eq 'Unmapped');

				if(($taxonID_nonMaster eq '0') and (($line{Name} eq 'Undefined') or ($line{Name} eq 'Unclassified')))
				{
					$taxonID_nonMaster = $line{Name};
				}
				if($taxonID_nonMaster eq '0')
				{
					die Dumper($line, $f);
				}
				
				my $taxonID_master = taxTree::findCurrentNodeID($extendedMaster, $masterTaxonomy_merged, $taxonID_nonMaster);
							
				die Dumper(\%line) unless(defined $taxonID_master);
				die Dumper("Unknown taxon ID $taxonID_master in file $f", $extendedMaster->{$taxonID_master}) unless(($taxonID_master eq 'Unclassified') or (defined $specificTaxonomy->{$taxonID_master}));
				die unless(defined $line{Absolute});
				die unless(defined $line{PotFrequency});
				$inference{$line{AnalysisLevel}}{$taxonID_master}[0] += $line{Absolute};
				$inference{$line{AnalysisLevel}}{$taxonID_master}[1] += $line{PotFrequency};
			}
			close(I);
			
			foreach my $level (@evaluateAccuracyAtLevels)
			{
				next unless(defined $inference{$level});
				die unless(defined $truth_allReads{$level});
				my $totalFreq = 0;
				my $totalFreqCorrect = 0;
				foreach my $inferredTaxonID (keys %{$inference{$level}})
				{
					my $isFreq = $inference{$level}{$inferredTaxonID}[1];
					my $shouldBeFreq = 0;
					if(exists $truth_allReads{$level}{$inferredTaxonID})
					{
						$shouldBeFreq = $truth_allReads{$level}{$inferredTaxonID};
					}
					
					$totalFreq += $isFreq;
					if($isFreq <= $shouldBeFreq)					
					{
						$totalFreqCorrect += $isFreq;
					}
					else
					{
						$totalFreqCorrect += $shouldBeFreq;
					}
				}
				die Dumper("Weird total freq", $name, $totalFreq, $f, $level) unless(abs(1 - $totalFreq) <= 1e-3);
				
				print join("\t", $varietyName, $name, $level, $totalFreqCorrect), "\n";
			}
		}
			
	}
}

sub getRank_phylogeny
{
	my $fullTaxonomy = shift;
	my $reducedTaxonomy = shift;
	my $taxonID = shift;
	
	my %forReturn = map {$_ => 'Unclassified'} @evaluateAccuracyAtLevels;
	
	unless(all {exists $fullTaxonomy->{$_}} keys %$reducedTaxonomy)
	{
		my @missing = grep {not exists $fullTaxonomy->{$_}} keys %$reducedTaxonomy;
		die Dumper("Reduced taxonomy doesn't seem to be a proper subset of full taxonomy!", \@missing);
	}
	die unless(defined $fullTaxonomy->{$taxonID});
	
	my @nodes_to_consider = ($taxonID, taxTree::get_ancestors($fullTaxonomy, $taxonID));
	
	my $inTaxonomicAgreement = 0;
	my $firstRankAssigned;
	foreach my $nodeID (@nodes_to_consider)
	{
		my $rank = $fullTaxonomy->{$nodeID}{rank};
		die unless(defined $rank);
		if(exists $forReturn{$rank})
		{
			if(exists $reducedTaxonomy->{$nodeID})
			{
				$forReturn{$rank} = $nodeID;
				$inTaxonomicAgreement = 1;
				$firstRankAssigned = $rank;
			}
			else
			{
				$forReturn{$rank} = 'Unclassified';
				die if($inTaxonomicAgreement);
			}
		}
	}
	
	if(defined $firstRankAssigned)
	{
		my $setToUndefined = 0;
		foreach my $rank (taxTree::getRelevantRanks())
		{
			next unless(exists $forReturn{$rank});
			if($rank eq $firstRankAssigned)
			{
				$setToUndefined = 1;
				die if($forReturn{$rank} eq 'Unclassified');
			}
			else
			{
				if($setToUndefined and ($forReturn{$rank} eq 'Unclassified'))
				{
					$forReturn{$rank} = 'Undefined';
				}
			}
		}
	}
	
	return \%forReturn;
}


sub addInferenceRoundsWithReducedDBs
{
	my $simulation_href = shift;
	my $taxonomy_href = shift;

	my @targetTaxons_shuffled = shuffle keys %{$simulation_href->{targetTaxons}};
	die unless(scalar(@targetTaxons_shuffled) > 1);

	my $removal_origin = $targetTaxons_shuffled[0];
	$removal_origin = 1496303; # todo remove
	
	my $ancestors_href = taxTree::get_ancestors_by_rank($taxonomy_href, $removal_origin);
	$ancestors_href->{self} = $removal_origin;
	
	foreach my $removeRank (qw/self species genus/)
	{
		next unless(exists $ancestors_href->{$removeRank});
		my $removeNode = $ancestors_href->{$removeRank};
		push(@{$simulation_href->{inferenceDBs}}, ['remove', [$removeNode], 'removeOne_' . $removeRank]);
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
	die unless(defined $taxonomy_href);
	
	die unless($n_species <= scalar(keys %$taxon_2_contig_href));
	
	my @availableTaxonIDs = shuffle(keys %$taxon_2_contig_href);
	
	@availableTaxonIDs = grep {die unless(defined $taxon_2_genomelength_href->{$_}); ($taxon_2_genomelength_href->{$_} > (5 * $simulation_read_length))} @availableTaxonIDs;
	
	my @selectedTaxa = @availableTaxonIDs[0 .. ($n_species-1)];
	
	my %targetTaxonIDs;
	if($equalCoverage)
	{
		my $fr_taxon = 1/scalar(@selectedTaxa);
		%targetTaxonIDs = map {$_ => $fr_taxon} @selectedTaxa;
	}
	else
	{
		die "Not implemented yet";
	}
	
	my $simulation_href = {
		coverageFactor => 1,
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
	
	setup_directory_and_simulate_reads($simulation_href);
	prepareDBs($simulation_href);
	
	die unless($simulation_href->{dbDirs_metamap});
	die unless($simulation_href->{readsFastq});
	
	store $simulation_href, $simulation_href->{outputDirectory} . '/simulationStore';
	
	inferenceOneSimulation($simulation_href);

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
		
	# todo perhaps reconsider minreads
	my $cmd_classify = qq(/usr/bin/time -v $metamap_bin classify --DB $DB --mappings $file_mappings --minreads 1 &> file_res_classification);
	$cmd_classify = qq(/usr/bin/time -v $metamap_bin classify --DB $DB --mappings $file_mappings --minreads 1); 
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
	
	foreach my $inference_variety (@{$simulation_href->{inferenceDBs}})
	{
		my $DB_target_dir = $simulation_href->{outputDirectory} . '/DB_' . $inference_variety->[2];
		my $outputFile_truth_inferenceDB = $DB_target_dir . '/truth_inferenceDB.txt';
		
		if($inference_variety->[0] eq '$SAME')
		{
			unless(-d $DB_target_dir)
			{
				mkdir($DB_target_dir) or die "Cannot mkdir $DB_target_dir";
			}
			
			Util::copyMetaMapDB($simulation_href->{DB_simulation}, $DB_target_dir);
			
			# todo perhaps resoncier
			# truthFileFromReadCounts($outputFile_truth_inferenceDB, \%reads_taxon, $taxonomy_simulation);
		}
		else
		{
			die unless($inference_variety->[0] eq 'remove');
			produceReducedDB($simulation_href->{DB_simulation}, $inference_variety->[1], $DB_target_dir);	
			
			# perhaps reconsider
			# truthFileForReducedInferenceDB();
		}
		
		push(@{$simulation_href->{dbDirs_metamap}}, $DB_target_dir);	
	}
}

sub setup_directory_and_simulate_reads
{
	my $simulation_href = shift;

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
	my $outputFile_truthOverCompleteTaxonomy = $simulation_href->{outputDirectory} . '/truth_frequencies_completeTaxonomy.txt';
	
	my %contigs_and_relative_coverage;
	
	my @targetTaxonIDs = sort keys %{$simulation_href->{targetTaxons}};
	die unless(all {$taxonomy_simulation->{$_}} @targetTaxonIDs );
	die unless(all {$taxonID_2_contigIDs_href->{$_}} @targetTaxonIDs );
	
	# check that target coverages are normalized
	{
		print "Target coverage proportions:\n" if($verbose);
		my $_s_taxons = 0;
		foreach my $taxonID (@targetTaxonIDs)
		{
			$_s_taxons += $simulation_href->{targetTaxons}{$taxonID};
			print "\t", $taxonID, ": ", $simulation_href->{targetTaxons}{$taxonID}, "\n" if ($verbose);
		}
		die unless(abs(1 - $_s_taxons) <= 1e-3);
		print "\n" if ($verbose);
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

	if($verbose)
	{
		print "Genome-size adjusted coverages:\n";
		foreach my $taxonID (keys %simulationRelativeCoverages)
		{
			print "\t", $taxonID, ": ", sprintf("%.3f", $simulationRelativeCoverages{$taxonID}), " (length ", sprintf("%.3f", $taxa_genome_lengths{$taxonID}/(1024**2)), "M)\n";
		}	
		print "\n";
	}
	
	my $genome_files_href = extractTaxonIDsFromREF($simulation_FASTA, $simulation_href->{outputDirectory} . '/forSimulation', \%simulationRelativeCoverages, $taxonID_2_contigIDs_href);
	
	my $totalReads_allTaxa = 0;
	my %reads_taxon;
	my %readID_2_taxon;
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
		my $reads_taxon = combineFASTQ(\@files_reads, $outputFile_combinedReads, $taxonI, $doAppend, $readIDs_aref);
		$reads_taxon{$taxonID} = $reads_taxon;
		$totalReads_allTaxa += $reads_taxon;
		foreach my $readID (@$readIDs_aref)
		{
			die if(defined $readID_2_taxon{$readID});
			$readID_2_taxon{$readID} = $taxonID;
		}
		system($cmd_rm_outputDir_reads) and die "Couldn't execute: $cmd_rm_outputDir_reads";
		
		unlink($fn_genome) or die "Cannot delete $fn_genome";
	}
	
	my %realizedTaxonProportions;
	foreach my $taxonID (@targetTaxonIDs)
	{
		$realizedTaxonProportions{$taxonID} = $reads_taxon{$taxonID} / $totalReads_allTaxa;
	}
	
	if($verbose)
	{
		print "Realized coverages:\n";
		foreach my $taxonID (@targetTaxonIDs)
		{
			print "\t", $taxonID, ": ", $reads_taxon{$taxonID}, " reads, proportion ", sprintf("%.3f", $realizedTaxonProportions{$taxonID}), "; vs target ", sprintf("%.3f", $simulation_href->{targetTaxons}{$taxonID}), "\n";
		}	
		print "\n";
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
	
	truthFileFromReadCounts($outputFile_truthOverCompleteTaxonomy, \%reads_taxon, $taxonomy_simulation);
	
	$simulation_href->{readsFastq} = $outputFile_combinedReads;
}

sub produceReducedDB
{
	my $baseDB = shift;
	my $removeNodes_aref = shift;
	my $targetDir = shift;
	
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
	warn  "Check command:\n$cmd_self_similarity\n\n";
	
	# todo activate
	# system($cmd_self_similarity) and die "Command $cmd_self_similarity failed";

	
	print "\t\tCreated reduced DB ", $targetDir, " with ", scalar(keys %$reducedTaxonomy), " nodes instead of ", scalar(keys %$taxonomy_base), " nodes\n";

}

sub truthFileForReducedInferenceDB
{
	die;
}

sub truthFileFromReadCounts
{
	my $outputFn = shift;
	my $readCounts_href = shift;
	my $taxonomy_href = shift;
	
	my %relevantLevels = map {$_ => 1} taxTree::getRelevantRanks();
	
	my %readCount_by_level;
	foreach my $taxonID (keys %$readCounts_href)
	{
		my %thisTaxon_levels;
		$thisTaxon_levels{'EqualCoverageUnit'} = $taxonID;
	
		my @relevantNodes = ($taxonID, taxTree::get_ancestors($taxonomy_href, $taxonID));
		
		
		foreach my $n (@relevantNodes)
		{
			my $rank = $taxonomy_href->{$n}{rank};
			die unless(defined $rank);
			if(exists $relevantLevels{$rank})
			{
				die if(defined $thisTaxon_levels{$rank});
				$thisTaxon_levels{$rank} = $n;
			}
		}
		
		foreach my $l ('EqualCoverageUnit', keys %relevantLevels)
		{
			if(not defined $thisTaxon_levels{$l})
			{
				$thisTaxon_levels{$l} = 'Undefined';
			}
			
			$readCount_by_level{$l}{$thisTaxon_levels{$l}} += $readCounts_href->{$taxonID};
		}
	}
	
	my $totalReads;
	open(O, '>', $outputFn) or die "Cannot open $outputFn";
	print O join("\t", qw/AnalysisLevel taxonID Name Absolute PotFrequency/), "\n";
	foreach my $level (keys %readCount_by_level)
	{
		my $totalReads_thisLevel = 0;
		foreach my $taxonID (keys %{$readCount_by_level{$level}})
		{
			$totalReads_thisLevel += $readCount_by_level{$level}{$taxonID};
		}
		
		if(not defined $totalReads)
		{
			$totalReads = $totalReads_thisLevel;
		}
		else
		{
			die "Discrepancy with total read counts: $totalReads vs $totalReads_thisLevel at level $level" unless($totalReads == $totalReads_thisLevel);
		}
		
		foreach my $taxonID (keys %{$readCount_by_level{$level}})
		{
			my $nReads = $readCount_by_level{$level}{$taxonID};
			my $f = $nReads / $totalReads;
			my $name = ($taxonID ne 'Undefined') ? taxTree::taxon_id_get_name($taxonID, $taxonomy_href) : 'Undefined';
			print O join("\t", $level, $taxonID, $name, $nReads, $f), "\n";
		}
		
	}
	close(O);
	
}
sub combineFASTQ
{
	my $in_aref = shift;
	my $out_fn = shift;
	my $readID_prefix = shift;
	my $append = shift;
	my $readIDs_aref = shift;
	
	my $totalReads = 0;
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
			my $line_plus = <F>;
			my $line_qual = <F>;
			die unless(substr($line_plus, 0, 1) eq '+');
			
			substr($line, 0, 1) = ('@' . $readID_prefix);
			
			print OUT $line, "\n", $line_seq, $line_plus, $line_qual;
			
			if(defined($readIDs_aref))
			{
				push(@$readIDs_aref, substr($line, 1));
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
