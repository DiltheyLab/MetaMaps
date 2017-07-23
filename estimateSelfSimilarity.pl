#!/usr/bin/perl -w

# Example command:
# ./estimateWithinTreeDistances.pl.pl --action dbDir ../db/refseq --taxonomyDir /data/projects/phillippy/projects/mashsim/NCBI/refseq/taxonomy/

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

$| = 1;

use taxTree;
use Util;

my @taxonomy_fields = qw/species genus family order phylum superkingdom/;
my $metamap_bin = './mashmap';
unless(-e $metamap_bin)
{
	die "Please execute me from the main MetaMap directory";
}

my $action = '';   

my $DB = '';
my $templateDB;
my $mode;
my $jobI;
my $jobITemplate;
my $readSimSize = 2000;
my $readSimDelta = 1000;
my $readSimSizeFrom = 2000;
my $readSimSizeTo = 50000;
my $readSimSizeStep = 1000;

my $target_max_simulatedChunks = 2000;

my $path_to_myself = $FindBin::Bin . '/estimateSelfSimilarity.pl';
my $path_to_cwd = getcwd();

GetOptions (
	'DB:s' => \$DB, 
	'templateDB:s' => \$templateDB, 
	'jobITemplate:s' => \$jobITemplate, 
	'mode:s' => \$mode, 
	'jobI:s' => \$jobI, 
);

unless($DB and (-d $DB))
{
	die "Please specify valid DB directory via --DB";
}
my $taxonomyDir = $DB . '/taxonomy';

my $file_ref = $DB . '/DB.fa';
die "File $file_ref not existing" unless(-e $file_ref);

my $outputDir = $DB . '/selfSimilarity';
my $outputDir_computation = $outputDir . '/computation';
my $outputDir_results = $outputDir . '/results';
my $outputFn_jobs = $outputDir . '/jobs';
my $outputFn_jobs_fromTemplate = $outputDir . '/jobs.fromTemplate';
my $outputFn_qsub = $outputDir . '/jobs.qsub';
my $outputFn_reads_results_many = $outputDir . '/results.reads.many.byNode';
my $finalResultsFile = $DB . '/selfSimilarities.txt';

my $outputDir_templateDB = ($templateDB // '') . '/selfSimilarity';
my $outputDir_results_templateDB = $outputDir_templateDB . '/results';
my $outputFn_jobs_inTemplateDir = $outputDir_templateDB . '/jobs';

if(not $mode)
{
	die "Please specify --mode; e.g.: prepareFromScratch, prepareFromTemplate";
}
elsif($mode eq 'prepareFromScratch')
{
	(mkdir($outputDir) or die "Cannot mkdir $outputDir") unless(-d $outputDir);
	(mkdir($outputDir_computation) or die "Cannot mkdir $outputDir") unless(-d $outputDir_computation);
	(mkdir($outputDir_results) or die "Cannot mkdir $outputDir") unless(-d $outputDir_results);

	print "Prepare self-similarity computation from scratch - output directory $outputDir\n";
	
	# read taxonID -> contigs
	my %taxonID_2_contigs;
	my %contigLength;
	read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLength);
	
	# read taxonomy
	my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

	# strip down to mappable components
	taxTree::removeUnmappableParts($taxonomy, \%taxonID_2_contigs);
	
	# get leaves that we want to estimate potentially new node distances for
	my @nodesForAttachment = taxTree::getNodesForPotentialAttachmentOfNovelSpecies($taxonomy);
	
	# sanity check that really all leaf nodes are mappable
	my @nodesForAttachment_leaves = map {taxTree::descendants_leaves($taxonomy, $_)} @nodesForAttachment;
	die unless(all {exists $taxonID_2_contigs{$_}} @nodesForAttachment_leaves);
	print "Have nodes ", scalar(@nodesForAttachment), " that hypothetical new species could be attached to.\n";
	
	# open jobs file
	open(JOBS, '>', $outputFn_jobs) or die "Canot open $outputFn_jobs";
	
	# derive subcomputations for each node
	my $total_computations = 0;
	my $computations_max_per_node = 0;
	my $computations_max_per_node_which;
	foreach my $nodeID (sort @nodesForAttachment)
	{
		my @thisNode_self_simililarity_subcomputations = taxTree::getSubComputationsForAttachment($taxonomy, $nodeID, \%taxonID_2_contigs);
		die Dumper("Fewer subcomputations than expected", scalar(@thisNode_self_simililarity_subcomputations), scalar(@{$taxonomy->{$nodeID}{children}})) unless(scalar(@thisNode_self_simililarity_subcomputations) >= scalar(@{$taxonomy->{$nodeID}{children}}));
		
		my $computations_thisNode = scalar(@thisNode_self_simililarity_subcomputations);
		
		$total_computations += $computations_thisNode;
		if($computations_thisNode > $computations_max_per_node)
		{
			$computations_max_per_node = $computations_thisNode;
			$computations_max_per_node_which = $nodeID;
		}
		
		# print "Node $nodeID, $computations_thisNode computations.\n";
		
		foreach my $node_computation (@thisNode_self_simililarity_subcomputations)
		{
			my $targetLeafID = $node_computation->[1];
			my $targetLeaf_genomeSize = Util::getGenomeLength($targetLeafID, \%taxonID_2_contigs, \%contigLength);
			my @targetLeaf_compareAgainst = sort @{$node_computation->[2]};
			my @targetLeaf_contigs = sort keys %{$taxonID_2_contigs{$targetLeafID}};
			
			my $targetLeaf_compareAgainst_combinedSize = 0;
			my @compareAgainst_contigs;
			my %ranks_compareAgainst;
			foreach my $taxonID (@targetLeaf_compareAgainst)
			{
				$targetLeaf_compareAgainst_combinedSize += Util::getGenomeLength($taxonID, \%taxonID_2_contigs, \%contigLength);
				push(@compareAgainst_contigs, keys %{$taxonID_2_contigs{$taxonID}});
				$ranks_compareAgainst{$taxonomy->{$taxonID}{rank}}++;
			}
			@compareAgainst_contigs = sort @compareAgainst_contigs;
			
			print JOBS join("\t",
				$nodeID,							# the node we're attaching new species todo
				$taxonomy->{$nodeID}{rank},		
				$node_computation->[0],						# the node's immediate child we're taking the to-map-genome from
				$taxonomy->{$node_computation->[0]}{rank},
				$targetLeafID,									# the node which we're now mapping against the genomes of the other children
				$taxonomy->{$targetLeafID}{rank},
				$targetLeaf_genomeSize,						
				scalar(@targetLeaf_compareAgainst),
				join(";", @targetLeaf_compareAgainst),
				$targetLeaf_compareAgainst_combinedSize,
				join(";", @targetLeaf_contigs),
				join(";", @compareAgainst_contigs),
				join(';', sort keys %ranks_compareAgainst)
			), "\n";		
		}
	}
	
	print "Need to carry out $total_computations computations in total; max $computations_max_per_node per node ($computations_max_per_node_which, " . taxTree::taxon_id_get_name($computations_max_per_node_which, $taxonomy), ").\n";
	
	close(JOBS);
	
	print "Generated file $outputFn_jobs\n";
	
open(QSUB, '>', $outputFn_qsub) or die "Cannot open $outputFn_qsub";
print QSUB qq(#!/bin/bash
#\$ -t 1-${total_computations}
jobID=\$(expr \$SGE_TASK_ID - 1)
cd $path_to_cwd
perl $path_to_myself --mode doJobI --DB $DB --jobI \$jobID
);
close(QSUB);

	print "\nIf you use an SGE environment, you can do 'qsub $outputFn_qsub' to submit inidividual computation jobs.\n\n";
}
elsif($mode eq 'prepareFromTemplate')
{
	die "Please specify --templateDB" unless(defined $templateDB);
	die "Template DB $templateDB does not have self-similarity results" unless(-e $templateDB . '/selfSimilarities.txt');
	
	(mkdir($outputDir) or die "Cannot mkdir $outputDir") unless(-d $outputDir);
	(mkdir($outputDir_computation) or die "Cannot mkdir $outputDir") unless(-d $outputDir_computation);
	(mkdir($outputDir_results) or die "Cannot mkdir $outputDir") unless(-d $outputDir_results);
	
	# read existing computation details
	
	my $get_computation_key = sub {
		my $contigs_A_str = shift;
		my $contigs_B_str = shift;
		return $contigs_A_str . '///' . $contigs_B_str;
	};
	my $fn_template_jobs = $templateDB . '/selfSimilarity/jobs';
	die "Template DB $templateDB does not have self-similarity jobs file" unless(-e $fn_template_jobs);
	my %templateDB_existingComputation_toJobID;
	my %templateDB_attachableNodes;
	my %templateDB_A_to_B;
	{
		open(JOBS, '<', $fn_template_jobs) or die "Cannot open $outputFn_jobs";
		my $jobI = -1;
		while(<JOBS>)
		{
			my $line = $_;
			$jobI++;
			chomp($line);
			my @fields = split(/\t/, $line);
			die unless($#fields == 12);
			my $nodeID = $fields[0];
			my $A_taxonID = $fields[4];
			my $B_taxonIDs = $fields[8];			
			my $contigs_A = $fields[10];
			my $contigs_B = $fields[11];
			
			my $computation_key = $get_computation_key->($contigs_A, $contigs_B);
			die if(exists $templateDB_existingComputation_toJobID{$computation_key});
			$templateDB_existingComputation_toJobID{$computation_key} = $jobI;
			$templateDB_attachableNodes{$nodeID} = 1;
			
			my %contigs_B_href = map {$_ => 1} split(/;/, $contigs_B);
			push(@{$templateDB_A_to_B{$contigs_A}}, [\%contigs_B_href, $jobI]);
		}
		close(JOBS);	
	}
	
	print "Prepare self-similarity for $DB, based on $templateDB \n";
	
	# details for reduced DB
	my %reduced_taxonID_2_contigs;
	my %reduced_contigLength;
	read_taxonIDs_and_contigs($DB, \%reduced_taxonID_2_contigs, \%reduced_contigLength);
	
	# read reduced taxonomy
	my $reduced_taxonomy = taxTree::readTaxonomy($taxonomyDir);
	taxTree::removeUnmappableParts($reduced_taxonomy, \%reduced_taxonID_2_contigs);
	my @nodesForAttachment = taxTree::getNodesForPotentialAttachmentOfNovelSpecies($reduced_taxonomy);
	
	
	# template data
	my %template_taxonID_2_contigs;
	my %template_contigLength;
	my $template_taxonomy = taxTree::readTaxonomy($templateDB . '/taxonomy');	
	read_taxonIDs_and_contigs($templateDB, \%template_taxonID_2_contigs, \%template_contigLength);
	taxTree::removeUnmappableParts($template_taxonomy, \%template_taxonID_2_contigs);
	
	# check that this is a valid template
	foreach my $contigID_in_reduced (keys %reduced_contigLength)
	{
		die unless(defined $template_contigLength{$contigID_in_reduced});
	}
	
	# sanity check that really all leaf nodes are mappable
	my @nodesForAttachment_leaves = map {taxTree::descendants_leaves($reduced_taxonomy, $_)} @nodesForAttachment;
	die unless(all {exists $reduced_taxonID_2_contigs{$_}} @nodesForAttachment_leaves);
	
	print "\nHave nodes ", scalar(@nodesForAttachment), " that hypothetical new species could be attached to.\n";
	
	# derive subcomputations for each node
	# open jobs file
	open(JOBS, '>', $outputFn_jobs_fromTemplate) or die "Canot open $outputFn_jobs_fromTemplate";
		
	my @missingData_computations;
	my $total_computations = 0;
	my $total_computations_haveResult = 0;
	my $jobI = -1;	
	foreach my $nodeID (sort @nodesForAttachment)
	{
		die unless($templateDB_attachableNodes{$nodeID});
		my @thisNode_self_simililarity_subcomputations = taxTree::getSubComputationsForAttachment($reduced_taxonomy, $nodeID, \%reduced_taxonID_2_contigs);
		die Dumper("Fewer subcomputations than expected", scalar(@thisNode_self_simililarity_subcomputations), scalar(@{$reduced_taxonomy->{$nodeID}{children}})) unless(scalar(@thisNode_self_simililarity_subcomputations) >= scalar(@{$reduced_taxonomy->{$nodeID}{children}}));
		
		foreach my $node_computation (@thisNode_self_simililarity_subcomputations)
		{
			$jobI++;

			my $targetLeafID = $node_computation->[1];
			my $targetLeaf_genomeSize = Util::getGenomeLength($targetLeafID, \%reduced_taxonID_2_contigs, \%reduced_contigLength);
			my @targetLeaf_compareAgainst = sort @{$node_computation->[2]};
			my @targetLeaf_contigs = sort keys %{$reduced_taxonID_2_contigs{$targetLeafID}};
			
			my $targetLeaf_compareAgainst_combinedSize = 0;
			my @compareAgainst_contigs;
			my %ranks_compareAgainst;
			foreach my $taxonID (@targetLeaf_compareAgainst)
			{
				$targetLeaf_compareAgainst_combinedSize += Util::getGenomeLength($taxonID, \%reduced_taxonID_2_contigs, \%reduced_contigLength);
				push(@compareAgainst_contigs, keys %{$reduced_taxonID_2_contigs{$taxonID}});
				$ranks_compareAgainst{$reduced_taxonomy->{$taxonID}{rank}}++;
			}
			
			@compareAgainst_contigs = sort @compareAgainst_contigs;
			
			my $contigs_A = join(';', @targetLeaf_contigs);
			my $computation_key = $get_computation_key->($contigs_A, join(';', @compareAgainst_contigs));
			$total_computations++;
			
			print JOBS join("\t",
				$nodeID,							# the node we're attaching new species todo
				$reduced_taxonomy->{$nodeID}{rank},		
				$node_computation->[0],						# the node's immediate child we're taking the to-map-genome from
				$reduced_taxonomy->{$node_computation->[0]}{rank},
				$targetLeafID,									# the node which we're now mapping against the genomes of the other children
				$reduced_taxonomy->{$targetLeafID}{rank},
				$targetLeaf_genomeSize,						
				scalar(@targetLeaf_compareAgainst),
				join(";", @targetLeaf_compareAgainst),
				$targetLeaf_compareAgainst_combinedSize,
				join(";", @targetLeaf_contigs),
				join(";", @compareAgainst_contigs),
				join(';', sort keys %ranks_compareAgainst)
			), "\n";			
			
			if(exists $templateDB_existingComputation_toJobID{$computation_key})
			{
				$total_computations_haveResult++;
				my $template_jobI = $templateDB_existingComputation_toJobID{$computation_key};
				my $existing_results_file = get_results_file_for_jobI_fromTemplateDB($template_jobI);
				my $new_results_file = get_results_file_for_jobI($jobI);
				# print "Have 1:1 template for $computation_key - copy $existing_results_file --> $new_results_file \n";
				copy($existing_results_file, $new_results_file);
			}
			else
			{
				die unless(exists $templateDB_A_to_B{$contigs_A});
				my %alternativeI_to_distance;
				for(my $alternativeI = 0; $alternativeI <= $#{$templateDB_A_to_B{$contigs_A}}; $alternativeI++)
				{
					my $alternative = $templateDB_A_to_B{$contigs_A}[$alternativeI][0];
					next unless(all {exists $alternative->{$_}} @compareAgainst_contigs);
					my $d = scalar(keys %$alternative) - scalar(@compareAgainst_contigs);
					$alternativeI_to_distance{$alternativeI} = $d; 
				}
				
				die unless(keys %alternativeI_to_distance);
				my @sorted_alternative_Is = sort {$alternativeI_to_distance{$a} <=> $alternativeI_to_distance{$a}} keys %alternativeI_to_distance;
				if(scalar(@sorted_alternative_Is) > 1)
				{
					die unless($alternativeI_to_distance{$sorted_alternative_Is[0]} <= $alternativeI_to_distance{$sorted_alternative_Is[1]});
				}
				my $best_template_idx = $sorted_alternative_Is[0];
				my $best_template_jobI = $templateDB_A_to_B{$contigs_A}[$best_template_idx][1];
				
				# print "Don't have 1:1 template for $computation_key \n";
				# print "\tClosest existing computation: $best_template_jobI with distance $alternativeI_to_distance{$sorted_alternative_Is[0]}\n";
				# print "\t\t", join(' ', @compareAgainst_contigs), "\n";
				# print "\t\t", join(' ', keys %{$templateDB_A_to_B{$contigs_A}->[$best_template_idx][0]}), "\n\n";
				
				# print "\t--mode doJobIFromTemplate --jobI $jobI --jobITemplate $best_template_jobI --DB $DB --templateDB $templateDB \n";
				
				push(@missingData_computations, [$jobI, $best_template_jobI]);
			}	
		}
	}
	close(JOBS);
	
	print "\nTotal required computations: $total_computations / of which we have results for from template: $total_computations_haveResult\n";
	
	print "\n";
	foreach my $computation (@missingData_computations)
	{
		my $jobI = $computation->[0];
		my $best_template_jobI = $computation->[1];
		# print "Do \t--mode doJobIFromTemplate --jobI $jobI --jobITemplate $best_template_jobI --DB $DB --templateDB $templateDB \n";	
		doJobIFromTemplate($jobI, $best_template_jobI);
	}
	
	doCollect(1);
}
elsif($mode eq 'doJobIFromTemplate')
{
	doJobIFromTemplate($jobI, $jobITemplate); 
}
elsif($mode eq 'doJobI')
{
	die "Please specify --jobI" unless(defined $jobI);
	
	my ($A_taxonID, $B_taxonIDs, $contigs_A, $contigs_B) = get_job_info($outputFn_jobs, $jobI);

	print "Job $jobI";
	print "\tMapping from: $A_taxonID\n";
	print "\tMapping to  : $B_taxonIDs\n";
	
	my @contigIDs_A = split(/;/, $contigs_A);
	my @contigIDs_B = split(/;/, $contigs_B);

	
	my %contigIDs_A_to_i = Util::get_index_hash(\@contigIDs_A);
	my %contigIDs_B_to_i = Util::get_index_hash(\@contigIDs_B);
	
	my $dir_computation = $outputDir_computation . '/' . $jobI;
	mkdir($dir_computation);

	my $file_A = $dir_computation . '/A';
	my $file_B = $dir_computation . '/B';
	
	construct_A_B_files($file_ref, $file_A, $file_B, \@contigIDs_A, \@contigIDs_B);
	
	my $results_fn_reads_many = get_results_file_for_jobI($jobI);
	my $readInfo_fn_reads_many = get_readInfo_file_for_jobI($jobI);
	my $readResults_fn_reads_many = get_readResults_file_for_jobI($jobI);

	my $seed_rand = rand(1024); # todo
	srand($seed_rand);
				
	open(READINFO, '>', $readInfo_fn_reads_many) or die "Cannot open $readInfo_fn_reads_many";
	print READINFO $contigs_A, "\t", $contigs_B, "\n";	
	print READINFO join("\t", $readSimSizeFrom, $readSimSizeTo, $readSimSizeStep), "\n";
	print READINFO $readSimDelta, "\n";
	print READINFO $seed_rand, "\n";
	
	open(READDETAILRESULTS, '>', $readResults_fn_reads_many) or die "Cannot open $readResults_fn_reads_many";

	my $sub_call_each_simulatedRead = sub { 
		my ($readID, $contigID, $posI, $chunkLength) = @_;
		my $contigI = $contigIDs_A_to_i{$contigID};
		die unless(defined $contigI);
		print READINFO $contigI, "\t", $posI, "\n";
	};
				
	my $sub_call_each_readResult = sub {
		my ($readID, $bestContig, $bestIdentity) = @_;
		my $bestContig_index = $contigIDs_B_to_i{$bestContig};
		die unless(defined $bestContig_index);
		print READDETAILRESULTS join("\t", $readID, $bestContig_index, $bestIdentity), "\n";	
	};
	
	my $identities_by_length_href = {};
	map_reads_keepTrack_similarity(
		$file_A,
		$file_B,
		$dir_computation,
		$seed_rand,
		\@contigIDs_A,
		undef,
		$sub_call_each_simulatedRead,
		$sub_call_each_readResult,
		$identities_by_length_href
	);
	
	open(RESULTS, '>', $results_fn_reads_many) or die "Cannot open $results_fn_reads_many";
	print RESULTS $contigs_A, "\t", $contigs_B, "\n";	
	foreach my $chunkLength (sort {$a <=> $b} keys %$identities_by_length_href)
	{
		foreach my $k (sort {$a <=> $b} keys %{$identities_by_length_href->{$chunkLength}})
		{
			print RESULTS join("\t", $chunkLength, $k, $identities_by_length_href->{$chunkLength}{$k}), "\n";
		}	
	}
	close(RESULTS);
	
	print "--------\nJob $jobI done. Produced $results_fn_reads_many.\n";
	close(READDETAILRESULTS);		
	close(READINFO);
}
elsif($mode eq 'collect')
{
	doCollect();
}
else
{
	die "Unknown mode";
}

sub map_reads_keepTrack_similarity
{
	my $file_A = shift;
	my $file_B = shift;
	my $tmpDir = shift;
	my $v_for_srand = shift;
	my $contigs_A_order = shift;
	my $limitToReadsIDs = shift;
	my $call_eachSimulatedRead = shift;
	my $call_eachReadResult = shift;
	my $readCounts_toPopulate_href = shift;
	
	die unless(defined $v_for_srand);
	die unless(defined $readCounts_toPopulate_href);
	
	my $A_contigs_href = Util::readFASTA($file_A);
	
	my $file_A_reads_many = $tmpDir . '/A.reads.many';
	my $outputFn_MetaMap = $tmpDir . '/MetaMap.many';	
	
	open(READS, '>', $file_A_reads_many) or die "Cannot open $file_A_reads_many";
	
	my @chunk_positions = getChunkPositions($file_A, $contigs_A_order, $v_for_srand);
	my %readID_2_chunkLength;
	my %generated_reads_chunkLength;
	foreach my $oneChunk (@chunk_positions)
	{
		my $chunkLength = $oneChunk->[0];
		my $readI_within_chunkLength_1based = $oneChunk->[1];
		my $contigID = $oneChunk->[2];
		my $posI = $oneChunk->[3];
		my $readID = $oneChunk->[4];
		
		if(defined $limitToReadsIDs)
		{
			next unless($limitToReadsIDs->{$readID});
		}			
		
		if(defined $call_eachSimulatedRead)
		{
			$call_eachSimulatedRead->($readID, $contigID, $posI, $chunkLength);
		}		

		die unless(defined $A_contigs_href->{$contigID});			
		my $S = substr($A_contigs_href->{$contigID}, $posI, $chunkLength);
		die unless(length($S) == $chunkLength);
		print READS '>', $readID, "\n";
		print READS $S, "\n";
		
		$generated_reads_chunkLength{$chunkLength}++;	
		$readID_2_chunkLength{$readID} = $chunkLength;
	}
	close(READS);
		
	print "\t\tReads file $file_A_reads_many with ", sum(values %generated_reads_chunkLength), " reads";

	# important - without this we get a different sort order!

	my $MetaMap_cmd = qq($metamap_bin mapDirectly -r $file_B -q $file_A_reads_many -m 2000 -o $outputFn_MetaMap --pi 80);
	print "\t\tExecuted command $MetaMap_cmd\n";		
	system($MetaMap_cmd) and die "MetaMap command $MetaMap_cmd failed";
	
	my %readCounts_toPopulate_local;
	my $currentReadID = '';
	my @currentReadLines;
	my $processAlignments_oneRead = sub {
		my $bestIdentity;
		my $bestIdentity_which;
		my $readID;
		foreach my $line (@currentReadLines)
		{
			my @fields = split(/ /, $line);
			die Dumper($fields[0], $currentReadID) unless($fields[0] eq $currentReadID);
			
			if(not defined $readID)
			{
				$readID = $fields[0];
			}
			else
			{
				die unless($readID eq $fields[0]);
			}
			my $targetContig = $fields[5];
			my $identity = $fields[9];
			die unless(($identity >= 0) and ($identity <= 100));
			if((not defined $bestIdentity) or ($bestIdentity < $identity))
			{
				$bestIdentity = $identity;
				$bestIdentity_which = $targetContig;
			}
		}
		$bestIdentity = int($bestIdentity + 0.5);
		
		die unless(defined $readID_2_chunkLength{$readID});
		
		my $chunkLength = $readID_2_chunkLength{$readID};
		$readCounts_toPopulate_local{$chunkLength}{$bestIdentity}++;
			
		if($call_eachReadResult)
		{
			$call_eachReadResult->($readID, $bestIdentity_which, $bestIdentity);
		}
	};	
		
	open(MetaMapOUTPUT, '<', $outputFn_MetaMap) or die "Cannot open $outputFn_MetaMap";
	while(<MetaMapOUTPUT>)
	{
		chomp;
		next unless($_);
		die "Weird input" unless($_ =~ /^(.+?) /);
		my $readID = $1;
		if($currentReadID ne $readID)
		{
			if(@currentReadLines)
			{
				$processAlignments_oneRead->();
			}
			$currentReadID = $readID;
			@currentReadLines = ();
		}
		push(@currentReadLines, $_);
	}
	if(@currentReadLines)
	{
		$processAlignments_oneRead->();
	}	
	
	close(MetaMapOUTPUT);

	foreach my $chunkLength (keys %generated_reads_chunkLength)
	{
		my $reads_for_mapping = $generated_reads_chunkLength{$chunkLength};
		my $n_read_alignments = 0;
		if(exists $readCounts_toPopulate_local{$chunkLength})
		{
			$n_read_alignments = sum values %{$readCounts_toPopulate_local{$chunkLength}};
		}	
		my $n_missing_reads = $reads_for_mapping - $n_read_alignments;
		die Dumper("Problem missing reads", $reads_for_mapping, $n_read_alignments) unless($n_missing_reads >= 0);
		$readCounts_toPopulate_local{$chunkLength}{0} += $n_missing_reads;				
	}
	
	foreach my $chunkLength (keys %readCounts_toPopulate_local)
	{
		foreach my $idty (keys %{$readCounts_toPopulate_local{$chunkLength}})
		{
			my $count = $readCounts_toPopulate_local{$chunkLength}{$idty};
			$readCounts_toPopulate_href->{$chunkLength}{$idty} += $count;
		}
	}
	
	system('rm ' . $tmpDir . '/*') and die "Cannot delete $tmpDir (I)";	
	system('rm -r ' . $tmpDir) and die "Cannot delete $tmpDir (II)"
}

sub getChunkPositions
{
	my $file_A = shift;
	my $contigs_A_order = shift;	
	my $srand_value = shift;
	
	my $A_contigs_href = Util::readFASTA($file_A);
		
	# todo
	srand(length(join(';', @$contigs_A_order)));
	
	# srand($srand_value);
	
	my @forReturn;
	
	my $totalReadI = 0;
	for(my $chunkLength = $readSimSizeFrom; $chunkLength <= $readSimSizeTo; $chunkLength += $readSimSizeStep)
	{						
		my $read_start_positions  = 0;
		die unless(scalar(keys %$A_contigs_href) == scalar(@$contigs_A_order));
		foreach my $contigID (@$contigs_A_order)
		{
			my $contigSequence = $A_contigs_href->{$contigID};
			die unless(defined $contigSequence);
			for(my $posI = 0; $posI < length($contigSequence); $posI += $readSimDelta)
			{
				my $lastPos = $posI + $chunkLength - 1;
				if($lastPos < length($contigSequence))
				{
					$read_start_positions++;
				}
			}  
		}	
		
		my $start_rate = 1;
		if($read_start_positions > $target_max_simulatedChunks)
		{
			$start_rate = $target_max_simulatedChunks / $read_start_positions;
			# print "\t\t\t(Read lengths; $chunkLength) Adjusted start rate to $start_rate (eligible start positions: $read_start_positions, want $target_max_simulatedChunks)\n";
		}
		die unless(($start_rate >= 0) and ($start_rate <= 1));
		
		my $readI_with_chunkLength = 0;
		for(my $contigI = 0; $contigI <= $#{$contigs_A_order}; $contigI++)
		{
			my $contigID = $contigs_A_order->[$contigI];
			my $contigSequence = $A_contigs_href->{$contigID};
			POS: for(my $posI = 0; $posI < length($contigSequence); $posI += $readSimDelta)
			{
				my $lastPos = $posI + $chunkLength - 1;
				if($lastPos < length($contigSequence))
				{
					if($start_rate != 1)
					{
						next POS if(rand(1) > $start_rate);
					}
					
					$readI_with_chunkLength++;
					$totalReadI++;
					my $readID = 'read' . $totalReadI;

					push(@forReturn, [$chunkLength, $readI_with_chunkLength, $contigID, $posI, $readID]);
				}
			}  
		}
	}
	
	return @forReturn;
}

sub doJobIFromTemplate
{
	my $jobI = shift;
	my $jobITemplate = shift;
	
	print "Carying out job $jobI (DB $DB) based on template $jobITemplate ($templateDB)\n";
	
	die "Please specify --jobI" unless(defined $jobI);
	die "Please specify --jobITemplate" unless(defined $jobITemplate);
	
	my ($A_taxonID, $B_taxonIDs, $contigs_A, $contigs_B) = get_job_info($outputFn_jobs_fromTemplate, $jobI);
	my ($template_A_taxonID, $template_B_taxonIDs, $template_contigs_A, $template_contigs_B) = get_job_info($outputFn_jobs_inTemplateDir, $jobITemplate);

	my @template_contigIDs_A = split(/;/, $template_contigs_A);
	my @template_contigIDs_B = split(/;/, $template_contigs_B);
	my %template_contigIDs_A_to_i = Util::get_index_hash(\@template_contigIDs_A);
	my %template_contigIDs_B_to_i = Util::get_index_hash(\@template_contigIDs_B);
	my %template_A_i_to_contigID = reverse %template_contigIDs_A_to_i;
	my %template_B_i_to_contigID = reverse %template_contigIDs_B_to_i;
	
	unless($contigs_A eq $template_contigs_A)
	{
		die "File $jobITemplate is not a valid template for job $jobI -- expect contigs_A = $contigs_A, but got $template_contigs_A";
	}
	
	my @contigIDs_A = split(/;/, $contigs_A);
	my @contigIDs_B = split(/;/, $contigs_B);
	my %contigs_A = map {$_ => 1} @contigIDs_A; die unless(scalar(keys %contigs_A));
	my %contigs_B = map {$_ => 1} @contigIDs_B; die unless(scalar(keys %contigs_B));
	my %contigIDs_A_to_i = Util::get_index_hash(\@contigIDs_A);
	my %contigIDs_B_to_i = Util::get_index_hash(\@contigIDs_B);
	
	foreach my $contigID (@contigIDs_B)
	{
		die Dumper("Invalid template $jobITemplate for job $jobI - mapping target (B) $contigID is not in list of template mapping targets", \@template_contigIDs_B)  unless(exists $template_contigIDs_B_to_i{$contigID});
	}
	
	my ($_tmp_A, $_tmp_B, $v_for_srand, $readPositions_aref) = check_readsInfo_compatible_and_get_contigsA_contigsB_srand_readPositions(get_readInfo_file_for_jobI_fromTemplateDB($jobITemplate));
	die Dumper("A contig discrepancy", $_tmp_A, $contigs_A) unless($_tmp_A eq $contigs_A);
	die Dumper("B contig discrepancy", $_tmp_B, $template_contigs_B) unless($_tmp_B eq $template_contigs_B);
	

	# get reads, and read info
	my $dir_computation = $outputDir_computation . '/' . $jobI;
	mkdir($dir_computation);
	
	my $file_A = $dir_computation . '/A';
	my $file_B = $dir_computation . '/B';
		
	construct_A_B_files($file_ref, $file_A, $file_B, \@contigIDs_A, \@contigIDs_B);
	
	my %n_reads_per_chunkLength;
	my %readID_2_info;
	my @chunk_positions = getChunkPositions($file_A, \@contigIDs_A, $v_for_srand);
	
	# todo remove
	# my @chunk_positions_2 = getChunkPositions($file_A, \@contigIDs_A, $v_for_srand);
	
	my $readI = -1;
	foreach my $oneChunk (@chunk_positions)
	{
		my $chunkLength = $oneChunk->[0];
		my $readI_within_chunkLength_1based = $oneChunk->[1];
		my $contigID = $oneChunk->[2];
		my $posI = $oneChunk->[3];
		my $readID = $oneChunk->[4];	
		
		$n_reads_per_chunkLength{$chunkLength}++;
		die if(defined $readID_2_info{$readID});
		$readID_2_info{$readID} = $oneChunk;
		
		$readI++;
			
		# todo remove
		# my $contigID_supposed = $template_A_i_to_contigID{$readPositions_aref->[$readI][0]};
		# print Dumper([$contigID, [$posI, $chunk_positions_2[$readI][3]]], [$readPositions_aref->[$readI][0], $contigID_supposed, $readPositions_aref->[$readI][1]], $contigs_A, $template_contigs_A), "\n";
		
		if(defined $readPositions_aref)
		{
			die unless($posI == $readPositions_aref->[$readI][1]);
		}
	}
	
	
	# find reads that need to be remapped, and fill histogram with non-remapped reads
	my %reads_for_remapping;
	my $identities_by_length_href = {};
	{
		my $fn_mappings_details = get_readResults_file_for_jobI_fromTemplateDB($jobITemplate);
		open(F, '<', $fn_mappings_details) or die "Cannot open $fn_mappings_details";
		while(<F>)
		{
			chomp;
			next unless($_);
			my @f = split(/\t/, $_);
			die Dumper("Weird line in $fn_mappings_details", $_) unless(scalar(@f) == 3);
			my $readID = $f[0];
			my $contigI = $f[1];
			my $bestIdentity = $f[2];
			
			my $contigID = $template_B_i_to_contigID{$contigI};
			die Dumper("Unknown contigI $contigI", $fn_mappings_details, $., $_) unless(defined $contigID);
			
			die unless(defined $readID_2_info{$readID});
			if(exists $contigIDs_B_to_i{$contigID})
			{
				my $chunkLength = $readID_2_info{$readID}[0];
				$identities_by_length_href->{$chunkLength}{$bestIdentity}++;
			}
			else
			{
				$reads_for_remapping{$readID} = 1;
			}
			
		}
		close(F);
	}
		
	my $sub_call_each_simulatedRead = sub { 
		my ($readID, $contigID, $posI, $chunkLength) = @_;
		if($readPositions_aref)
		{
			die unless($readID =~ /read(\d+)/);
			my $readI = $1 - 1;
			die "Weird readI $readI $readID" unless($readI >= 0);
			die "Weird readI $readI $readID " . scalar(@$readPositions_aref) unless($readI < scalar(@$readPositions_aref));
			my $contigI = $contigIDs_A_to_i{$contigID};
			my $contigII = $template_contigIDs_A_to_i{$contigID};
			
			die unless($contigI == $contigII);
			die unless($readPositions_aref->[$readI][0] == $contigI);
			die Dumper("Mismatch in generated reads", [$readID, $contigID, $contigI, $posI, $chunkLength], $readPositions_aref->[$readI], $v_for_srand) unless($readPositions_aref->[$readI][1] == $posI);
		}
	};

				
	my $sub_call_each_readResult = sub {	
	};
	
	print "\tNow remapping ", scalar(keys %reads_for_remapping), " / ", scalar(@chunk_positions), " genome chunks.\n";
	
	map_reads_keepTrack_similarity(
		$file_A,
		$file_B,
		$dir_computation,
		$v_for_srand,
		\@contigIDs_A,
		\%reads_for_remapping,
		$sub_call_each_simulatedRead,
		$sub_call_each_readResult,
		$identities_by_length_href
	);
	
	foreach my $chunkLength (keys %n_reads_per_chunkLength)
	{
		my $n_trackedReads = 0;
		if(defined $identities_by_length_href->{$chunkLength})
		{
			$n_trackedReads = sum values %{$identities_by_length_href->{$chunkLength}};
		}
		my $n_missingReads = $n_reads_per_chunkLength{$chunkLength} - $n_trackedReads;
		die unless($n_missingReads >= 0);
		$identities_by_length_href->{$chunkLength}{0} += $n_missingReads;
	}
	
	my $results_fn_reads_many = get_results_file_for_jobI($jobI);
	
	open(RESULTS, '>', $results_fn_reads_many) or die "Cannot open $results_fn_reads_many";
	print RESULTS $contigs_A, "\t", $contigs_B, "\n";		
	foreach my $chunkLength (sort {$a <=> $b} keys %$identities_by_length_href)
	{
		foreach my $k (sort {$a <=> $b} keys %{$identities_by_length_href->{$chunkLength}})
		{
			print RESULTS join("\t", $chunkLength, $k, $identities_by_length_href->{$chunkLength}{$k}), "\n";
		}	
	}
	close(RESULTS);	
	
	print "\tDone. Produced file ", $results_fn_reads_many, "\n";
}

sub get_results_file_for_jobI
{
	my $jobI = shift;
	return  $outputDir_results . '/' . $jobI . '.results.reads.many';;
}

sub get_results_file_for_jobI_fromTemplateDB
{
	my $jobI = shift;
	return  $outputDir_results_templateDB . '/' . $jobI . '.results.reads.many';;
}

sub get_readInfo_file_for_jobI
{
	my $jobI = shift;
	return  $outputDir_results . '/' . $jobI . '.results.reads.many.readInfo';
}

sub get_readInfo_file_for_jobI_fromTemplateDB
{
	my $jobI = shift;
	return  $outputDir_results_templateDB . '/' . $jobI . '.results.reads.many.readInfo';
}


sub get_readResults_file_for_jobI
{
	my $jobI = shift;
	return  $outputDir_results . '/' . $jobI . '.results.reads.many.readResults';
}

sub get_readResults_file_for_jobI_fromTemplateDB
{
	my $jobI = shift;
	return  $outputDir_results_templateDB . '/' . $jobI . '.results.reads.many.readResults';
}


sub read_taxonIDs_and_contigs
{
	my $DB = shift;
	my $taxonID_2_contigs_href = shift;
	my $contigLength_href = shift;
	
	my $file_taxonGenomes = $DB . '/taxonInfo.txt';
	
	open(GENOMEINFO, '<', $file_taxonGenomes) or die "Cannot open $file_taxonGenomes";
	while(<GENOMEINFO>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/ /, $line);
		die unless(scalar(@line_fields) == 2);
		my $taxonID = $line_fields[0];
		my $contigs = $line_fields[1];
		die if(exists $taxonID_2_contigs_href->{$taxonID});

		my @components = split(/;/, $contigs);
		foreach my $component (@components)
		{
			my @p = split(/=/, $component);
			die unless(scalar(@p) == 2);
			die if(exists $taxonID_2_contigs_href->{$taxonID}{$p[0]});
			$taxonID_2_contigs_href->{$taxonID}{$p[0]} = $p[1];
			$contigLength_href->{$p[0]} = $p[1];
		}
	}
	close(GENOMEINFO);
}


sub get_job_info
{
	my $fn = shift;
	my $jobI = shift;
	
	my ($A_taxonID, $B_taxonIDs, $contigs_A, $contigs_B);
	open(JOBS, '<', $fn) or die "Cannot open $fn";
	my $lineI = -1;
	while(<JOBS>)
	{
		my $line = $_;
		$lineI++;
		if($lineI == $jobI)
		{
			chomp($line);
			my @fields = split(/\t/, $line);
			die unless($#fields == 12);
			$A_taxonID = $fields[4];
			$B_taxonIDs = $fields[8];			
			$contigs_A = $fields[10];
			$contigs_B = $fields[11];
		}
	}
	close(JOBS);
	die "Can't find job data for $jobI" unless(defined $contigs_A);
	
	return ($A_taxonID, $B_taxonIDs, $contigs_A, $contigs_B);
}

sub check_readsInfo_compatible_and_get_contigsA_contigsB_srand_readPositions
{
	my $fn = shift;
	
	my $aref_read_positions = undef;
	
	open(F, '<', $fn) or die "Cannot open $fn";
	my $firstLine = <F>; chomp($firstLine);
	my $secondLine = <F>; chomp($secondLine);
	my $thirdLine = <F>; chomp($thirdLine);
	my $fourthLine = <F>; chomp($fourthLine);
	
	my @firstLine_fields = split(/\t/, $firstLine); die unless(scalar(@firstLine_fields) == 2);
	my @secondLine_fields = split(/\t/, $secondLine);
	die unless($secondLine_fields[0] == $readSimSizeFrom);
	die unless($secondLine_fields[1] == $readSimSizeTo);
	die unless($secondLine_fields[2] == $readSimSizeStep);
	
	die unless($thirdLine == $readSimDelta);
	
	while(<F>)
	{
		chomp;
		next unless($_);
		my @f = split(/\t/, $_);
		die unless(scalar(@f) == 2);
		$aref_read_positions = [] unless(defined $aref_read_positions);
		push(@$aref_read_positions, \@f);
	}
	close(F);
	
	return ($firstLine_fields[0], $firstLine_fields[1], $fourthLine, $aref_read_positions);
}


sub construct_A_B_files
{
	my $ref = shift;
	my $file_A = shift;
	my $file_B = shift;
	my $contigs_A_aref = shift;
	my $contigs_B_aref = shift;

	my %contigs_A = map {$_ => 1} @$contigs_A_aref; die unless(scalar(keys %contigs_A));
	my %contigs_B = map {$_ => 1} @$contigs_B_aref; die unless(scalar(keys %contigs_B));	
	
	open(A, '>', $file_A) or die;
	open(B, '>', $file_B) or die;
	open(REF, '<', $file_ref) or die "Cannot open $file_ref";
	
	my $in_A = 0;
	my $in_B = 0;
	while(<REF>)
	{
		my $line = $_;
		chomp($line);
		if(substr($line, 0, 1) eq '>')
		{
			$in_A = 0;	
			$in_B = 0;
			my $contigID = substr($line, 1);
			if($contigs_A{$contigID})
			{
				$in_A = 1;
				$contigs_A{$contigID}--;
			}
			if($contigs_B{$contigID})
			{
				$in_B = 1;
				$contigs_B{$contigID}--;
			}
			die if($in_A and $in_B);
		}

		if($in_A)
		{
			print A $line, "\n";
		}

		if($in_B)
		{
			print B $line, "\n";
		}
	}
	close(REF);
	close(A);
	close(B);

	foreach my $contigID (keys %contigs_A)
	{
		die "Missed contig $contigID" if($contigs_A{$contigID});
	}

	foreach my $contigID (keys %contigs_B)
	{
		die "Missed contig $contigID" if($contigs_B{$contigID});
	}	
}

sub doCollect
{
	my $collectFromTemplate = shift;
	
	print "Collect results from folder $outputDir_results\n";
	
	my %taxonID_2_contigs;
	my %contigLength;
	read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLength);
	
	# read taxonomy
	my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

	# strip down to mappable components
	taxTree::removeUnmappableParts($taxonomy, \%taxonID_2_contigs);
	
	my %results_reads_many_per_node;
	my $total_jobs = 0;
	my $total_jobs_reads_many_ok = 0;	
	my $jobs_fn = ($collectFromTemplate) ? $outputFn_jobs_fromTemplate : $outputFn_jobs;
	open(JOBS, '<', $jobs_fn) or die "Canot open $jobs_fn";
	my $lineI = -1;
	while(<JOBS>)
	{
		my $line = $_;
		chomp($line);
		my @fields = split(/\t/, $line);
		my $nodeID = $fields[0];
		my $contigs_A = $fields[10];
		my $contigs_B = $fields[11];		
		
		$lineI++;
		my $jobI = $lineI; 
		
		my $results_reads_many_fn = get_results_file_for_jobI($jobI);
		$total_jobs++; 
		
		if(-e $results_reads_many_fn)
		{
			my $S = 0;
			my %S_by_readLength;
			my %h;
			open(R, '<', $results_reads_many_fn) or die;
			my $firstLine = <R>;
			chomp($firstLine);
			my @firstLine_fields = split(/\t/, $firstLine);
			die Dumper("Weird first line", $results_reads_many_fn, \@firstLine_fields) unless(scalar(@firstLine_fields) == 2);
			die unless($firstLine_fields[0] eq $contigs_A);
			die unless($firstLine_fields[1] eq $contigs_B);
			while(<R>)
			{
				my $l = $_;
				chomp($l);
				next unless($l);
				my @fields = split(/\t/, $l);
				die Dumper("Weird field count in line $. of file $results_reads_many_fn", \@fields) unless(scalar(@fields) == 3);
				$h{$fields[0]}{$fields[1]} = $fields[2];
				$S += $fields[2];
				$S_by_readLength{$fields[0]} += $fields[2];
			}
			close(R);
			if($S == 0)
			{
				# warn "Problem with $results_reads_fn";
			}
			else
			{
				foreach my $readLength (keys %h)
				{
					next if(not $S_by_readLength{$readLength});
					foreach my $k (keys %{$h{$readLength}})
					{
						$h{$readLength}{$k} /= $S_by_readLength{$readLength};
					}
					push(@{$results_reads_many_per_node{$readLength}{$nodeID}}, $h{$readLength});
				}
			}
			
			$total_jobs_reads_many_ok++;			
		}
						
	}
	close(JOBS);
	
	print "\nTotal jobs: $total_jobs\n";
	print "\tof which with results: $total_jobs_reads_many_ok \n";

	open(RESULTSREADSMANY, '>', $outputFn_reads_results_many) or die "Cannot open $outputFn_reads_results_many";

	foreach my $readLength (keys %results_reads_many_per_node)
	{
		foreach my $nodeID (keys %{$results_reads_many_per_node{$readLength}})
		{
			my $nodeRank = $taxonomy->{$nodeID}{rank};
			my $nodeName = join('; ', @{$taxonomy->{$nodeID}{names}});

			my @descendants = taxTree::descendants($taxonomy, $nodeID);
			my @descendants_with_genomes = grep {exists $taxonID_2_contigs{$_}} @descendants;
			my %combinedHistogram;
			my $S = 0;		 
			foreach my $histogram (@{$results_reads_many_per_node{$readLength}{$nodeID}})
			{
				foreach my $k (keys %$histogram)
				{	
					$combinedHistogram{$k} += $histogram->{$k};
					$S += $histogram->{$k};
				}
			}
			
			my $firstKey = 1;				
			foreach my $k (keys %combinedHistogram)
			{
				$combinedHistogram{$k} /= $S;
				my $string_source_genomes = ($firstKey) ? join(';', @descendants_with_genomes) : '';
				my $string_rank = ($firstKey) ? $nodeRank : '';
				my $string_name = ($firstKey) ? $nodeName : '';
				$firstKey = 0;
				print RESULTSREADSMANY join("\t", $nodeID, $readLength, $k, $combinedHistogram{$k}, $string_source_genomes, $string_rank, $string_name), "\n";
			}
		}	
	}
	
	close(RESULTSREADSMANY);
	

	copy($outputFn_reads_results_many, $finalResultsFile) or die "Cannot copy $outputFn_reads_results_many -> $finalResultsFile";

	print "\n\nProduced results file:\n";
	print " - $finalResultsFile \n\n";
}
__END__


	my @nodes_consider_hypothetical_new_child = grep {scalar(@{$taxonomy->{$_}{children}}) > 1} keys %nodes_upTo_family;

	print "Consider ", scalar(@nodes_consider_hypothetical_new_child), " nodes for hypothetical-new-child.\n";

	open(JOBS, '>', $outputFn_jobs) or die "Canot open $outputFn_jobs";
	my $total_comparisons = 0;
	my $total_storage = 0;
	my $immediateDescendants_rank_heterogeneity = 0;
	for(my $nodeI = 0; $nodeI <= $#nodes_consider_hypothetical_new_child; $nodeI++)
	{
		print "\r\tProcessing $nodeI   ";
		my $nodeID = $nodes_consider_hypothetical_new_child[$nodeI];

		my @descendants = taxTree::descendants($taxonomy, $nodeID);
		my @descendants_with_genomes = grep {exists $taxonID_2_contigs{$_}} @descendants;
		die unless(scalar(@descendants_with_genomes) >= 2);
		
		my @immediate_descendants = @{$taxonomy->{$nodeID}{children}};
		my @immediate_descendants_ranks = map {$taxonomy->{$_}{rank}} @immediate_descendants;
		my %_descendants_ranks = map {$_ => 1} @immediate_descendants_ranks;

		if(scalar(keys %_descendants_ranks) > 1)
		{
			$immediateDescendants_rank_heterogeneity++;
			# warn Dumper(\%_descendants_ranks, \@immediate_descendants, \@immediate_descendants_ranks);
		}

		die unless(scalar(@immediate_descendants >= 2));

		my %leave_to_immediate_descendant_i;
		for(my $descendantI = 0; $descendantI <= $#immediate_descendants; $descendantI++)
		{
			my $descendantID = $immediate_descendants[$descendantI];
			
			if(exists $taxonID_2_contigs{$descendantID})
			{
				$leave_to_immediate_descendant_i{$descendantID} = $descendantI;
			}
			
			if(scalar(@{$taxonomy->{$descendantID}{children}}) == 0)
			{
				# $leave_to_immediate_descendant_i{$descendantID} = $descendantI;
			}
			else
			{
				my @immedidateDescendant_descendants = taxTree::descendants($taxonomy, $descendantID);
				my @immedidateDescendant_descendants_with_genomes = grep {exists $taxonID_2_contigs{$_}} @immedidateDescendant_descendants;
				
				foreach my $leafID (@immedidateDescendant_descendants_with_genomes)
				{
					die if(exists $leave_to_immediate_descendant_i{$leafID});
					$leave_to_immediate_descendant_i{$leafID} = $descendantI;
				}
			}
		}
		die Dumper("Weird", scalar(@descendants_with_genomes), scalar(keys %leave_to_immediate_descendant_i), \@descendants_with_genomes, \%leave_to_immediate_descendant_i, $taxonomy->{$nodeID}) unless(scalar(keys %leave_to_immediate_descendant_i) == scalar(@descendants_with_genomes));

		for(my $leafI = 0; $leafI <= $#descendants_with_genomes; $leafI++)
		{
			my $leafID = $descendants_with_genomes[$leafI];
			my $leaf_immediateDescendantI = $leave_to_immediate_descendant_i{$leafID}; die unless(defined $leaf_immediateDescendantI);
			my $leaf_immediateDescendantNodeID = $immediate_descendants[$leaf_immediateDescendantI];
			my $leaf_immediateDescendantNodeRank = $taxonomy->{$leaf_immediateDescendantNodeID}{rank};

			my @compare_against_taxonIDs;
			for(my $leafII = 0; $leafII <= $#descendants_with_genomes; $leafII++)
			{
				my $leafIID = $descendants_with_genomes[$leafII];
				my $leafII_immediateDescendantI = $leave_to_immediate_descendant_i{$leafIID}; die unless(defined $leafII_immediateDescendantI);	
				next if($leaf_immediateDescendantI == $leafII_immediateDescendantI);
				push(@compare_against_taxonIDs, $leafIID);
			}

			die unless(exists $taxonID_2_contigs{$leafID});
				
			my $thisLeaf_genomeSize = taxonID_get_genome_length($leafID, \%taxonID_2_contigs);
			my $sum_genomeSize_comparisons = 0;
			my @thisLeaf_contigs = keys %{$taxonID_2_contigs{$leafID}};
			my @comparison_contigs;
			foreach my $comparisonID (@compare_against_taxonIDs)
			{
				die unless(exists $taxonID_2_contigs{$comparisonID});
				push(@comparison_contigs, keys %{$taxonID_2_contigs{$comparisonID}});

				$sum_genomeSize_comparisons += taxonID_get_genome_length($comparisonID, \%taxonID_2_contigs);
			}
			# print "$leafID vs ", scalar(@compare_against_taxonIDs), " genomes.\n";
			# print "\tDescendants: ", join(", ", @immediate_descendants_ranks), "\n";
			$total_comparisons++;
			$total_storage += ($thisLeaf_genomeSize + $sum_genomeSize_comparisons);

			print JOBS join("\t", $nodeID, $taxonomy->{$nodeID}{rank}, $leaf_immediateDescendantNodeID, $leaf_immediateDescendantNodeRank, $leafID, $taxonomy->{$leafID}{rank}, $thisLeaf_genomeSize, scalar(@compare_against_taxonIDs), join(";", @compare_against_taxonIDs), $sum_genomeSize_comparisons, join(";", @thisLeaf_contigs), join(";", @comparison_contigs)), "\n";
		}
	}

	close(JOBS);
	print "\n";
	print "\nExpect to carry out ", $total_comparisons, " total comparisons.\n";
	print "\tRank heterogeneity in immediate descendants: $immediateDescendants_rank_heterogeneity \n";
	print "\tTotal storage: ", $total_storage/1e9, " GB.\n";

	print "\nJobs file: $outputFn_jobs \n";

	my $path_to_script = $FindBin::Bin.'/'.$FindBin::Script;


	my $qsub_file = $outputDir_computation . '/compute.qsub';
	open(QSUB, '>', $qsub_file) or die "Cannot open $qsub_file";
	print QSUB qq(#!/bin/bash
#\$ -t 1-${total_comparisons}
#\$ -l mem_free=1G
#\$ -N W_T_D
jobID=\$(expr \$SGE_TASK_ID - 1)
cd $FindBin::Bin
perl ${FindBin::Script} --dbDir $dbDir --taxonomyDir $taxonomyDir --mode doJob --jobI \$jobID
);
	close(QSUB);

	print "\n\nqsub $qsub_file\n\n";
	
elsif($mode eq 'collect')
{
	my %taxonID_2_contigs;
	my %contigLength;
	my $file_taxonGenomes = $dbDir . '/ref.fa.taxonGenomes';
	open(GENOMEINFO, '<', $file_taxonGenomes) or die "Cannot open $file_taxonGenomes";
	my $gI_headerLine = <GENOMEINFO>;
	chomp($gI_headerLine);
	my @gI_headerFields = split(/\t/, $gI_headerLine);
	while(<GENOMEINFO>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line);
		die unless($#line_fields == $#gI_headerFields);
		my %line_hash = (mesh @gI_headerFields, @line_fields);
		my $taxonID = $line_hash{taxonID};
		die if(exists $taxonID_2_contigs{$taxonID});
		my $contigs = $line_hash{contigs};
		my @components = split(/;/, $contigs);
		foreach my $component (@components)
		{
			my @p = split(/=/, $component);
			die unless(scalar(@p) == 2);
			die if(exists $taxonID_2_contigs{$taxonID}{$p[0]});
			$taxonID_2_contigs{$taxonID}{$p[0]} = $p[1];
			$contigLength{$p[0]} = $p[1];
		}
	}
	close(GENOMEINFO);
	my $taxonomy = taxTree::readTaxonomy($taxonomyDir);
		
	my %results_per_node;
	my %results_reads_per_node;
	my %results_reads_many_per_node;
	my $total_jobs = 0;
	my $total_jobs_ok = 0;
	my $total_jobs_reads_ok = 0;	
	my $total_jobs_reads_many_ok = 0;	
	open(JOBS, '<', $outputFn_jobs) or die "Canot open $outputFn_jobs";
	my $lineI = -1;
	open(SUMMARYINDIVIDUAL, '>', $outputFn_summary_individual) or die "Cannot open $outputFn_summary_individual";	
	while(<JOBS>)
	{
		my $line = $_;
		chomp($line);
		my @fields = split(/\t/, $line);
		my $nodeID = $fields[0];
		
		$lineI++;
		my $jobI = $lineI; 
		my $results_fn = $outputDir_results . '/' . $jobI . '.results';
		my $results_reads_fn = $outputDir_results . '/' . $jobI . '.results.reads';
		my $results_reads_many_fn = $outputDir_results . '/' . $jobI . '.results.reads.many';
		$total_jobs++;
		my $result;	
		# print join("\t", $results_fn, (-e $results_fn)), "\n";
		if(-e $results_fn)
		{
			open(R, '<', $results_fn) or die;
			while(<R>)
			{
				my $l = $_;
				chomp($l);
				# print $l, "\n";
				$l =~ s/[\n\r\s]//g;
				next if($_ =~ /^\s*$/);
				if($_ =~ /^([\d\.]+)$/)
				{
					$result = $1;
				}
				else
				{
					warn "Weird number  -" . $_ . "- in $results_fn";
				}
			}
			close(R);
		}
		
		if(-e $results_reads_fn)
		{
			my $S = 0;
			my %h;
			open(R, '<', $results_reads_fn) or die;
			while(<R>)
			{
				my $l = $_;
				chomp($l);
				next unless($l);
				my @fields = split(/\t/, $l);
				die unless(scalar(@fields) == 2);
				$h{$fields[0]} = $fields[1];
				$S += $fields[1];
			}
			close(R);
			if($S == 0)
			{
				# warn "Problem with $results_reads_fn";
			}
			else
			{
				foreach my $k (keys %h)
				{
					$h{$k} /= $S;
				}
				push(@{$results_reads_per_node{$nodeID}}, \%h);
				$total_jobs_reads_ok++;
			}
		}
		
		if(-e $results_reads_many_fn)
		{
			my $S = 0;
			my %S_by_readLength;
			my %h;
			open(R, '<', $results_reads_many_fn) or die;
			while(<R>)
			{
				my $l = $_;
				chomp($l);
				next unless($l);
				my @fields = split(/\t/, $l);
				die Dumper("Weird field count in line $. of file $results_reads_many_fn", \@fields) unless(scalar(@fields) == 3);
				$h{$fields[0]}{$fields[1]} = $fields[2];
				$S += $fields[2];
				$S_by_readLength{$fields[0]} += $fields[2];
			}
			close(R);
			if($S == 0)
			{
				# warn "Problem with $results_reads_fn";
			}
			else
			{
				foreach my $readLength (keys %h)
				{
					next if(not $S_by_readLength{$readLength});
					foreach my $k (keys %{$h{$readLength}})
					{
						$h{$readLength}{$k} /= $S_by_readLength{$readLength};
					}
					push(@{$results_reads_many_per_node{$readLength}{$nodeID}}, $h{$readLength});

				}

				$total_jobs_reads_many_ok++;
			}
		}
						
		
		if(defined $result)
		{
			$total_jobs_ok++;
			
			die unless(defined $taxonomy->{$nodeID});
			
			my $nodeRank = $taxonomy->{$nodeID}{rank};
			my $nodeName = join('; ', @{$taxonomy->{$nodeID}{names}});
			# warn "What $nodeRank $fields[0]" unless($nodeRank eq $fields[1]);
			
			push(@{$results_per_node{$nodeID}}, $result);
			print SUMMARYINDIVIDUAL join("\t", $nodeID, $nodeName, $nodeRank, $result), "\n";			
			
		}
	}
	close(JOBS);
	close(SUMMARYINDIVIDUAL);
	
	print "\nTotal jobs: $total_jobs\n";
	print "\tof which results: $total_jobs_ok \n";
	print "\tof which results for reads: $total_jobs_reads_ok \n";
	print "\tof which results for many reads: $total_jobs_reads_many_ok \n";
	
	open(RESULTS, '>', $outputFn_results) or die "Cannot open $outputFn_results";
	open(RESULTSREADS, '>', $outputFn_reads_results) or die "Cannot open $outputFn_reads_results";
	open(RESULTSREADSMANY, '>', $outputFn_reads_results_many) or die "Cannot open $outputFn_reads_results_many";
	my @taxonomy_leaves = taxTree::get_leave_ids($taxonomy);
	foreach my $leaveID (@taxonomy_leaves)
	{
		# print RESULTS join("\t", $leaveID, 1, $leaveID), "\n";
	}
	
	open(SUMMARY, '>', $outputFn_summary) or die "Cannot open $outputFn_summary";
	foreach my $nodeID (keys %results_per_node)
	{
		my $nodeRank = $taxonomy->{$nodeID}{rank};
		my $nodeName = join('; ', @{$taxonomy->{$nodeID}{names}});
		my @r = @{$results_per_node{$nodeID}};
		my $m = mean(@r);
		my $sd = sd(@r);
		
		print SUMMARY join("\t", $nodeID, $nodeName, $nodeRank, join(";", @r), $m, $sd), "\n";
		
		my @descendants = taxTree::descendants($taxonomy, $nodeID);
		my @descendants_with_genomes = grep {exists $taxonID_2_contigs{$_}} @descendants;
		
		print RESULTS join("\t", $nodeID, $m, join(';', @descendants_with_genomes)), "\n";
	}	
	
	foreach my $nodeID (keys %results_reads_per_node)
	{
		my $nodeRank = $taxonomy->{$nodeID}{rank};
		my $nodeName = join('; ', @{$taxonomy->{$nodeID}{names}});

		my @descendants = taxTree::descendants($taxonomy, $nodeID);
		my @descendants_with_genomes = grep {exists $taxonID_2_contigs{$_}} @descendants;
				
		my %combinedHistogram;
		my $S = 0;		
		foreach my $histogram (@{$results_reads_per_node{$nodeID}})
		{
			foreach my $k (keys %$histogram)
			{	
				$combinedHistogram{$k} += $histogram->{$k};
				$S += $histogram->{$k};
			}
		}
		
		my $firstKey = 1;
		foreach my $k (keys %combinedHistogram)
		{
			$combinedHistogram{$k} /= $S;
			my $string_source_genomes = '';
			if($firstKey)
			{
				$string_source_genomes = join(';', @descendants_with_genomes);
				$firstKey = 0;
			}
			print RESULTSREADS join("\t", $nodeID, $k, $combinedHistogram{$k}, $string_source_genomes), "\n";
		}
	}		
	
	foreach my $readLength (keys %results_reads_many_per_node)
	{
	
		foreach my $nodeID (keys %{$results_reads_many_per_node{$readLength}})
		{
			my $nodeRank = $taxonomy->{$nodeID}{rank};
			my $nodeName = join('; ', @{$taxonomy->{$nodeID}{names}});

			my @descendants = taxTree::descendants($taxonomy, $nodeID);
			my @descendants_with_genomes = grep {exists $taxonID_2_contigs{$_}} @descendants;
			my %combinedHistogram;
			my $S = 0;		
			foreach my $histogram (@{$results_reads_many_per_node{$readLength}{$nodeID}})
			{
				foreach my $k (keys %$histogram)
				{	
					$combinedHistogram{$k} += $histogram->{$k};
					$S += $histogram->{$k};
				}
			}
			
			my $firstKey = 1;				
			foreach my $k (keys %combinedHistogram)
			{
				$combinedHistogram{$k} /= $S;
				my $string_source_genomes = ($firstKey) ? join(';', @descendants_with_genomes) : '';
				my $string_rank = ($firstKey) ? $nodeRank : '';
				my $string_name = ($firstKey) ? $nodeName : '';
				$firstKey = 0;
				print RESULTSREADSMANY join("\t", $nodeID, $readLength, $k, $combinedHistogram{$k}, $string_source_genomes, $string_rank, $string_name), "\n";
			}
		}	
	}
	
	close(SUMMARY);  
	
	
	close(RESULTS);
	close(RESULTSREADS);
	close(RESULTSREADSMANY);
	
	print "Produced, amongst others:\n";
	print " - $outputFn_results \n";
	print " - $outputFn_reads_results \n";
}
elsif($mode eq 'doJob')
{
	die unless(defined $jobI);
	
	my ($contigs_A, $contigs_B);
	open(JOBS, '<', $outputFn_jobs) or die "Cannot open $outputFn_jobs";
	my $lineI = -1;
	while(<JOBS>)
	{
		my $line = $_;
		$lineI++;
		if($lineI == $jobI)
		{
			chomp($line);
			my @fields = split(/\t/, $line);
			die unless($#fields == 11);
			$contigs_A = $fields[10];
			$contigs_B = $fields[11];
		}
	}
	close(JOBS);
	die "Can't find job data for $jobI" unless(defined $contigs_A);

	my %contigs_A = map {$_ => 1} split(/;/, $contigs_A); die unless(scalar(keys %contigs_A));
	my %contigs_B = map {$_ => 1} split(/;/, $contigs_B); die unless(scalar(keys %contigs_B));

	my $dir_computation = $outputDir_computation . '/' . $jobI;
	mkdir($dir_computation);

	my $file_A = $dir_computation . '/A';
	my $file_A_reads = $dir_computation . '/A.reads';
	my $file_A_reads_many = $dir_computation . '/A.reads.many';
	my $file_B = $dir_computation . '/B';
	open(A, '>', $file_A) or die;
	open(B, '>', $file_B) or die;
	open(REF, '<', $file_ref) or die "Cannot open $file_ref";
	
	my $in_A = 0;
	my $in_B = 0;
	while(<REF>)
	{
		my $line = $_;
		chomp($line);
		if(substr($line, 0, 1) eq '>')
		{
			$in_A = 0;	
			$in_B = 0;
			my $contigID = substr($line, 1);
			if($contigs_A{$contigID})
			{
				$in_A = 1;
				$contigs_A{$contigID}--;
			}
			if($contigs_B{$contigID})
			{
				$in_B = 1;
				$contigs_B{$contigID}--;
			}
			die if($in_A and $in_B);
		}

		if($in_A)
		{
			print A $line, "\n";
		}

		if($in_B)
		{
			print B $line, "\n";
		}
	}
	close(REF);
	close(A);
	close(B);

	foreach my $contigID (keys %contigs_A)
	{
		die "Missed contig $contigID" if($contigs_A{$contigID});
	}

	foreach my $contigID (keys %contigs_B)
	{
		die "Missed contig $contigID" if($contigs_B{$contigID});
	}

	print "Two files:\n\t$file_A\n\t$file_B\n\n";

	my $mash_sketch_cmd = qq($mash_bin sketch -s 20000 $file_A);
	system($mash_sketch_cmd) and die "Mash sketching failed: $mash_sketch_cmd";
	my $mash_sketch = $file_A . '.msh';
	die unless(-e $mash_sketch);
	
	my $outputFn = $dir_computation . '.mashWithin';
	my $mash_comparison = qq($mash_bin screen $mash_sketch $file_B > $outputFn);
	system($mash_comparison) and die "Mash sketching failed: $mash_comparison";
	
	my $distance;
	open(OUTPUT, '<', $outputFn) or die "Cannot open $outputFn";
	while(<OUTPUT>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die "Problem with file $outputFn" if(defined $distance);
		die "Weird line $. in file $outputFn\n\n$line  " unless($line =~ /^([\d\.]+)\s/);
		$distance = $1;
	}
	close(OUTPUT);
	unless(defined $distance)
	{
		die "Distance problem -- $distance -- command $mash_comparison";
	}
	# print "Distance: $distance\n\n";

	my $results_fn = $outputDir_results . '/' . $jobI . '.results';
	open(RESULTS, '>', $results_fn) or die "Cannot open $results_fn";
	print RESULTS $distance, "\n";
	close(RESULTS);

	print "Produced $results_fn\n\n";
	
	{
		open(READS, '>', $file_A_reads) or die "Cannot open $file_A_reads";
		my $A_contigs_href = readFASTA($file_A);
		my $read_start_positions  = 0;
		foreach my $contigID (keys %$A_contigs_href)
		{
			my $contigSequence = $A_contigs_href->{$contigID};
			die unless(defined $contigSequence);
			for(my $posI = 0; $posI < length($contigSequence); $posI += $readSimDelta)
			{
				my $lastPos = $posI + $readSimSize - 1;
				if($lastPos < length($contigSequence))
				{
					$read_start_positions++;
				}
			}  
		}	
		
		my $start_rate = 1;
		if($read_start_positions > $target_max_simulatedChunks)
		{
			$start_rate = $target_max_simulatedChunks / $read_start_positions;
			print "(One read length) Adjusted start rate to $start_rate (eligible start positions: $read_start_positions, want $target_max_simulatedChunks)\n";
		}
		die unless(($start_rate >= 0) and ($start_rate <= 1));
		
		my $n_reads = 0;
		foreach my $contigID (keys %$A_contigs_href)
		{
			my $contigSequence = $A_contigs_href->{$contigID};
			POS: for(my $posI = 0; $posI < length($contigSequence); $posI += $readSimDelta)
			{
				my $lastPos = $posI + $readSimSize - 1;
				if($lastPos < length($contigSequence))
				{
					if($start_rate != 1)
					{
						next POS if(rand(1) > $start_rate);
					}
					
					$n_reads++;
					my $readID = 'read' . $n_reads;
					my $S = substr($contigSequence, $posI, $readSimSize);
					print READS '>', $readID, "\n";
					print READS $S, "\n";
				}
			}  
		}
		close(READS);
		
		print "Reads file $file_A_reads";
		
		my $outputFn_MetaMap = $dir_computation . '.MetaMap';
		# my $MetaMap_cmd = qq($metamap_bin -s $file_B -q $file_A_reads -m $readSimSize -o $outputFn_MetaMap); # todo reinstate
		my $MetaMap_cmd = qq($metamap_bin -s $file_B -q $file_A_reads -m 2000 -o $outputFn_MetaMap);
		system($MetaMap_cmd) and die "MetaMap command $MetaMap_cmd failed";
		print "Executed command $MetaMap_cmd \n";
		
		my $n_read_alignments = 0;
		my %read_alignment_histogram;
		my $currentReadID = '';
		my @currentReadLines;
		my $processAlignments_oneRead = sub {
			my $bestIdentity;
			foreach my $line (@currentReadLines)
			{
				my @fields = split(/ /, $line);
				die Dumper($fields[0], $currentReadID) unless($fields[0] eq $currentReadID);
				
				my $readID = $fields[0];
				my $identity = $fields[9];
				die unless(($identity >= 0) and ($identity <= 100));
				if((not defined $bestIdentity) or ($bestIdentity < $identity))
				{
					$bestIdentity = $identity;
				}
			}
			$bestIdentity = int($bestIdentity + 0.5);
			$read_alignment_histogram{$bestIdentity}++;
			$n_read_alignments++;
		};	
		
		open(MetaMapOUTPUT, '<', $outputFn_MetaMap) or die "Cannot open $outputFn_MetaMap";
		while(<MetaMapOUTPUT>)
		{
			chomp;
			next unless($_);
			die "Weird input" unless($_ =~ /^(.+?) /);
			my $readID = $1;
			if($currentReadID ne $readID)
			{
				if(@currentReadLines)
				{
					$processAlignments_oneRead->();
				}
				$currentReadID = $readID;
				@currentReadLines = ();
			}
			push(@currentReadLines, $_);
			# last if ($processed_reads > 10); 
		}
		if(@currentReadLines)
		{
			$processAlignments_oneRead->();
		}	
		
		close(MetaMapOUTPUT);

		my $n_missing_reads = $n_reads - $n_read_alignments;
		die unless($n_missing_reads >= 0);
		$read_alignment_histogram{0} += $n_missing_reads;
		
		my $results_fn_reads = $outputDir_results . '/' . $jobI . '.results.reads';
		
		open(RESULTSMetaMap, '>', $results_fn_reads) or die "Cannot open $results_fn_reads";
		foreach my $k (sort {$a <=> $b} keys %read_alignment_histogram)
		{
			print RESULTSMetaMap join("\t", $k, $read_alignment_histogram{$k}), "\n";
		}
		close(RESULTSMetaMap);
		
		print "Produced $results_fn_reads\n";
	}
	
	{
		my $results_fn_reads_many = $outputDir_results . '/' . $jobI . '.results.reads.many';
		open(F, '>', $results_fn_reads_many) or die "Cannot open $results_fn_reads_many";
		close(F);

		for(my $chunkLength = $readSimSizeFrom; $chunkLength <= $readSimSizeTo; $chunkLength += $readSimSizeStep)
		{
			print "Chunk length $chunkLength\n";
			open(READS, '>', $file_A_reads_many) or die "Cannot open $file_A_reads_many";
			my $A_contigs_href = readFASTA($file_A);
			
			my $read_start_positions  = 0;
			foreach my $contigID (keys %$A_contigs_href)
			{
				my $contigSequence = $A_contigs_href->{$contigID};
				die unless(defined $contigSequence);
				for(my $posI = 0; $posI < length($contigSequence); $posI += $readSimDelta)
				{
					my $lastPos = $posI + $chunkLength - 1;
					if($lastPos < length($contigSequence))
					{
						$read_start_positions++;
					}
				}  
			}	
			
			my $start_rate = 1;
			if($read_start_positions > $target_max_simulatedChunks)
			{
				$start_rate = $target_max_simulatedChunks / $read_start_positions;
				print "(Many read lengths; $chunkLength) Adjusted start rate to $start_rate (eligible start positions: $read_start_positions, want $target_max_simulatedChunks)\n";
			}
			die unless(($start_rate >= 0) and ($start_rate <= 1));
			
		
			
			my $n_reads = 0;
			foreach my $contigID (keys %$A_contigs_href)
			{
				my $contigSequence = $A_contigs_href->{$contigID};
				POS: for(my $posI = 0; $posI < length($contigSequence); $posI += $readSimDelta)
				{
					my $lastPos = $posI + $chunkLength - 1;
					if($lastPos < length($contigSequence))
					{
					
						if($start_rate != 1)
						{
							next POS if(rand(1) > $start_rate);
						}
						
						$n_reads++;
						my $readID = 'read' . $n_reads;
						my $S = substr($contigSequence, $posI, $chunkLength);
						print READS '>', $readID, "\n";
						print READS $S, "\n";
					}
				}  
			}
			close(READS);
			
			print "Reads file $file_A_reads_many";
			
			my $outputFn_MetaMap = $dir_computation . '.MetaMap.many';
			my $MetaMap_cmd = qq($metamap_bin -s $file_B -q $file_A_reads_many -m $chunkLength -o $outputFn_MetaMap);
			system($MetaMap_cmd) and die "MetaMap command $MetaMap_cmd failed";
			print "Executed command $MetaMap_cmd \n";
			
			my $n_read_alignments = 0;
			my %read_alignment_histogram;
			my $currentReadID = '';
			my @currentReadLines;
			my $processAlignments_oneRead = sub {
				my $bestIdentity;
				foreach my $line (@currentReadLines)
				{
					my @fields = split(/ /, $line);
					die Dumper($fields[0], $currentReadID) unless($fields[0] eq $currentReadID);
					
					my $readID = $fields[0];
					my $identity = $fields[9];
					die unless(($identity >= 0) and ($identity <= 100));
					if((not defined $bestIdentity) or ($bestIdentity < $identity))
					{
						$bestIdentity = $identity;
					}
				}
				$bestIdentity = int($bestIdentity + 0.5);
				$read_alignment_histogram{$bestIdentity}++;
				$n_read_alignments++;
			};	
			
			open(MetaMapOUTPUT, '<', $outputFn_MetaMap) or die "Cannot open $outputFn_MetaMap";
			while(<MetaMapOUTPUT>)
			{
				chomp;
				next unless($_);
				die "Weird input" unless($_ =~ /^(.+?) /);
				my $readID = $1;
				if($currentReadID ne $readID)
				{
					if(@currentReadLines)
					{
						$processAlignments_oneRead->();
					}
					$currentReadID = $readID;
					@currentReadLines = ();
				}
				push(@currentReadLines, $_);
				# last if ($processed_reads > 10); 
			}
			if(@currentReadLines)
			{
				$processAlignments_oneRead->();
			}	
			
			close(MetaMapOUTPUT);

			my $n_missing_reads = $n_reads - $n_read_alignments;
			die unless($n_missing_reads >= 0);
			$read_alignment_histogram{0} += $n_missing_reads;
			
			open(RESULTSMetaMap, '>>', $results_fn_reads_many) or die "Cannot open $results_fn_reads_many";
			foreach my $k (sort {$a <=> $b} keys %read_alignment_histogram)
			{
				print RESULTSMetaMap join("\t", $chunkLength, $k, $read_alignment_histogram{$k}), "\n";
			}
			close(RESULTSMetaMap);
		}
		
		print "Produced $results_fn_reads_many\n";
	}
		
	unlink($file_A);
	unlink($file_B);
	unlink($mash_sketch);
	unlink($outputFn);
	system('rm -r ' . $dir_computation);
}
else
{
	die "Unknown mode: $mode";
}

