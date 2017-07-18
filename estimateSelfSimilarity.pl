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
my $mode;
my $jobI;
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
my $outputFn_qsub = $outputDir . '/jobs.qsub';
my $outputFn_reads_results_many = $outputDir . '/results.reads.many.byNode';
my $finalResultsFile = $DB . '/selfSimilarities.txt';

if(not $mode)
{
	die "Please specify --mode; e.g.: prepareFromScratch";
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
	foreach my $nodeID (@nodesForAttachment)
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
			my @targetLeaf_compareAgainst = @{$node_computation->[2]};
			my @targetLeaf_contigs = keys %{$taxonID_2_contigs{$targetLeafID}};
			
			my $targetLeaf_compareAgainst_combinedSize = 0;
			my @compareAgainst_contigs;
			my %ranks_compareAgainst;
			foreach my $taxonID (@targetLeaf_compareAgainst)
			{
				$targetLeaf_compareAgainst_combinedSize += Util::getGenomeLength($taxonID, \%taxonID_2_contigs, \%contigLength);
				push(@compareAgainst_contigs, keys %{$taxonID_2_contigs{$taxonID}});
				$ranks_compareAgainst{$taxonomy->{$taxonID}{rank}}++;
			}
			
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
elsif($mode eq 'doJobI')
{
	die "Please specify --jobI" unless(defined $jobI);
	
	my ($A_taxonID, $B_taxonIDs, $contigs_A, $contigs_B);
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
			die unless($#fields == 12);
			$A_taxonID = $fields[4];
			$B_taxonIDs = $fields[8];			
			$contigs_A = $fields[10];
			$contigs_B = $fields[11];
		}
	}
	close(JOBS);
	die "Can't find job data for $jobI" unless(defined $contigs_A);

	print "Job $jobI";
	print "\tMapping from: $A_taxonID\n";
	print "\tMapping to  : $B_taxonIDs\n";
	
	my %contigs_A = map {$_ => 1} split(/;/, $contigs_A); die unless(scalar(keys %contigs_A));
	my %contigs_B = map {$_ => 1} split(/;/, $contigs_B); die unless(scalar(keys %contigs_B));
	
	my $dir_computation = $outputDir_computation . '/' . $jobI;
	mkdir($dir_computation);

	my $file_A = $dir_computation . '/A';
	my $file_A_reads_many = $dir_computation . '/A.reads.many';
	my $file_B = $dir_computation . '/B';
	my $outputFn_MetaMap = $dir_computation . '.MetaMap.many';
	
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

	my $results_fn_reads_many = get_results_file_for_jobI($jobI);
	
	# make sure results file is writable and empty
	open(F, '>', $results_fn_reads_many) or die "Cannot open $results_fn_reads_many";
	close(F);

	for(my $chunkLength = $readSimSizeFrom; $chunkLength <= $readSimSizeTo; $chunkLength += $readSimSizeStep)
	{
		print "\tSimulate chunk length $chunkLength\n";
		open(READS, '>', $file_A_reads_many) or die "Cannot open $file_A_reads_many";
		my $A_contigs_href = Util::readFASTA($file_A);
		
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
			print "\t\t\t(Read lengths; $chunkLength) Adjusted start rate to $start_rate (eligible start positions: $read_start_positions, want $target_max_simulatedChunks)\n";
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
					die unless(length($S) == $chunkLength);
					print READS '>', $readID, "\n";
					print READS $S, "\n";
				}
			}  
		}
		close(READS);
		
		print "\t\tReads file $file_A_reads_many with $n_reads reads";
				
		my $MetaMap_cmd = qq($metamap_bin mapDirectly -r $file_B -q $file_A_reads_many -m $readSimSizeFrom -o $outputFn_MetaMap --pi 80);
		print "\t\tExecuted command $MetaMap_cmd\n";		
		system($MetaMap_cmd) and die "MetaMap command $MetaMap_cmd failed";
		
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
	
	print "--------\nJob $jobI done. Produced $results_fn_reads_many.\n";
		
	unlink($file_A);
	unlink($file_A_reads_many);
	unlink($file_B);
	unlink($outputFn_MetaMap);
	system('rm ' . $dir_computation . '/*') and die "Cannot delete $dir_computation (I)";	
	system('rm -r ' . $dir_computation) and die "Cannot delete $dir_computation (II)";	
}
elsif($mode eq 'collect')
{
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
	open(JOBS, '<', $outputFn_jobs) or die "Canot open $outputFn_jobs";
	my $lineI = -1;
	while(<JOBS>)
	{
		my $line = $_;
		chomp($line);
		my @fields = split(/\t/, $line);
		my $nodeID = $fields[0];
		
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


else
{
	die "Unknown mode";
}


sub get_results_file_for_jobI
{
	my $jobI = shift;
	return  $outputDir_results . '/' . $jobI . '.results.reads.many';;
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

