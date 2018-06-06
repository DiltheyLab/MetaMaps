package validation;

use strict;
use Data::Dumper;
use List::Util qw/all sum min max/;
use List::MoreUtils qw/mesh/;
use taxTree;
use Statistics::Basic qw/correlation/;
use Storable qw/dclone store retrieve/;

my @evaluateAccuracyAtLevels = qw/species genus family superkingdom/;
{
	my %_knowRank = map {$_ => 1} taxTree::getRelevantRanks();
	die unless(all {$_knowRank{$_}} @evaluateAccuracyAtLevels);
}

sub getEvaluationLevels
{
	return @evaluateAccuracyAtLevels;
}

my %_cache_prepare_masterTaxonomy_withX_taxonomies;
sub prepare_masterTaxonomy_withX
{
	my $masterTaxonomy_dir = shift; 
	my $taxonomyWithX = shift;

	my $masterTaxonomy;
	if($_cache_prepare_masterTaxonomy_withX_taxonomies{$masterTaxonomy_dir})
	{
		$masterTaxonomy = $_cache_prepare_masterTaxonomy_withX_taxonomies{$masterTaxonomy_dir};
	}
	else
	{
		$masterTaxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
		$_cache_prepare_masterTaxonomy_withX_taxonomies{$masterTaxonomy_dir} = $masterTaxonomy;
	}
	
	my $masterTaxonomy_merged = taxTree::readMerged($masterTaxonomy_dir);
	
	my $extendedMaster = taxTree::cloneTaxonomy_integrateX($masterTaxonomy, $masterTaxonomy_merged, $taxonomyWithX);

	return ($extendedMaster, $masterTaxonomy_merged);
}

sub readTruthFileReads
{
	my $masterTaxonomy = shift;
	my $masterTaxonomy_merged = shift;
	my $file = shift;
	
	die unless(defined $file);
	my %truth_raw_reads;
	
	open(T, '<', $file) or die "Cannot open $file";
	while(<T>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		my $readID = $f[0];
		my $taxonID_original = $f[1];
		
		# get the taxon ID in the master taxonomy
		my $taxonID_master;
		
		if($taxonID_original ne 0)
		{
			$taxonID_master = taxTree::findCurrentNodeID($masterTaxonomy, $masterTaxonomy_merged, $taxonID_original);
			die "Unknown ID '$taxonID_master' (from $taxonID_original)" unless(exists $masterTaxonomy->{$taxonID_master});
		}
		else
		{
			$taxonID_master = 0;
		}	
		
		die "Read ID $readID has multiple truth definitions in $file" if(defined $truth_raw_reads{$readID});
		
		$truth_raw_reads{$readID} = $taxonID_master;
	}
	close(T);	
	
	return \%truth_raw_reads;
}

sub translateID2ReducedTaxonomy
{
	my $masterTaxonomy = shift;
	my $reducedTaxonomy = shift;
	my $taxonID = shift;
	
	die unless(defined $taxonID);
	die unless(exists $masterTaxonomy->{$taxonID});
	
	my @nodes_to_consider = ($taxonID, taxTree::get_ancestors($masterTaxonomy, $taxonID));
	
	foreach my $nodeID (@nodes_to_consider)
	{
		if(exists $reducedTaxonomy->{$nodeID})
		{
			return $nodeID;
		}
	}

	checkTaxonomyIsProperReduction();
	
	die "Can't translate taxon $taxonID into the space of reduced taxonomy.";
}

sub translateReadsTruthToReducedTaxonomy
{
	my $masterTaxonomy = shift;
	my $reducedTaxonomy = shift;
	my $reads_href = shift;
	
	checkTaxonomyIsProperReduction();

	my %forReturn_reads;
	my %taxonID_translation;	
	$taxonID_translation{'0'} = '0';
	$taxonID_translation{0} = '0';
	
	foreach my $readID (keys %$reads_href)
	{
		my $taxonID = $reads_href->{$readID};
		die unless(($taxonID eq '0') or (exists $masterTaxonomy->{$taxonID}));
		unless(defined $taxonID_translation{$taxonID})
		{
			$taxonID_translation{$taxonID} = translateID2ReducedTaxonomy($masterTaxonomy, $reducedTaxonomy, $taxonID);
		}
		$forReturn_reads{$readID} = $taxonID_translation{$taxonID};
		die unless(defined $forReturn_reads{$readID});
	}
	
	return \%forReturn_reads;
}

sub checkTaxonomyIsProperReduction
{
	my $masterTaxonomy = shift;
	my $reducedTaxonomy = shift;
	
	
	unless(all {exists $masterTaxonomy->{$_}} keys %$reducedTaxonomy)
	{
		my @missing = grep {not exists $masterTaxonomy->{$_}} keys %$reducedTaxonomy;
		die Dumper("Reduced taxonomy doesn't seem to be a proper subset of full taxonomy!", \@missing);
	}	
}

# returns hash:
# taxonID -> [
#	total attached reads
#	newly attached reads
#	total newly attached reads

sub truthReadsTree
{
	my $taxonomy_href = shift;
	my $truthReads_href = shift;
	my $taxonIDs_in_direct_truth_href = shift;

	my %counts_per_taxonID;
	foreach my $readID (keys %$truthReads_href)
	{
		my $taxonID = $truthReads_href->{$readID};
		$counts_per_taxonID{$taxonID}++;
	}
	
	my %treeCounts;
	my $n_indirect = 0;
	foreach my $taxonID (keys %counts_per_taxonID)
	{
		my $n_reads = $counts_per_taxonID{$taxonID};
		my @ancestor_nodes = taxTree::get_ancestors($taxonomy_href, $taxonID);
		foreach my $nodeID ($taxonID, @ancestor_nodes)
		{
			unless(defined $treeCounts{$nodeID})
			{
				$treeCounts{$nodeID} = [0, 0, 0];
			}
			$treeCounts{$nodeID}[0] += $n_reads;
		}
		unless($taxonIDs_in_direct_truth_href->{$taxonID})
		{
			$n_indirect += $n_reads;
		}
		
	}
	
	foreach my $taxonID (keys %treeCounts)
	{
		unless($taxonIDs_in_direct_truth_href->{$taxonID})
		{
			my $reads_from_below = 0;
			foreach my $childID (@{$taxonomy_href->{$taxonID}{children}})
			{
				if(exists $treeCounts{$childID})
				{
					$reads_from_below += $treeCounts{$childID}[0];
				}
			}
			die Dumper($reads_from_below, $treeCounts{$taxonID}[0]) unless($reads_from_below <= $treeCounts{$taxonID}[0]);
			my $new_reads = $treeCounts{$taxonID}[0] - $reads_from_below;
			die unless($new_reads >= 0);
			$treeCounts{$taxonID}[1] = $new_reads;
			$treeCounts{$taxonID}[2] += $new_reads;
			
			my @ancestor_nodes = taxTree::get_ancestors($taxonomy_href, $taxonID);
			foreach my $ancestorNodeID (@ancestor_nodes)
			{
				$treeCounts{$ancestorNodeID}[2] += $new_reads;
			}		
		}
	}
	
	die Dumper("Read counts don't add up", $treeCounts{1}, scalar(keys %$truthReads_href)) unless(($treeCounts{1}[0] + $treeCounts{1}[1]) == scalar(keys %$truthReads_href));
	die unless($treeCounts{1}[2] == $n_indirect);
	foreach my $treeNode (keys %treeCounts)
	{
		die unless($treeCounts{$treeNode}[1] <= $treeCounts{$treeNode}[2]);
		die unless($treeCounts{$treeNode}[2] <= $treeCounts{$treeNode}[0]);
		
		my $children_sum_total = 0;
		my $children_sum_novel = 0;
		foreach my $childID (@{$taxonomy_href->{$treeNode}{children}})
		{
			if(exists $treeCounts{$childID})
			{
				$children_sum_total += $treeCounts{$childID}[0];
				$children_sum_novel += $treeCounts{$childID}[2];
			}
		}
		my $expected_total_sum = $children_sum_total + ((exists $counts_per_taxonID{$treeNode}) ? $counts_per_taxonID{$treeNode} : 0);
		die unless($treeCounts{$treeNode}[0] == $expected_total_sum);
		
		my $expected_novel_sum_I = ($taxonIDs_in_direct_truth_href->{$treeNode}) ? 0 : ((exists $counts_per_taxonID{$treeNode}) ? $counts_per_taxonID{$treeNode} : 0);
		my $expected_novel_sum_II = $children_sum_novel + $expected_novel_sum_I;
		
		die unless($treeCounts{$treeNode}[1] == $expected_novel_sum_I);
		die unless($treeCounts{$treeNode}[2] == $expected_novel_sum_II);
		
	}
	return \%treeCounts;
}

sub truthReadsToTruthSummary
{
	my $taxonomy = shift;
	my $truthReads_href = shift;
	my $mappableTaxonIDs_href = shift;
	die unless(defined $mappableTaxonIDs_href);
	
	my %taxonID_translation;	
	my %truth_allReads;
	foreach my $readID (keys %$truthReads_href)
	{
		my $taxonID = $truthReads_href->{$readID};
		unless(defined $taxonID_translation{$taxonID})
		{
			$taxonID_translation{$taxonID} = getAllRanksForTaxon_withUnclassified($taxonomy, $taxonID, $mappableTaxonIDs_href);
			$taxonID_translation{$taxonID}{definedAndHypotheticalGenomes} = $taxonID;
			$taxonID_translation{$taxonID}{absolute} = $taxonID;
		}
		
		foreach my $rank (keys %{$taxonID_translation{$taxonID}})
		{
			my $v = $taxonID_translation{$taxonID}{$rank};
			$truth_allReads{$rank}{$v}++;
		}
	}
	
		
	foreach my $rank (keys %truth_allReads)
	{
		foreach my $taxonID (keys %{$truth_allReads{$rank}})
		{
			$truth_allReads{$rank}{$taxonID} /= scalar(keys %$truthReads_href);
		}
	}	
	
	return \%truth_allReads;
}

sub getAllRanksForTaxon_withUnclassified
{
	my $taxonomy = shift;
	my $taxonID = shift;
	my $mappableTaxonIDs_href = shift;
	die unless(defined $mappableTaxonIDs_href);
	
	my %forReturn = map {$_ => 'Unclassified'} @evaluateAccuracyAtLevels;
	if(defined $mappableTaxonIDs_href->{$taxonID})
	{
		$forReturn{definedGenomes} = $taxonID; 
	}
	else
	{
		$forReturn{definedGenomes} = 'Unclassified';
	}
	return \%forReturn if($taxonID eq '0');
	
	die unless(defined $taxonID);
	die Dumper("Undefined taxon ID $taxonID", $taxonID, ($taxonID eq '0'), ($taxonID == 0)) unless(defined $taxonomy->{$taxonID});
	
	my @nodes_to_consider = ($taxonID, taxTree::get_ancestors($taxonomy, $taxonID));
	
	my $firstRankAssigned;
	my $inDefinedRanks = 0;
	my %sawRank;
	foreach my $nodeID (@nodes_to_consider)
	{
		my $rank = $taxonomy->{$nodeID}{rank};
		die unless(defined $rank);
		$sawRank{$rank}++;
		
		
		if(exists $forReturn{$rank})
		{
			$forReturn{$rank} = $nodeID;
			$firstRankAssigned = $rank if not(defined $firstRankAssigned);
		}
	}
	
	die Dumper("Unexpected behaviour", $taxonID) if(($sawRank{'subspecies'} and not $sawRank{'species'}) or ($sawRank{'strain'} and not $sawRank{'species'}));
	
	if(defined $firstRankAssigned)
	{
		# we don't set stuff to undefined / NotLabelledAtLevel anymore
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
					$forReturn{$rank} = 'Unclassified'; # previously: NotLabelledAtLevel
				}
			}
		}
	}		
	
	return \%forReturn;
}

# sub getRank_phylogeny
# {
	# my $fullTaxonomy = shift;
	# my $reducedTaxonomy = shift;
	# my $taxonID = shift;
	
	# my %forReturn = map {$_ => 'Unclassified'} @evaluateAccuracyAtLevels;
	
	# unless(all {exists $fullTaxonomy->{$_}} keys %$reducedTaxonomy)
	# {
		# my @missing = grep {not exists $fullTaxonomy->{$_}} keys %$reducedTaxonomy;
		# die Dumper("Reduced taxonomy doesn't seem to be a proper subset of full taxonomy!", \@missing);
	# }
	# die unless(defined $fullTaxonomy->{$taxonID});
	
	# my @nodes_to_consider = ($taxonID, taxTree::get_ancestors($fullTaxonomy, $taxonID));
	
	# my $inTaxonomicAgreement = 0;
	# my $firstRankAssigned;
	# foreach my $nodeID (@nodes_to_consider)
	# {
		# my $rank = $fullTaxonomy->{$nodeID}{rank};
		# die unless(defined $rank);
		# if(exists $forReturn{$rank})
		# {
			# if(exists $reducedTaxonomy->{$nodeID})
			# {
				# $forReturn{$rank} = $nodeID;
				# $inTaxonomicAgreement = 1;
				# $firstRankAssigned = $rank;
			# }
			# else
			# {
				# $forReturn{$rank} = 'Unclassified';
				# die if($inTaxonomicAgreement);
			# }
		# }
	# }
	
	# if(defined $firstRankAssigned)
	# {
		# my $setToUndefined = 0;
		# foreach my $rank (taxTree::getRelevantRanks())
		# {
			# next unless(exists $forReturn{$rank});
			# if($rank eq $firstRankAssigned)
			# {
				# $setToUndefined = 1;
				# die if($forReturn{$rank} eq 'Unclassified');
			# }
			# else
			# {
				# if($setToUndefined and ($forReturn{$rank} eq 'Unclassified'))
				# {
					# $forReturn{$rank} = 'Undefined';
				# }
			# }
		# }
	# }
	
	# return \%forReturn;
# }


sub readLevelComparison
{
	my $masterTaxonomy = shift;
	my $reads_truth_absolute = shift;
	my $reads_truth_mappingDB = shift;
	my $reads_inferred = shift;
	my $label = shift;
	my $external_reads_correct = shift;
	my $external_reads_correct_byLevel = shift;
	my $external_reads_correct_byLevel_byLength = shift;
	my $mappableTaxonIDs = shift;
	my $readLengths_href = shift;
	my $external_n_reads_unknownStats_byLevel_href = shift;
	
	die unless(defined $label);
	die unless(defined $mappableTaxonIDs);
	die unless(defined $readLengths_href);
	
	my @readIDs = keys %$reads_truth_absolute;
	die "readLevelComparison(..): We don't have absolute truth for some reads" unless(all {exists $reads_truth_absolute->{$_}} keys %$reads_inferred);
	
	#unless((scalar(@readIDs) == scalar(keys %$reads_inferred)) and (all {exists $reads_inferred->{$_}} @readIDs))
	#{
		#die Dumper("Read ID problem - truth and inference sets are not congruent.", $label, [grep {not exists $reads_inferred->{$_}} @readIDs]);
	#}
	
	die unless(all {($_ eq '0') or (exists $masterTaxonomy->{$_})} values %$reads_truth_absolute);
	die unless(all {($_ eq '0') or (exists $masterTaxonomy->{$_})} values %$reads_inferred);

	# a 'lightning' is the way from a taxon ID to the top of the tree, and back
	# ... setting levels below the ID to 'unclassified'
	my %_getLightning_cache;
	my $getLightning = sub {
		my $taxonID = shift;
		if(exists $_getLightning_cache{$taxonID})
		{
			return $_getLightning_cache{$taxonID};
		}
		else
		{
			my $lightning = getAllRanksForTaxon_withUnclassified($masterTaxonomy, $taxonID, $mappableTaxonIDs);
			$_getLightning_cache{$taxonID} = $lightning;
			return $lightning;
		}
	};
	
	my $get_read_categories = sub {
		my $readID = shift;
		
		my @categories = ('ALL');
		if(($reads_truth_mappingDB->{$readID} eq $reads_truth_absolute->{$readID}))
		{
			if($reads_truth_absolute->{$readID})
			{
				push(@categories, 'truthLeafInDB');
			}
			else
			{
				push(@categories, 'truthUnclassified');
			}
		}
		else
		{
			push(@categories, 'novel');
			my $lightning_truth_inDB = $getLightning->($reads_truth_mappingDB->{$readID});
			my $shouldBeAssignedTo;
			RANK: foreach my $rank (@evaluateAccuracyAtLevels)
			{
				die unless(defined $lightning_truth_inDB->{$rank});
				die if($lightning_truth_inDB->{$rank} eq 'NotLabelledAtLevel');
				if(($lightning_truth_inDB->{$rank} ne 'Unclassified') and ($lightning_truth_inDB->{$rank} ne 'NotLabelledAtLevel'))
				{
					die if($lightning_truth_inDB->{$rank} eq 'Undefined');
					$shouldBeAssignedTo = $rank;
					last RANK;
				}
			}
			die unless(defined $shouldBeAssignedTo);
			push(@categories, 'novel_to_' . $shouldBeAssignedTo);
		}
		
		if($reads_truth_absolute->{$readID})
		{
			# push(@categories, 'truthClassified');
		}	
		
		return @categories;
	};
	
	my $get_read_assignedToLevel = sub {
		my $readID = shift;
		my $taxonID_inferred =  $reads_inferred->{$readID};
		die unless(defined $taxonID_inferred);
		
		my $lightning_inferred = $getLightning->($taxonID_inferred);
		
		if($taxonID_inferred eq '0')
		{
			return 'Unclassified';
		}
		my $assignedToRank;
		RANK: foreach my $rank (@evaluateAccuracyAtLevels)
		{
			die unless(defined $lightning_inferred->{$rank});
			if(($lightning_inferred->{$rank} ne 'Unclassified') and ($lightning_inferred->{$rank} ne 'NotLabelledAtLevel'))
			{
				$assignedToRank = $rank;
				last RANK;
			}
		}
		unless(defined $assignedToRank)
		{
			$assignedToRank = '>' . $evaluateAccuracyAtLevels[$#evaluateAccuracyAtLevels];
		}
		die unless(defined $assignedToRank);
		# die "Cannot assign rank of ID $readID -- assigned to $taxonID_inferred" unless(defined $assignedToRank);#
		return $assignedToRank;
	};
	
	my $getLengthBinForLength = sub {
		my $length = shift;
		return int(($length / 500) + 0.5) * 500;
	};
	die unless($getLengthBinForLength->(200) == 0);
	die unless($getLengthBinForLength->(300) == 500);
	die unless($getLengthBinForLength->(500) == 500);
	die unless($getLengthBinForLength->(700) == 500);
	die unless($getLengthBinForLength->(800) == 1000);
	die unless($getLengthBinForLength->(1000) == 1000);

	my %n_reads_correct;
	my %n_reads_correct_byLevel;
	my %n_reads_correct_byLevel_byLength;
	my %n_reads_unknownStats_byLevel;
	# my %taxonID_across_ranks;
	foreach my $readID (@readIDs)
	{
		my @read_categories = $get_read_categories->($readID);
		my %_read_categories = map {$_ => 1} @read_categories;
		
		die unless(defined $readLengths_href->{$readID});
		my $readLengthBin = $getLengthBinForLength->($readLengths_href->{$readID});
		
		foreach my $category (@read_categories)
		{
			unless(defined $n_reads_correct{$category}{N})
			{
				$n_reads_correct{$category}{missing} = 0;
				$n_reads_correct{$category}{N} = 0;
				$n_reads_correct{$category}{correct} = 0;
			}	
			
			unless(defined $n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{N})
			{
				$n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{missing} = 0;
				$n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{N} = 0;
				$n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{correct} = 0;			
			}	
			
			unless(defined $n_reads_unknownStats_byLevel{$category}{'genome'}{N_unclassified_should0_is0})
			{
				$n_reads_unknownStats_byLevel{$category}{'genome'}{N_unclassified_should0_is0} = 0;
				$n_reads_unknownStats_byLevel{$category}{'genome'}{N_unclassified_should0_is1} = 0;
				$n_reads_unknownStats_byLevel{$category}{'genome'}{N_unclassified_should1_is0} = 0;
				$n_reads_unknownStats_byLevel{$category}{'genome'}{N_unclassified_should1_is1} = 0;			
			}

			foreach my $level ('absolute', @evaluateAccuracyAtLevels)
			{
				unless(defined $n_reads_correct_byLevel{$category}{$level}{N})
				{
					$n_reads_correct_byLevel{$category}{$level}{missing} = 0;
					$n_reads_correct_byLevel{$category}{$level}{N} = 0;
					$n_reads_correct_byLevel{$category}{$level}{N_truthDefined} = 0;
					$n_reads_correct_byLevel{$category}{$level}{correct} = 0;
					$n_reads_correct_byLevel{$category}{$level}{correct_exactlyAtLevel} = 0;
					$n_reads_correct_byLevel{$category}{$level}{correct_truthDefined} = 0;	
					
					$n_reads_unknownStats_byLevel{$category}{$level}{N_unclassified_should0_is0} = 0;
					$n_reads_unknownStats_byLevel{$category}{$level}{N_unclassified_should0_is1} = 0;
					$n_reads_unknownStats_byLevel{$category}{$level}{N_unclassified_should1_is0} = 0;
					$n_reads_unknownStats_byLevel{$category}{$level}{N_unclassified_should1_is1} = 0;				
				}
				
				unless(defined $n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{N})
				{
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{missing} = 0;
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{N} = 0;
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{correct} = 0;
				}
			}
		}
		
		my $trueTaxonID_inUsedDB = $reads_truth_mappingDB->{$readID};		
		my $lightning_truth = $getLightning->($trueTaxonID_inUsedDB);
		$lightning_truth->{absolute} = $trueTaxonID_inUsedDB;

		if(exists $reads_inferred->{$readID})
		{
			my $inferredTaxonID = $reads_inferred->{$readID};
			
			my $lightning_inferred = $getLightning->($inferredTaxonID);
			my $attachedTo_inInference = $get_read_assignedToLevel->($readID);
			
			$lightning_inferred->{absolute} = $inferredTaxonID;
			
			foreach my $category (@read_categories)
			{
				$n_reads_correct{$category}{N}++; 
				$n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{N}++;

				$n_reads_correct{$category}{'attachedTo_' . $attachedTo_inInference}++;
				$n_reads_correct{$category}{'attachedToDirectlyMappable'} += ((exists $mappableTaxonIDs->{$inferredTaxonID}) ? 1 : 0);
					
				if($inferredTaxonID eq $trueTaxonID_inUsedDB)
				{
					$n_reads_correct{$category}{correct}++;
					$n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{correct}++;
				}
				else
				{
					if($category eq 'novel_to_genus')
					{
						# print join("\t", $category, $readID, $inferredTaxonID, $trueTaxonID_inUsedDB, $reads_truth_absolute->{$readID}), "\n";
					}
				}
			}
			
			my $printRead = (scalar(grep {$_ eq 'novel_to_superkingdom'} @read_categories) > 0);
			$printRead = 0;
			if($printRead)
			{
				print "Read $readID label $label\n";
				print "\ttrueTaxonID_inUsedDB: ", $trueTaxonID_inUsedDB, "\n";
				print "\tinferredTaxonID     : ", $inferredTaxonID, "\n";
				print "\tabsolute truth      : ", $reads_truth_absolute->{$readID}, "\n";
			}
			
			# genome-level unclassifiedness
			{
				my $mappable_absolute_truth = ($mappableTaxonIDs->{$trueTaxonID_inUsedDB}) ? 1 : 0;
				my $mappable_inferred =  ($mappableTaxonIDs->{$inferredTaxonID}) ? 1 : 0;

				foreach my $category (@read_categories)
				{					
					
					my $unclassified_key;
					if(($mappable_absolute_truth) and ($mappable_inferred))
					{
						$unclassified_key = 'N_unclassified_should0_is0';
					}
					elsif(($mappable_absolute_truth) and (not $mappable_inferred))
					{
						$unclassified_key = 'N_unclassified_should0_is1';
					}
					elsif((not $mappable_absolute_truth) and ($mappable_inferred))
					{
						$unclassified_key = 'N_unclassified_should1_is0';
					}
					elsif((not $mappable_absolute_truth) and (not $mappable_inferred))
					{
						$unclassified_key = 'N_unclassified_should1_is1'; 
					} 
										
					if($_read_categories{'truthLeafInDB'})
					{
						# open(MISSING, '>>missingTaxa') or die;
						# print MISSING $reads_truth_absolute->{$readID}, "\n";
						# close(MISSING);
						 
						# die Dumper(
							# "Weird - read $readID", "Not sure whether mappable or not?",
							# "trueTaxonID_inUsedDB $trueTaxonID_inUsedDB",
							# "reads_truth_mappingDB->{readID} <=> reads_truth_absolute->{readID}: $reads_truth_mappingDB->{$readID} //  $reads_truth_absolute->{$readID}"
						# ) unless($mappable_absolute_truth);
					}
					
					$n_reads_unknownStats_byLevel{$category}{'genome'}{$unclassified_key}++;
				}
			}
			
			foreach my $level ('absolute', @evaluateAccuracyAtLevels)
			{
				die unless((defined $lightning_truth->{$level}) and (defined defined $lightning_inferred->{$level}));
				if($printRead)
				{
					print "\t\t$level\n";
					print "\t\t\ttruth    : $lightning_truth->{$level}   \n";
					print "\t\t\tinference: $lightning_inferred->{$level}  \n";
				}
				

				foreach my $category (@read_categories)
				{					
					$n_reads_correct_byLevel{$category}{$level}{N}++;
					$n_reads_correct_byLevel{$category}{$level}{N_truthDefined}++ if($lightning_truth->{$level} ne 'NotLabelledAtLevel');
					if($level ne 'absolute')
					{
						$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{N}++;
					}
					
					$n_reads_correct_byLevel{$category}{$level}{'attachedTo_' . $attachedTo_inInference}++;
					$n_reads_correct_byLevel{$category}{$level}{'attachedToDirectlyMappable'} += ((exists $mappableTaxonIDs->{$inferredTaxonID}) ? 1 : 0);
					
					my $unclassified_key;
										
					if(($lightning_truth->{$level} ne 'Unclassified') and ($lightning_inferred->{$level} ne 'Unclassified'))
					{
						$unclassified_key = 'N_unclassified_should0_is0';
					}
					elsif(($lightning_truth->{$level} ne 'Unclassified') and ($lightning_inferred->{$level} eq 'Unclassified'))
					{
						$unclassified_key = 'N_unclassified_should0_is1';
					}
					elsif(($lightning_truth->{$level} eq 'Unclassified') and ($lightning_inferred->{$level} ne 'Unclassified'))
					{
						$unclassified_key = 'N_unclassified_should1_is0';
					}
					elsif(($lightning_truth->{$level} eq 'Unclassified') and ($lightning_inferred->{$level} eq 'Unclassified'))
					{
						$unclassified_key = 'N_unclassified_should1_is1';
					}
					die unless(defined $unclassified_key);

					$n_reads_unknownStats_byLevel{$category}{$level}{$unclassified_key}++;
					
					if($lightning_truth->{$level} eq $lightning_inferred->{$level})
					{
						# if(($level eq 'absolute') and ($category eq 'novel') and ($label =~ /Metamap-EM-Reads/))
						# {
							# die Dumper("How can this be?", $readID, \@read_categories, $level, $lightning_truth->{$level}, $lightning_inferred->{$level}, $lightning_truth, $lightning_inferred);
						# }
						
						$n_reads_correct_byLevel{$category}{$level}{correct}++;
						if($level ne 'absolute')
						{						
							$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{correct}++;						
						}
						if($inferredTaxonID eq $trueTaxonID_inUsedDB)
						{
							$n_reads_correct_byLevel{$category}{$level}{correct_exactly}++;
						}					

											
						$n_reads_correct_byLevel{$category}{$level}{correct_truthDefined}++ if($lightning_truth->{$level} ne 'NotLabelledAtLevel');
					}		
					else
					{
						if(($level eq 'family') and (scalar(grep {$_ eq 'truthInDB'} @read_categories )))
						{
							# print join("\t", "Debug read", $readID, $level, $lightning_truth->{$level}, $lightning_inferred->{$level}, $inferredTaxonID), "\n"; 
						}
					}
				}
			}
		}
		else
		{
			my $mappable_absolute_truth = ($mappableTaxonIDs->{$trueTaxonID_inUsedDB}) ? 1 : 0;
		
			my $unclassified_key_absolute;
			if($mappable_absolute_truth)
			{
				$unclassified_key_absolute = 'N_unclassified_should0_is0';
			}
			elsif(not $mappable_absolute_truth) 
			{
				$unclassified_key_absolute = 'N_unclassified_should1_is0';
			}

			foreach my $category (@read_categories)
			{
				$n_reads_correct{$category}{missing}++;			
				die unless($n_reads_correct{$category}{missing} <= scalar(@readIDs));
				$n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{missing}++;
				
				foreach my $level ('absolute', @evaluateAccuracyAtLevels)
				{
					$n_reads_correct_byLevel{$category}{$level}{missing}++;
					die unless($n_reads_correct_byLevel{$category}{$level}{missing} <= scalar(@readIDs));
					
					if($level ne 'absolute')
					{
						$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{missing}++;
					}
				}
				
				$n_reads_unknownStats_byLevel{$category}{'genome'}{$unclassified_key_absolute}++;		

				foreach my $level ('absolute', @evaluateAccuracyAtLevels)
				{
					die Dumper("No truth value for level $level", $lightning_truth, $trueTaxonID_inUsedDB) unless(defined $lightning_truth->{$level});
		
					my $unclassified_key;
											
					if($lightning_truth->{$level} ne 'Unclassified')
					{
						$unclassified_key = 'N_unclassified_should0_is0';
					}
					elsif($lightning_truth->{$level} eq 'Unclassified')
					{
						$unclassified_key = 'N_unclassified_should1_is0';
					}
		
					$n_reads_unknownStats_byLevel{$category}{$level}{$unclassified_key}++;		
				}
			}
		}
	}
	
	foreach my $category (keys %n_reads_correct)
	{
		my $prop_absoluteCorrectness = ($n_reads_correct{$category}{N} > 0) ? ($n_reads_correct{$category}{correct} / $n_reads_correct{$category}{N}) : -1;
		print join("\t", $label, $category, "absoluteCorrectness", $n_reads_correct{$category}{N}, $n_reads_correct{$category}{correct}, sprintf("%.2f", $prop_absoluteCorrectness)), "\n";
		foreach my $level (@evaluateAccuracyAtLevels)
		{
			my $prop_levelCorrectness = ($n_reads_correct_byLevel{$category}{$level}{N} > 0) ? ($n_reads_correct_byLevel{$category}{$level}{correct} / $n_reads_correct_byLevel{$category}{$level}{N}) : -1;
			print join("\t", $label, $category, $level, $n_reads_correct_byLevel{$category}{$level}{N}, $n_reads_correct_byLevel{$category}{$level}{correct}, sprintf("%.2f", $prop_levelCorrectness)), "\n";

		}
	}	

	if(defined $external_reads_correct)
	{
		foreach my $category (keys %n_reads_correct)
		{
			foreach my $key (keys %{$n_reads_correct{$category}})
			{
				$external_reads_correct->{$label}{$category}{$key} += $n_reads_correct{$category}{$key};
			}
		}
	}

	if(defined $external_reads_correct_byLevel)
	{
		foreach my $category (keys %n_reads_correct_byLevel)
		{
			foreach my $level (keys %{$n_reads_correct_byLevel{$category}})
			{
				foreach my $key (keys %{$n_reads_correct_byLevel{$category}{$level}})
				{
					$external_reads_correct_byLevel->{$label}{$category}{$level}{$key} += $n_reads_correct_byLevel{$category}{$level}{$key};
				}
				
				if(exists $external_reads_correct_byLevel->{$label}{$category}{$level}{missing})
				{
					#die unless($external_reads_correct_byLevel->{$label}{$category}{$level}{missing} <= $external_reads_correct_byLevel->{$label}{$category}{$level}{N});
				}
			}
		}
	}	
	
	if(defined $external_reads_correct_byLevel_byLength)
	{
		foreach my $category (keys %n_reads_correct_byLevel_byLength)
		{
			foreach my $level (keys %{$n_reads_correct_byLevel_byLength{$category}})
			{
				foreach my $key (keys %{$n_reads_correct_byLevel_byLength{$category}{$level}})
				{
					foreach my $rL (keys %{$n_reads_correct_byLevel_byLength{$category}{$level}{$key}})
					{
						$external_reads_correct_byLevel_byLength->{$label}{$category}{$level}{$key}{$rL} += $n_reads_correct_byLevel_byLength{$category}{$level}{$key}{$rL};
					}
				}
			}
		}	
	}
	
	if(defined $external_n_reads_unknownStats_byLevel_href)
	{
		foreach my $category (keys %n_reads_unknownStats_byLevel)
		{
			foreach my $level (keys %{$n_reads_unknownStats_byLevel{$category}})
			{
				foreach my $key (keys %{$n_reads_unknownStats_byLevel{$category}{$level}})
				{
					$external_n_reads_unknownStats_byLevel_href->{$label}{$category}{$level}{$key} += $n_reads_unknownStats_byLevel{$category}{$level}{$key};
				}
			}
		}	
	}
}

sub getEmptyGlobalResltsStore
{
	my $allSimulations_data_href = {};
	$allSimulations_data_href->{n_reads_correct_byVariety_bySimulation} = [];
	$allSimulations_data_href->{n_reads_correct_byVariety_byLevel_bySimulation} = [];
	$allSimulations_data_href->{n_reads_correct_byVariety_byLevel_byLength_bySimulation} = [];	
	$allSimulations_data_href->{freq_byVariety_byLevel_bySimulation} = [];
	$allSimulations_data_href->{directlyMappable_bySimulation} = [];
	$allSimulations_data_href->{n_reads_correct_byVariety} = {};
	$allSimulations_data_href->{n_reads_correct_byVariety_byLevel} = {};
	$allSimulations_data_href->{freq_byVariety_byLevel} = {};
	$allSimulations_data_href->{frequencyComparisons_bySimulation} = [];
	$allSimulations_data_href->{frequencyComparisons_details_bySimulation} = [];
	$allSimulations_data_href->{highLevel_stats_keptSeparate_bySimulation} = [];
	$allSimulations_data_href->{callRate_and_accuracy_byReadCategory} = {};
	$allSimulations_data_href->{callRate_and_accuracy_byReadCategory_byLength} = {};
	$allSimulations_data_href->{attachedTo_byReadCategory} = {};
	$allSimulations_data_href->{realizedN} = 0;
	return $allSimulations_data_href;
}

sub analyseAndAddOneExperiment
{
	my $extendedMaster = shift;
	my $reduced_taxonID_master_2_contigs_href = shift;
	
	my $truth_reads_href = shift;
	my $readLengths_href = shift;
	
	my $methodNames_aref = shift;	
	my $inferred_reads_aref = shift;
	my $inferred_distributions_aref = shift;
	
	my $varietyName_forStorage = shift;
	
	die unless(defined $extendedMaster);
	die unless(defined $varietyName_forStorage);
	
	my $allSimulations_data_href = getEmptyGlobalResltsStore();

	my $frequencyComparison = {};
	my %n_reads_correct_byVariety_local;
	my %n_reads_correct_byVariety_byLevel_local;	
	my %n_reads_correct_byVariety_byLevel_byLength_local;		
	my %n_reads_unknownStats_byLevel;		
	my %freq_byVariety_byLevel_local;		
	#my %freq_details_byVariety_byLevel_local;	
	my %directlyMappable;	

	foreach my $tID (keys %$reduced_taxonID_master_2_contigs_href)
	{
		$directlyMappable{$varietyName_forStorage}{$tID} = 1;
	}
		
	$n_reads_correct_byVariety_local{$varietyName_forStorage} = {};
	$n_reads_correct_byVariety_byLevel_local{$varietyName_forStorage} = {};
	$n_reads_correct_byVariety_byLevel_byLength_local{$varietyName_forStorage} = {};
	$n_reads_unknownStats_byLevel{$varietyName_forStorage} = {};
	$freq_byVariety_byLevel_local{$varietyName_forStorage} = {};
	$frequencyComparison->{$varietyName_forStorage} = {};
	
	my $mappableTaxonomy = dclone $extendedMaster;
	taxTree::removeUnmappableParts($mappableTaxonomy, $reduced_taxonID_master_2_contigs_href);
		
	my $truth_reads_mappable = validation::translateReadsTruthToReducedTaxonomy($extendedMaster, $mappableTaxonomy, $truth_reads_href);
	
	my $truth_reads_href_noUnknown = { map {$_ => $truth_reads_mappable->{$_}} grep {$truth_reads_mappable->{$_}} keys %$truth_reads_mappable };
	die Dumper("Did you input undefined reads?", scalar(keys %{$truth_reads_href_noUnknown}), scalar(keys %{$truth_reads_href})) unless(scalar(keys %{$truth_reads_href_noUnknown}) == scalar(keys %{$truth_reads_href}));
	
	die unless(all {exists $mappableTaxonomy->{$_}} values %$truth_reads_href_noUnknown);
	my $truth_mappingDatabase_distribution = validation::truthReadsToTruthSummary($mappableTaxonomy, $truth_reads_href_noUnknown, $reduced_taxonID_master_2_contigs_href);

	die unless($#{$inferred_reads_aref} == $#{$methodNames_aref});
	die unless($#{$inferred_reads_aref} == $#{$inferred_distributions_aref});
	
	for(my $methodI = 0; $methodI <= $#{$inferred_reads_aref}; $methodI++)
	{
		if(defined $inferred_reads_aref->[$methodI])
		{
			validation::readLevelComparison(
				$extendedMaster,
				$truth_reads_href,
				$truth_reads_mappable,
				$inferred_reads_aref->[$methodI],
				$methodNames_aref->[$methodI],
				$n_reads_correct_byVariety_local{$varietyName_forStorage},
				$n_reads_correct_byVariety_byLevel_local{$varietyName_forStorage},
				$n_reads_correct_byVariety_byLevel_byLength_local{$varietyName_forStorage},
				$reduced_taxonID_master_2_contigs_href,
				$readLengths_href,
				$n_reads_unknownStats_byLevel{$varietyName_forStorage}
			);
		}
		
		if(defined $inferred_distributions_aref->[$methodI])
		{
			validation::distributionLevelComparison(
				$extendedMaster,
				$truth_mappingDatabase_distribution,
				$inferred_distributions_aref->[$methodI],
				$methodNames_aref->[$methodI],
				$freq_byVariety_byLevel_local{$varietyName_forStorage},
				$frequencyComparison->{$varietyName_forStorage}
			);		
		}			
	}

		
	my %taxonIDs_in_direct_truth = map {$_ => 1} values %$truth_reads_href;
	foreach my $taxonID (keys %taxonIDs_in_direct_truth)
	{
		next if($taxonID eq '0');
		my @descendants = taxTree::descendants($extendedMaster, $taxonID);
		foreach my $descendantID (@descendants)
		{
			die unless(exists $extendedMaster->{$descendantID});
			die if(exists $taxonIDs_in_direct_truth{$descendantID});
		}
	}

	my $truth_reads_novelTree = validation::truthReadsTree($extendedMaster, $truth_reads_mappable, \%taxonIDs_in_direct_truth);
	my $unknown_and_frequencyContributions_href = {};
	$unknown_and_frequencyContributions_href->{$varietyName_forStorage}{mappable} = {};
	foreach my $taxonID (keys %$reduced_taxonID_master_2_contigs_href)
	{
		$unknown_and_frequencyContributions_href->{$varietyName_forStorage}{mappable}{$taxonID} = 1;
	}
	$unknown_and_frequencyContributions_href->{$varietyName_forStorage}{truthReadsNovelTree} = $truth_reads_novelTree;
	$unknown_and_frequencyContributions_href->{$varietyName_forStorage}{nReads} = scalar(keys %$truth_reads_href);		
	my $frequencyComparison_details = $unknown_and_frequencyContributions_href;
		
	validation::addResultsToGlobalStore(
		0,
		$allSimulations_data_href, 
		$frequencyComparison,
		$frequencyComparison_details,
		\%n_reads_correct_byVariety_local,
		\%n_reads_correct_byVariety_byLevel_local,
		\%n_reads_correct_byVariety_byLevel_byLength_local,
		\%n_reads_unknownStats_byLevel,
		\%freq_byVariety_byLevel_local,
		undef, # apparently not necessary: \%freq_details_byVariety_byLevel_local,
		\%directlyMappable,
	);	

	return $allSimulations_data_href;

}

sub addResultsToGlobalStore
{
	my $jobI = shift;
	my $allSimulations_data_href = shift;
	my $frequencyComparison = shift;
	my $frequencyComparison_details = shift;
	my $n_reads_correct_byVariety_local_href = shift;
	my $n_reads_correct_byVariety_byLevel_local_href = shift;
	my $n_reads_correct_byVariety_byLevel_byLength_local_href = shift;
	my $n_reads_unknownStats_byLevel_href = shift;
	my $freq_byVariety_byLevel_local_href = shift;
	#my $freq_details_byVariety_byLevel_local_href = shift;
	shift;
	my $directlyMappable_href = shift;
		
		
	# variety = fullDB/removeOne_genus ...
	# label = MetaMap / Kraken ...
	# category = read category ...

	
	# my %n_reads_correct_byVariety_byLevel_byLength;
	foreach my $variety (keys %$n_reads_unknownStats_byLevel_href)
	{		
		my $variety_forStore = $variety;
		$variety_forStore =~ s/_\d+$//;		
		foreach my $label (keys %{$n_reads_unknownStats_byLevel_href->{$variety}}) 
		{
			foreach my $category (keys %{$n_reads_unknownStats_byLevel_href->{$variety}{$label}})
			{
				foreach my $level (keys %{$n_reads_unknownStats_byLevel_href->{$variety}{$label}{$category}})
				{
					foreach my $key (keys %{$n_reads_unknownStats_byLevel_href->{$variety}{$label}{$category}{$level}})
					{
						my $value = $n_reads_unknownStats_byLevel_href->{$variety}{$label}{$category}{$level}{$key};
						die unless(not ref($value));
						push(@{$allSimulations_data_href->{n_reads_unknownStats_byLevel}->{$variety_forStore}{$label}{$category}{$level}{$key}}, $value);
					}						 
				}
			}
		}
	}
		
	foreach my $variety (keys %$n_reads_correct_byVariety_local_href)
	{		
		foreach my $label (keys %{$n_reads_correct_byVariety_local_href->{$variety}}) 
		{
			foreach my $category (keys %{$n_reads_correct_byVariety_local_href->{$variety}{$label}})
			{
				# die Dumper([keys %{$n_reads_correct_byVariety_local_href->{$variety}{$label}{$category}}]);
					# 'attachedTo_species',
					# 'N',
					# 'missing',
					# 'correct',
					# 'attachedToDirectlyMappable'
				
				#die Dumper($jobI, $variety, $label, $category);
				my $variety_forStore = $variety;
				$variety_forStore =~ s/_\d+$//;
				
				if($variety_forStore ne $variety)
				{
					# die Dumper($variety, $variety_forStore);
				}
				
				foreach my $key (keys %{$n_reads_correct_byVariety_local_href->{$variety}{$label}{$category}})
				{
					my $value = $n_reads_correct_byVariety_local_href->{$variety}{$label}{$category}{$key};
					die unless(not ref($value));
					push(@{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety_forStore}{$label}{$category}{$key}}, $value); # hopefully fixed
				}
				
				my $d = $n_reads_correct_byVariety_local_href->{$variety}{$label}{$category};
				die unless(exists $d->{N});
				die unless(exists $d->{missing});
				die unless(exists $d->{correct});
				my $N = $d->{N} + $d->{missing};
				die unless($N > 0);
				die unless($d->{N} > 0);
			
				my $CR = $d->{N} / $N; die unless(($CR >= 0) and ($CR <= 1));
				my $accuracy = $d->{correct} / $d->{N}; die unless(($accuracy >= 0) and ($accuracy <= 1));
				
				push(@{$allSimulations_data_href->{highLevel_stats_keptSeparate_bySimulation}->[$jobI]{$variety_forStore}{$label}{$category}{absolute}{CR}}, $CR); # hopefully ok
				push(@{$allSimulations_data_href->{highLevel_stats_keptSeparate_bySimulation}->[$jobI]{$variety_forStore}{$label}{$category}{absolute}{Accuracy}}, $accuracy); # hopefully ok
				
				if(($variety !~ /allCombined/) and ($variety !~ /incompleteCombined/))
				{
					
					push(@{$allSimulations_data_href->{callRate_and_accuracy_byReadCategory}->{$category}{$label}{absolute}}, [$CR, $accuracy, $accuracy]); # hopefully ok
					
					my @keys_attachedTo = grep {$_ =~ /^attachedTo/} keys %$d;
					die unless(scalar(@keys_attachedTo));
					
					my %localAttachedHash;  
					foreach my $key (@keys_attachedTo)
					{
						$localAttachedHash{$key} = $d->{$key} / $d->{N}; # hopefully ok
					}	

					push(@{$allSimulations_data_href->{attachedTo_byReadCategory}->{$category}{$label}{absolute}}, \%localAttachedHash); # hopefully ok

					# t odo
					# warn Dumper("Add $variety to attachedTo $category $label $variety", $allSimulations_data_href->{attachedTo_byReadCategory}->{$category}{$label}{absolute});

				}
			}
		}
	}
	
	# die Dumper($allSimulations_data_href->{attachedTo_byReadCategory}
	# exit;

	foreach my $variety (keys %$n_reads_correct_byVariety_byLevel_local_href)
	{		
		my $variety_forStore = $variety;
		$variety_forStore =~ s/_\d+$//;							
				
		foreach my $label (keys %{$n_reads_correct_byVariety_byLevel_local_href->{$variety}})
		{		
			foreach my $category (keys %{$n_reads_correct_byVariety_byLevel_local_href->{$variety}{$label}})
			{
				foreach my $level (keys %{$n_reads_correct_byVariety_byLevel_local_href->{$variety}{$label}{$category}})
				{
					foreach my $key (keys %{$n_reads_correct_byVariety_byLevel_local_href->{$variety}{$label}{$category}{$level}})
					{
						my $value = $n_reads_correct_byVariety_byLevel_local_href->{$variety}{$label}{$category}{$level}{$key};
						die unless(not ref($value));
						push(@{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety_forStore}{$label}{$category}{$level}{$key}}, $value); # hopefully fixed
					}


					my $d = $n_reads_correct_byVariety_byLevel_local_href->{$variety}{$label}{$category}{$level};
					die unless(exists $d->{N});
					die unless(exists $d->{missing});
					die unless(exists $d->{correct});
					die unless(exists $d->{correct_exactlyAtLevel}); 
					my $N = $d->{N} + $d->{missing};
					die unless($N > 0);
					die unless($d->{N} > 0);
					 
					my $CR = $d->{N} / $N; die unless(($CR >= 0) and ($CR <= 1));
					my $accuracy = $d->{correct} / $d->{N}; die unless(($accuracy >= 0) and ($accuracy <= 1));
					my $accuracy_exactlyAtLevel = $d->{correct_exactlyAtLevel} / $d->{N}; die unless(($accuracy_exactlyAtLevel >= 0) and ($accuracy_exactlyAtLevel <= 1));
					
					push(@{$allSimulations_data_href->{highLevel_stats_keptSeparate_bySimulation}->[$jobI]{$variety_forStore}{$label}{$category}{$level}{CR}}, $CR); # hopefully ok
					push(@{$allSimulations_data_href->{highLevel_stats_keptSeparate_bySimulation}->[$jobI]{$variety_forStore}{$label}{$category}{$level}{Accuracy}}, $accuracy);  # hopefully ok
					push(@{$allSimulations_data_href->{highLevel_stats_keptSeparate_bySimulation}->[$jobI]{$variety_forStore}{$label}{$category}{$level}{AccuracyExactlyAtLevel}}, $accuracy_exactlyAtLevel); # hopefully ok	
					 
					if(($variety !~ /allCombined/) and ($variety !~ /incompleteCombined/) and ($level ne 'absolute'))
					{
						push(@{$allSimulations_data_href->{callRate_and_accuracy_byReadCategory}->{$category}{$label}{$level}}, [$CR, $accuracy, $accuracy_exactlyAtLevel]);	 # hopefully ok			
						
						my @keys_attachedTo = grep {$_ =~ /^attachedTo/} keys %$d;
						die unless(scalar(@keys_attachedTo));
						
						my %localAttachedHash;
						foreach my $key (@keys_attachedTo)
						{
							$localAttachedHash{$key} = $d->{$key} / $d->{N};
						}	
						#warn "Attach I $variety $variety_forStore $category $label $level";
						#warn "\t", scalar(@{$allSimulations_data_href->{attachedTo_byReadCategory}->{$category}{$label}{$level}}), "\n";						
						push(@{$allSimulations_data_href->{attachedTo_byReadCategory}->{$category}{$label}{$level}}, \%localAttachedHash); # hopefully ok
						#warn "\t", scalar(@{$allSimulations_data_href->{attachedTo_byReadCategory}->{$category}{$label}{$level}}), "\n";
					}
				}
				
				die unless(defined $n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category});
				foreach my $level (keys %{$n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}})
				{
					foreach my $rL (keys %{$n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}{$level}})
					{
						foreach my $key (keys %{$n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}{$level}{$rL}})
						{
							my $value = $n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}{$level}{$rL}{$key};
							die Dumper($level, $rL, $key, $n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}) unless(defined $value);
							die unless(not ref($value));
						}
						
						my $d = $n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}{$level}{$rL};
						die Dumper('N', $level, $rL, $d) unless(exists $d->{N});
						die Dumper('missing', $level, $rL, $d) unless(exists $d->{missing});
						die unless(exists $d->{correct});				
						
						my $N = $d->{N} + $d->{missing}; die unless($d > 0);
						my $CR = $d->{N} / $N; die unless(($CR >= 0) and ($CR <= 1));
						my $accuracy = ($d->{N} > 0) ? ($d->{correct} / $d->{N}) : -1; die unless(($accuracy >= -1) and ($accuracy <= 1));

						if(($variety !~ /allCombined/) and ($variety !~ /incompleteCombined/))
						{	
							push(@{$allSimulations_data_href->{callRate_and_accuracy_byReadCategory_byLength}->{$category}{$label}{$level}{$rL}}, [$CR, $accuracy]); # seems OK						
						}
					}				
				}
			}
		}
	}
	
	# the following is no x/y data, but summary stats of accuracy
	foreach my $variety (keys %$freq_byVariety_byLevel_local_href)
	{			
		foreach my $label (keys %{$freq_byVariety_byLevel_local_href->{$variety}})
		{
			foreach my $level (keys %{$freq_byVariety_byLevel_local_href->{$variety}{$label}})
			{
			
				my $variety_forStore = $variety;
				$variety_forStore =~ s/_\d+$//;
								
				#die Dumper("freqOK", $freq_byVariety_byLevel_local_href->{$variety}{$label}) unless(exists $freq_byVariety_byLevel_local_href->{$variety}{$label}{freqOK});
				#die Dumper("L1", $freq_byVariety_byLevel_local_href->{$variety}{$label}) unless(exists $freq_byVariety_byLevel_local_href->{$variety}{$label}{L1});
				foreach my $key (keys %{$freq_byVariety_byLevel_local_href->{$variety}{$label}{$level}})
				{
					my $value = $freq_byVariety_byLevel_local_href->{$variety}{$label}{$level}{$key};
					die unless(defined $value);
					die unless((not ref($value)) or (ref($value) eq 'ARRAY'));
					if(not ref($value))
					{
						die "This should not happen!";
						$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety_forStore}{$label}{$level}{$key} += $value; # this should be OK
					}
					else
					{ 
						push(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety_forStore}{$label}{$level}{$key}}, @$value); # hopefully fixed
					}
				}
			}
		}
	}				
		
	push(@{$allSimulations_data_href->{n_reads_correct_byVariety_bySimulation}}, $n_reads_correct_byVariety_local_href);
	push(@{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel_bySimulation}}, $n_reads_correct_byVariety_byLevel_local_href);
	push(@{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel_byLength_bySimulation}}, $n_reads_correct_byVariety_byLevel_byLength_local_href);
	push(@{$allSimulations_data_href->{freq_byVariety_byLevel_bySimulation}}, $freq_byVariety_byLevel_local_href);			
	
	push(@{$allSimulations_data_href->{frequencyComparisons_bySimulation}}, $frequencyComparison);
	push(@{$allSimulations_data_href->{frequencyComparisons_details_bySimulation}}, $frequencyComparison_details);
	push(@{$allSimulations_data_href->{directlyMappable_bySimulation}}, $directlyMappable_href);	
}

sub distributionLevelComparison
{	
	my $masterTaxonomy = shift;
	my $distribution_truth = shift;
	my $distribution_inferred = shift;
	my $label = shift;
	my $external_comparison = shift;
	my $frequencyComparison_href = shift;
	# my $distributionLevelComparison = shift;

	foreach my $level ('absolute', 'definedGenomes', 'definedAndHypotheticalGenomes', @evaluateAccuracyAtLevels)
	{
		my $lookupKey_level_InInference = $level;
		if(($level eq 'definedAndHypotheticalGenomes') and (not defined $distribution_inferred->{$lookupKey_level_InInference}))
		{
			$lookupKey_level_InInference = 'definedGenomes';
		}
		
		if($level eq 'absolute')
		{
			if(defined $distribution_inferred->{'definedAndHypotheticalGenomes'})
			{	
				$lookupKey_level_InInference = 'definedAndHypotheticalGenomes';
				# die unless(defined $distribution_inferred->{'definedGenomes'});
			}
			elsif(defined $distribution_inferred->{'definedGenomes'})
			{
				$lookupKey_level_InInference = 'definedGenomes';
			}
			else
			{
				die if(defined $distribution_inferred->{$level});
			}
		}
		next unless(defined $distribution_inferred->{$lookupKey_level_InInference});
		die unless(defined $distribution_truth->{$level});
		
		my $totalFreq = 0;
		my $totalFreqCorrect = 0;
		foreach my $inferredTaxonID (keys %{$distribution_inferred->{$lookupKey_level_InInference}})
		{
			my $isFreq = $distribution_inferred->{$lookupKey_level_InInference}{$inferredTaxonID}[1];
			my $shouldBeFreq = 0;
			if(exists $distribution_truth->{$level}{$inferredTaxonID})
			{
				$shouldBeFreq = $distribution_truth->{$level}{$inferredTaxonID};
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
			
			if($level eq 'genus')
			{
				# print join("\t", "Debug output level $level", $inferredTaxonID, $isFreq, $shouldBeFreq), "\n";
			}	
		}

		die Dumper("Weird total freq", $label, $totalFreq, $level) unless(abs(1 - $totalFreq) <= 1e-3);
		
		my $S_AVGRE = 0;
		my $S_RRMSE = 0;
		foreach my $trueTaxonID (keys %{$distribution_truth->{$level}})
		{
			my $shouldBeFreq = $distribution_truth->{$level}{$trueTaxonID};
			die Dumper("Problem with shouldBeFreq", $shouldBeFreq, $distribution_truth) unless($shouldBeFreq > 0);
			my $isFreq = (exists $distribution_inferred->{$lookupKey_level_InInference}{$trueTaxonID}) ? $distribution_inferred->{$lookupKey_level_InInference}{$trueTaxonID}[1] : 0;
			$S_AVGRE += ( abs($shouldBeFreq - $isFreq) / $shouldBeFreq);
			$S_RRMSE += (($shouldBeFreq - $isFreq) / $shouldBeFreq)**2;
		}
		
		my $L1_sum = 0;
		my $L2_sum = 0;
		my @is_bigger0;
		my @should_bigger0;
		
				
		# my $unclassified_is = 0;
		# my $unclassified_shouldBe = 0;
		
		my %joint_taxonIDs = map {$_ => 1} ((keys %{$distribution_inferred->{$lookupKey_level_InInference}}), (keys %{$distribution_truth->{$level}}));
		foreach my $taxonID (keys %joint_taxonIDs)
		{
			my $isFreq = (exists $distribution_inferred->{$lookupKey_level_InInference}{$taxonID}) ? $distribution_inferred->{$lookupKey_level_InInference}{$taxonID}[1] : 0;
			my $shouldBeFreq = (exists $distribution_truth->{$level}{$taxonID}) ? $distribution_truth->{$level}{$taxonID} : 0;
			my $L1_diff = abs($isFreq - $shouldBeFreq);
			my $L2_diff = ($isFreq - $shouldBeFreq)**2;
			$L1_sum += $L1_diff;
			$L2_sum += $L2_diff;
			
			if($taxonID eq '1280')
			{
				# print join("\t", "VALUE IN COMPARSION for $label // $level: ", $taxonID, $shouldBeFreq, $isFreq), "\n";
			}
			
			die if(defined $frequencyComparison_href->{$label}{$level}{$taxonID});
			$frequencyComparison_href->{$label}{$level}{$taxonID} = [$shouldBeFreq, $isFreq];
			
			if(($shouldBeFreq > 0) or ($isFreq > 0))
			{
				push(@should_bigger0, $shouldBeFreq);			
				push(@is_bigger0, $isFreq);
			}
			
			# if($taxonID eq 'Unclassified')
			# {
				# $unclassified_shouldBe = $shouldBeFreq;
				# $unclassified_is = $isFreq;
			# }
		}
		
		die unless(scalar(@should_bigger0) == scalar(@is_bigger0));
		my $L1 = $L1_sum;
		my $L2 = sqrt($L2_sum);
		my $r = 0+correlation( \@should_bigger0, \@is_bigger0 );
		die Dumper("Weird r", $r, \@should_bigger0, \@is_bigger0) unless(($r >= -1*(1+1e-5)) and ($r <= 1+1e-5));
		my $r2 = $r ** 2;
		
		my $AVGRE *= (1 / scalar(keys %{$distribution_inferred->{$lookupKey_level_InInference}}));
		my $RRMSE *= (1 / scalar(keys %{$distribution_inferred->{$lookupKey_level_InInference}}));
		$RRMSE = sqrt($RRMSE);
		
		print join("\t", $label, $level, $totalFreqCorrect), "\n";
		
		
		if(defined $external_comparison)
		{
			# $external_comparison->{$label}{$level}{total} += $totalFreq; 
			# $external_comparison->{$label}{$level}{correct} += $totalFreqCorrect; 
			my $freqOK = $totalFreqCorrect / $totalFreq;
			die unless(($freqOK >= 0) and ($freqOK <= 1));
			
			push(@{$external_comparison->{$label}{$level}{freqOK}}, $freqOK);
			push(@{$external_comparison->{$label}{$level}{AVGRE}}, $AVGRE);
			push(@{$external_comparison->{$label}{$level}{RRMSE}}, $RRMSE); 
			push(@{$external_comparison->{$label}{$level}{L1}}, $L1);
			push(@{$external_comparison->{$label}{$level}{L2}}, $L2);
			push(@{$external_comparison->{$label}{$level}{r2}}, $r2);
		}
	}		
}			
			

sub readInferredDistribution
{
	my $taxonomy = shift;
	my $taxonomy_merged = shift;
	my $f = shift;
	die unless(defined $f);
	
	my %inference;
	open(I, '<', $f) or die "Cannot open $f";
	my $header_line = <I>;
	chomp($header_line);
	my @header_fields = split(/\t/, $header_line);
	my %freq_per_level;
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
		next if($line{Name} eq 'totalReads');
		next if($line{Name} eq 'readsLongEnough');
		next if($line{Name} eq 'readsLongEnough_unmapped');
		
		if(((substr($taxonID_nonMaster, 0, 1) ne 'x') and ($taxonID_nonMaster <= 0)) and (($line{Name} eq 'Undefined') or ($line{Name} eq 'Unclassified') or ($line{Name} eq 'NotLabelledAtLevel')or ($line{Name} eq 'NotLabelledAtLevel')))
		{
			$taxonID_nonMaster = $line{Name};
		}
		if($taxonID_nonMaster eq '0')
		{
			Dumper($line, $f);
		}
		
		my $taxonID_master = taxTree::findCurrentNodeID($taxonomy, $taxonomy_merged, $taxonID_nonMaster);
		if(($taxonID_master eq 'Undefined') or ($taxonID_master eq 'NotLabelledAtLevel'))
		{
			$taxonID_master = 'Unclassified';
		}
		die Dumper(\%line) unless(defined $taxonID_master);
		die Dumper("Unknown taxon ID $taxonID_master in file $f", $taxonomy->{$taxonID_master}, $taxonID_nonMaster, \%line) unless(($taxonID_master eq 'Unclassified') or (defined $taxonomy->{$taxonID_master}));
		die Dumper("Weird line in $f", \%line) unless(defined $line{Absolute});
		die unless(defined $line{PotFrequency});
		$inference{$line{AnalysisLevel}}{$taxonID_master}[0] += $line{Absolute};
		die unless(exists $line{PotFrequency});
		$inference{$line{AnalysisLevel}}{$taxonID_master}[1] += $line{PotFrequency};
		$freq_per_level{$line{AnalysisLevel}} += $line{PotFrequency};  
		
		# if(exists $line{EMFrequency})
		# {
			# $inference{$line{AnalysisLevel}}{$taxonID_master}[1] += $line{EMFrequency};
		# }
		# else
		# {
			# $inference{$line{AnalysisLevel}}{$taxonID_master}[1] += $line{PotFrequency};
		# }
	}
	close(I);	
	
	foreach my $level (keys %freq_per_level)
	{
		my $S = $freq_per_level{$level};
		unless(abs($S - 1) <= 1e-3)
		{
			die "Weird sum of frequencies in file $f for level $level";
		}
	}	
	
	return \%inference;
}

sub readInferredFileReads
{
	my $masterTaxonomy = shift;
	my $masterTaxonomy_merged = shift;
	my $file = shift;
	
	die unless(defined $file);
	my %inferred_raw_reads;
	
	open(I, '<', $file) or die "Cannot open $file";
	while(<I>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		die "Weird line $line in $file $." unless(scalar(@f) == 2);
		my $readID = $f[0];
		my $taxonID_original = $f[1];
		
		# get the taxon ID in the master taxonomy
		my $taxonID_master = taxTree::findCurrentNodeID($masterTaxonomy, $masterTaxonomy_merged, $taxonID_original);
		
		if($taxonID_master eq 'Unclassified')
		{
			die;
			$taxonID_master = 1;
		}
		
		die "Undefined taxon ID $taxonID_master in file $file $." unless(($taxonID_master eq '0') or (exists $masterTaxonomy->{$taxonID_master}));
		die if(defined $inferred_raw_reads{$readID});
		
		$inferred_raw_reads{$readID} = $taxonID_master;
	}
	close(I);	
	
	return \%inferred_raw_reads;	
}

 
sub produceValidationOutputFiles
{
	my $allSimulations_data_href = shift;
	my $fullTaxonomy_simulation = shift;
	my $prefix_for_outputFiles = shift;
	my $outputDir_allSimulations = shift;
	my $suffix = shift;
	
	my $assumeHaveHeader_READSCORRECTBYLEVEL_ALL = ((-e $outputDir_allSimulations . '/_readsCorrectByLevel') && (-s $outputDir_allSimulations . '/_readsCorrectByLevel'));
	my $assumeHaveHeader_FREQEVALUATION_ALL = ((-e $outputDir_allSimulations . '/_frequenciesCorrectByLevel') && (-s $outputDir_allSimulations . '/_frequenciesCorrectByLevel'));
	
	open(READSCORRECTBYLEVEL_ALL, '>>', $outputDir_allSimulations . '/_all_readsCorrectByLevel') or die;
	open(FREQEVALUATION_ALL, '>>', $outputDir_allSimulations . '/_all_frequenciesCorrectByLevel') or die;
	open(UNCLASSIFIED_ALL, '>>', $outputDir_allSimulations . '/_all_unclassifiedSummary_reads') or die;
	open(UNCLASSIFIED_FREQ_ALL, '>>', $outputDir_allSimulations . '/_all_unclassifiedSummary_frequencies') or die;
	

	my @varieties = qw/allCombined fullDB incompleteCombined removeOne_self removeOne_species removeOne_genus/;
	@varieties = grep {exists $allSimulations_data_href->{n_reads_correct_byVariety}->{$_}} @varieties;
	#die Dumper(\@varieties, \%n_reads_correct_byVariety) unless(scalar(@varieties) == scalar(keys %n_reads_correct_byVariety));
	die Dumper(\@varieties, [keys %{$allSimulations_data_href->{n_reads_correct_byVariety}}], "Issue I") unless(all {exists $allSimulations_data_href->{n_reads_correct_byVariety}->{$_}} @varieties);

	my @levels_ordered = validation::getEvaluationLevels();
	my %level_to_i;
	for(my $levelI = 0; $levelI <= $#levels_ordered; $levelI++)
	{
		$level_to_i{$levels_ordered[$levelI]} = $levelI;
	}
	
	{					
		my %_methods;
		my %_readStratification;
		my %_evaluationLevels;
		foreach my $variety (@varieties)
		{
			foreach my $methodName (keys %{$allSimulations_data_href->{n_reads_unknownStats_byLevel}->{$variety}})
			{
				$_methods{$methodName}++;
				foreach my $category (keys %{$allSimulations_data_href->{n_reads_unknownStats_byLevel}->{$variety}{$methodName}})
				{
					$_readStratification{$category}++;
					
					if(exists $allSimulations_data_href->{n_reads_unknownStats_byLevel}->{$variety}{$methodName}{$category})
					{
						foreach my $k (keys %{$allSimulations_data_href->{n_reads_unknownStats_byLevel}->{$variety}{$methodName}{$category}})
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
		}} keys %_readStratification;
		
		# my @evaluationLevels = sort keys %_evaluationLevels;
		my @evaluationLevels = qw/genome species genus family/;
		die Dumper("Missing evaluation levels I", \@evaluationLevels, \%_evaluationLevels, [\@varieties]) unless(all {exists $_evaluationLevels{$_}} @evaluationLevels);
		
			
		my $fn_unclassified =  $prefix_for_outputFiles . '_unclassifiedSummary_reads';
		open(UNCLASSIFIED, '>', $fn_unclassified) or die;		
		
		my @header_fields_1_byLevelUnknown = ('ReadLevel', 'EvaluationLevel');
		my @header_fields_2_byLevelUnknown = ('', '');
		my @header_fields_3_byLevelUnknown = ('', '');	
		
		foreach my $variety (@varieties)
		{
			my $hf2_before = $#header_fields_2_byLevelUnknown;
			foreach my $method (@methods)
			{
				push(@header_fields_2_byLevelUnknown, $method, '', '', '', '', '', '');							
				push(@header_fields_3_byLevelUnknown, 'averagedOver', 'readsAvg', 'shouldBeUnclassifiedAvg', 'isUnclassifiedAvg', 'sensitivityAvg', 'PPVAvg', 'specificityAvg');
			}
			my $hf2_after = $#header_fields_2_byLevelUnknown;
			my $requiredFields = $hf2_after - $hf2_before;
			die unless($requiredFields > 0);
			my @addToHeader1 = ($variety, (('') x ($requiredFields - 1)));
			die unless(scalar(@addToHeader1) == $requiredFields);
			push(@header_fields_1_byLevelUnknown, @addToHeader1);
		}
		
		print UNCLASSIFIED join("\t", @header_fields_1_byLevelUnknown), "\n";
		print UNCLASSIFIED join("\t", @header_fields_2_byLevelUnknown), "\n";
		print UNCLASSIFIED join("\t", @header_fields_3_byLevelUnknown), "\n";
		
		# unless($assumeHaveHeader_READSCORRECTBYLEVEL_ALL)
		{
			print UNCLASSIFIED_ALL join("\t", 'Experiment', @header_fields_1_byLevelUnknown), "\n";
			print UNCLASSIFIED_ALL join("\t", '', @header_fields_2_byLevelUnknown), "\n";
			print UNCLASSIFIED_ALL join("\t", '', @header_fields_3_byLevelUnknown), "\n";
		}
							
		foreach my $readLevel (@readLevels)
		{
			foreach my $evaluationLevel (@evaluationLevels)
			{
				my @output_fields_byLevelUnknown = ($readLevel, $evaluationLevel);
				foreach my $variety (@varieties)
				{
					foreach my $methodName (@methods)
					{		
						my $v = $allSimulations_data_href->{n_reads_unknownStats_byLevel}->{$variety}{$methodName}{$readLevel}{$evaluationLevel};
						
						my @Ns;
						my @shouldBeUnclassified;
						my @isUnclassified;
						my @sensitivities;
						my @PPVs;
						my @specificities;

						my $n_experiments = $#{$v->{N_unclassified_should1_is1}}+1;
						
						for(my $i = 0; $i < $n_experiments; $i++)
						{
							die unless(defined $v->{N_unclassified_should0_is0}[$i]);
							die unless(defined $v->{N_unclassified_should0_is1}[$i]);
							die unless(defined $v->{N_unclassified_should1_is0}[$i]);
							die unless(defined $v->{N_unclassified_should1_is1}[$i]);
							die unless(scalar(keys %{$v}) == 4);
							my $shouldBeUnclassified = $v->{N_unclassified_should1_is0}[$i] + $v->{N_unclassified_should1_is1}[$i];
							my $isUnclassified = $v->{N_unclassified_should0_is1}[$i] + $v->{N_unclassified_should1_is1}[$i];
							my $N = $v->{N_unclassified_should0_is0}[$i] + $v->{N_unclassified_should0_is1}[$i] + $v->{N_unclassified_should1_is0}[$i] + $v->{N_unclassified_should1_is1}[$i];
							my $sensitivity = -1;
							if(($v->{N_unclassified_should1_is0}[$i]+$v->{N_unclassified_should1_is1}[$i]) > 0)
							{
								$sensitivity = $v->{N_unclassified_should1_is1}[$i] / ($v->{N_unclassified_should1_is0}[$i]+$v->{N_unclassified_should1_is1}[$i]);
							}
							my $PPV = -1;
							if(($v->{N_unclassified_should1_is1}[$i]+$v->{N_unclassified_should0_is1}[$i]) > 0)
							{
								$PPV = $v->{N_unclassified_should1_is1}[$i] / ($v->{N_unclassified_should1_is1}[$i]+$v->{N_unclassified_should0_is1}[$i]);
							}
							my $specificity = -1;
							if(($v->{N_unclassified_should0_is1}[$i]+$v->{N_unclassified_should0_is0}[$i]) > 0)
							{
								$specificity = $v->{N_unclassified_should0_is0}[$i] / ($v->{N_unclassified_should0_is1}[$i]+$v->{N_unclassified_should0_is0}[$i]);
							}	
							push(@Ns, $N);
							push(@shouldBeUnclassified, $shouldBeUnclassified);
							push(@isUnclassified, $isUnclassified);
							push(@sensitivities, $sensitivity) if($sensitivity != -1);
							push(@PPVs, $PPV) if($PPV != -1);
							push(@specificities, $specificity) if($specificity != -1);											
						}	

						my $averagedOver = scalar(@Ns);
						my $avg_N = (scalar(@Ns)) ? Util::mean(@Ns) : 'NA';
						my $avg_shouldBeUnclassified = (scalar(@shouldBeUnclassified)) ? Util::mean(@shouldBeUnclassified) : 'NA';
						my $avg_isUnclassified = (scalar(@isUnclassified)) ? Util::mean(@isUnclassified) : 'NA';
						my $avg_sensitivity = (scalar(@sensitivities)) ? Util::mean(@sensitivities) : 'NA';
						my $avg_PPVs = (scalar(@PPVs)) ? Util::mean(@PPVs) : 'NA';
						my $avg_specificity = (scalar(@specificities)) ? Util::mean(@specificities) : 'NA';
						
						push(@output_fields_byLevelUnknown, $averagedOver, $avg_N, $avg_shouldBeUnclassified, $avg_isUnclassified, $avg_sensitivity, $avg_PPVs, $avg_specificity);
					}
				}
				print UNCLASSIFIED join("\t", @output_fields_byLevelUnknown), "\n";
				print UNCLASSIFIED_ALL join("\t", $suffix, @output_fields_byLevelUnknown), "\n";
			}
		}
				
		close(UNCLASSIFIED);
	}
	
	{
		open(BARPLOTSREADCAT, '>', $prefix_for_outputFiles . '_forPlot_barplots_readCategory') or die;
		print BARPLOTSREADCAT join("\t", qw/readCategory evaluationLevel method N callRateAvg accuracyAvg accuracyAvgExactltAtLevel callRate_raw accuracy_raw accuracy_raw_exactltyAtLevel/), "\n";
		
		#open(BYREADLENGTH, '>', $prefix_for_outputFiles . '_forPlot_byReadLength') or die;
		#print BYREADLENGTH join("\t", qw/readCategory evaluationLevel method readLength callRateAvg accuracyAvg accuracyAvgExactltAtLevel/), "\n";

		my %categories_attachment_forPrint = map {'attachedTo_' . $_ => 1} qw/species genus family superfamily Unclassified/;
		my %values_attachment_forPrint;
		my %values_attachment_forPrint_N;
		
		foreach my $readCategory (sort keys %{$allSimulations_data_href->{callRate_and_accuracy_byReadCategory}})
		{
			foreach my $label (sort keys %{$allSimulations_data_href->{callRate_and_accuracy_byReadCategory}->{$readCategory}})
			{
				foreach my $level (sort keys %{$allSimulations_data_href->{callRate_and_accuracy_byReadCategory}->{$readCategory}{$label}})
				{
					my $v = $allSimulations_data_href->{callRate_and_accuracy_byReadCategory}->{$readCategory}{$label}{$level};
					my @callRates;
					my @accuracies;
					my @accuracies_exactlyAtLevel;
					foreach my $e (@$v)
					{
						push(@callRates, $e->[0]);
						push(@accuracies, $e->[1]);
						push(@accuracies_exactlyAtLevel, $e->[2]);
					}	
					die unless(scalar(@callRates));
					die unless(scalar(@accuracies));
					die unless(scalar(@accuracies_exactlyAtLevel));
					die unless(scalar(@callRates) == scalar(@accuracies));
					die unless(scalar(@callRates) == scalar(@accuracies_exactlyAtLevel));
					my $avg_callRate = Util::mean(@callRates);
					my $avg_accuracy = Util::mean(@accuracies);
					my $avg_accuracy_exactlyAtLevel = Util::mean(@accuracies_exactlyAtLevel);
					print BARPLOTSREADCAT join("\t", $readCategory, $level, $label, scalar(@callRates), $avg_callRate, $avg_accuracy, $avg_accuracy_exactlyAtLevel, join(';', @callRates), join(';', @accuracies), join(';', @accuracies_exactlyAtLevel)), "\n";
					
					my @attachmentHashes = @{$allSimulations_data_href->{attachedTo_byReadCategory}->{$readCategory}{$label}{$level}};				
					foreach my $h (@attachmentHashes)
					{
						foreach my $k (keys %$h)
						{
							$categories_attachment_forPrint{$k}++;						
						}
					}
					
					my %valuesInHashes;
					foreach my $h (@attachmentHashes)
					{
						my $s_nonDirectlyAttached = 0;
						foreach my $k (keys %categories_attachment_forPrint)
						{
							my $v = (exists $h->{$k}) ? $h->{$k} : 0;
														
							push(@{$valuesInHashes{$k}}, $v);
							if($k ne 'attachedToDirectlyMappable')
							{
								$s_nonDirectlyAttached += $v;
							}
						}
						die unless(abs(1 - $s_nonDirectlyAttached) <= 1e-3); 
					}
					
					foreach my $k (keys %valuesInHashes)
					{
						my $v = Util::mean(@{$valuesInHashes{$k}});
						$values_attachment_forPrint{$readCategory}{$label}{$level}{$k} = $v;
						$values_attachment_forPrint_N{$readCategory}{$label}{$level}{$k} = scalar(@{$valuesInHashes{$k}});
						die unless(exists $categories_attachment_forPrint{$k});
					}
					
					# my @k_attachedTo = grep {$_ =~ /attachedTo_/) keys %valuesInHashes;
					# my $n_attachedTo = scalar(@{$valuesInHashes{$k_attachedTo[0]}});
					# foreach my $k (keys %valuesInHashes)
					# {
						# my $n_attachedTo_this = 
					# }
					
					foreach my $rL (sort {$a <=> $b} keys %{$allSimulations_data_href->{callRate_and_accuracy_byReadCategory_byLength}->{$readCategory}{$label}{$level}})
					{
						my $v_rL = $allSimulations_data_href->{callRate_and_accuracy_byReadCategory_byLength}->{$readCategory}{$label}{$level}{$rL};
						my @callRates_rL ;
						my @accuracies_rL ;
						foreach my $e (@$v_rL)
						{
							push(@callRates_rL, $e->[0]);
							push(@accuracies_rL, $e->[1]);
						}			

						die unless(scalar(@callRates_rL));
						die unless(scalar(@accuracies_rL));
						my $avg_callRate_rL = Util::mean(@callRates_rL);
						my $avg_accuracy_rL = Util::mean(@accuracies_rL);
						#print BYREADLENGTH join("\t", $readCategory, $level, $label, $rL, $avg_callRate_rL, $avg_accuracy_rL), "\n";
						
					}
				}
			}
		}
		close(BARPLOTSREADCAT);
		#close(BYREADLENGTH);
		
		open(BARPLOTS_ATTACHEDTO, '>', $prefix_for_outputFiles . '_forPlot_barplots_attachedTo') or die;
		my @keys_attachedTo = sort keys %categories_attachment_forPrint;
		# print BARPLOTS_ATTACHEDTO join("\t", qw/readCategory method/, @keys_attachedTo), "\n";
		print BARPLOTS_ATTACHEDTO join("\t", qw/readCategory method/, map {$_, 'N_' . $_ } @keys_attachedTo), "\n";
		foreach my $readCategory (sort keys %values_attachment_forPrint)
		{
			foreach my $label (sort keys %{$values_attachment_forPrint{$readCategory}})
			{
				foreach my $level ((sort keys %{$values_attachment_forPrint{$readCategory}{$label}})[0])
				{
					# my @values_forPrint = ($readCategory, $level, $label);
					my @values_forPrint = ($readCategory, $label);
					foreach my $k (@keys_attachedTo)
					{
						my $v = 0;
						my $N = 0;
						if(exists $values_attachment_forPrint{$readCategory}{$label}{$level}{$k})
						{
							$v = $values_attachment_forPrint{$readCategory}{$label}{$level}{$k};
							$N = $values_attachment_forPrint_N{$readCategory}{$label}{$level}{$k};
						}
						push(@values_forPrint, $v);
						push(@values_forPrint, $N);
					}
					print BARPLOTS_ATTACHEDTO join("\t", @values_forPrint), "\n";
				}
			}
		}
		close(BARPLOTS_ATTACHEDTO);	
	}
	
	# output data for bar plots (full DB?)
	{
		my %_methods;
		my %_readStratification;
		my %_evaluationLevels;
		foreach my $variety (@varieties)
		{
			foreach my $methodName (keys %{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}})
			{
				$_methods{$methodName}++;
				foreach my $category (keys %{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}})
				{
					$_readStratification{$category}++;
					
					if(exists $allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category})
					{
						foreach my $k (keys %{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category}})
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
		my @evaluationLevels = qw/absolute species genus family/;
		die Dumper("Missing evaluation levels I", \@evaluationLevels, \%_evaluationLevels, [\@varieties]) unless(all {exists $_evaluationLevels{$_}} @evaluationLevels);
		
		open(BARPLOTSFULLDB, '>', $prefix_for_outputFiles . '_forPlot_barplots_fullDB') or die;
		print BARPLOTSFULLDB join("\t", qw/readLevel variety method level callRate accuracy averagedOver/), "\n";
			
		{
			open(READSABSOLUTELYCORRECT, '>', $prefix_for_outputFiles . '_readsAbsolutelyCorrect') or die;
			my @header_fields_1_absolutelyCorrect = ('ReadLevel');
			my @header_fields_2_absolutelyCorrect = ('');
			my @header_fields_3_absolutelyCorrect = ('');

			foreach my $variety (@varieties)
			{
				my $hf2_before = $#header_fields_2_absolutelyCorrect;
				foreach my $method (@methods)
				{
					push(@header_fields_2_absolutelyCorrect, $method, '', '', '', '', '');		
					push(@header_fields_3_absolutelyCorrect, 'nExperiments', 'Ntotal_avg', 'OKtotal_avg', 'NmadeCall_avg', 'OKmadeCall_avg', 'noCall_avg');
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
						my @outer_missing;
						my @outer_NmadeCall;
						my @outer_correct;
						
						my @Ntotal;
						my @percOK_madeCall;
						my @percOK_madeCall_fullAccuracy;
						my @percOK_total;
						my @perc_missing;
						my @callRate;
						
						if(exists $allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}{$readLevel})
						{
							my @components_missing =  @{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}{$readLevel}{missing}};
							my @components_NmadeCall =  @{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}{$readLevel}{N}};
							my @components_correct =  @{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}{$readLevel}{correct}};
							
							die unless(scalar(@components_missing) == scalar(@components_NmadeCall));
							die unless(scalar(@components_NmadeCall) == scalar(@components_correct));
							
							for(my $componentI = 0; $componentI <= $#components_missing; $componentI++)
							{
								my $i_missing = $components_missing[$componentI];
								my $i_NmadeCall = $components_NmadeCall[$componentI];
								my $i_correct = $components_correct[$componentI];
								
								my $i_Ntotal =  $i_missing + $i_NmadeCall;
					
								my $i_percOK_madeCall = sprintf("%.2f", ($i_correct / $i_NmadeCall)) if($i_NmadeCall > 0);
								my $i_percOK_madeCall_fullAccuracy = ($i_correct / $i_NmadeCall) if($i_NmadeCall > 0);
								my $i_percOK_total = sprintf("%.2f", ($i_correct / $i_Ntotal)) if($i_Ntotal > 0);
								my $i_perc_missing = sprintf("%.2f", ($i_missing / ($i_Ntotal))) if($i_Ntotal > 0);						
								my $i_callRate = $i_NmadeCall / $i_Ntotal if($i_Ntotal > 0);	

								if($i_Ntotal > 0)
								{
									push(@Ntotal, $i_Ntotal);
									push(@outer_NmadeCall, $i_NmadeCall);
									push(@percOK_madeCall, $i_percOK_madeCall);
									push(@percOK_total, $i_percOK_total);
									push(@perc_missing, $i_perc_missing);
									push(@callRate, $i_callRate);
									push(@percOK_madeCall_fullAccuracy, $i_percOK_madeCall_fullAccuracy);
								}
							}
				
						}
						
						if(scalar(@Ntotal))
						{
							die unless(scalar(@outer_NmadeCall) == scalar(@Ntotal));
							die unless(scalar(@percOK_madeCall) == scalar(@Ntotal));
							die unless(scalar(@percOK_total) == scalar(@Ntotal));
							die unless(scalar(@perc_missing) == scalar(@Ntotal));
							die unless(scalar(@callRate) == scalar(@Ntotal));
							die unless(scalar(@percOK_madeCall_fullAccuracy) == scalar(@Ntotal));
						}
						
						my $Ntotal = (scalar(@Ntotal)) ? sprintf("%.2f", Util::mean(@Ntotal)) : 'NA';
						my $NmadeCall = (scalar(@Ntotal)) ? sprintf("%.2f", Util::mean(@outer_NmadeCall)) : 'NA';
						my $percOK_madeCall = (scalar(@Ntotal)) ? sprintf("%.2f", Util::mean(@percOK_madeCall)) : 'NA';
						my $percOK_total = (scalar(@Ntotal)) ? sprintf("%.2f", Util::mean(@percOK_total)) : 'NA';
						my $perc_missing = (scalar(@Ntotal)) ? sprintf("%.2f", Util::mean(@perc_missing)) : 'NA';
						my $callRate = (scalar(@Ntotal)) ? sprintf("%.2f", Util::mean(@callRate)) : 'NA';
						my $percOK_madeCall_fullAccuracy = (scalar(@Ntotal)) ? sprintf("%.2f", Util::mean(@percOK_madeCall_fullAccuracy)) : 'NA';
						
						push(@output_fields_absolutelyCorrect, scalar(@Ntotal), $Ntotal, $percOK_total, $NmadeCall, $percOK_madeCall, $perc_missing);
						
						print BARPLOTSFULLDB join("\t", $readLevel, $variety, $methodName, 'absolute', $callRate, $percOK_madeCall_fullAccuracy, scalar(@Ntotal)), "\n";

						print "Generating ", $prefix_for_outputFiles . '_readsAbsolutelyCorrect', " $readLevel $variety $methodName: averaging over ", scalar(@Ntotal), " iterations.\n";
					}
				}
			
				print READSABSOLUTELYCORRECT join("\t", @output_fields_absolutelyCorrect), "\n";			}	
			

			close(READSABSOLUTELYCORRECT);
		}
		
		{
			open(READSCORRECTBYLEVEL, '>', $prefix_for_outputFiles . '_readsCorrectByLevel') or die;
			my @header_fields_1_byLevelCorrect = ('ReadLevel', 'EvaluationLevel');
			my @header_fields_2_byLevelCorrect = ('', '');
			my @header_fields_3_byLevelCorrect = ('', '');	
			
			foreach my $variety (@varieties)
			{
				my $hf2_before = $#header_fields_2_byLevelCorrect;
				foreach my $method (@methods)
				{
					push(@header_fields_2_byLevelCorrect, $method, '', '', '', '', '');							
					push(@header_fields_3_byLevelCorrect, 'nExperiments', 'Ntotal_avg', 'OKtotal_avg', 'NmadeCall_avg', 'OKmadeCall_avg', 'noCall_avg');
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
			
			# unless($assumeHaveHeader_READSCORRECTBYLEVEL_ALL)
			{
				print READSCORRECTBYLEVEL_ALL join("\t", 'Experiment', @header_fields_1_byLevelCorrect), "\n";
				print READSCORRECTBYLEVEL_ALL join("\t", '', @header_fields_2_byLevelCorrect), "\n";
				print READSCORRECTBYLEVEL_ALL join("\t", '', @header_fields_3_byLevelCorrect), "\n";
			}
						
			foreach my $readLevel (@readLevels)
			{
				foreach my $evaluationLevel (@evaluationLevels)
				{
					my @output_fields_byLevelCorrect = ($readLevel, $evaluationLevel);
					foreach my $variety (@varieties)
					{
						foreach my $methodName (@methods)
						{											
							my @N_total;
							my @N_total_truthDefined;
							
							my @N_madeCall;
							my @N_madeCall_truthDefined;
							
							my @correct;
							my @correct_truthDefined;
							my @callRate;
							
							my @percOK_total;
							my @percOK_madeCall_fullAccuracy;
							my @percOK_total_truthDefined;							
							my @percOK_madeCall;
							my @percOK_madeCall_truthDefined;
							my @percMissing;
							
							if(exists $allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$readLevel}{$evaluationLevel})
							{
																
								my @components_N_madeCall = @{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$readLevel}{$evaluationLevel}{N}};
								my @components_N_truthDefined = @{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$readLevel}{$evaluationLevel}{N_truthDefined}};
								my @components_correct = @{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$readLevel}{$evaluationLevel}{correct}};
								my @components_correct_truthDefined = @{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$readLevel}{$evaluationLevel}{correct_truthDefined}};
								my @components_missing = @{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$readLevel}{$evaluationLevel}{missing}};
								die unless(all{scalar(@$_) == scalar(@components_N_madeCall)} (\@components_N_madeCall, \@components_N_truthDefined, \@components_correct, \@components_N_truthDefined, \@components_missing));
								
								for(my $componentI = 0; $componentI <= $#components_missing; $componentI++)
								{
									my $i_N_madeCall = $components_N_madeCall[$componentI];
									my $i_N_madeCall_truthDefined = $components_N_truthDefined[$componentI];
									my $i_correct = $components_correct[$componentI];
									my $i_correct_truthDefined = $components_correct_truthDefined[$componentI];
									my $i_missing = $components_missing[$componentI];
									
									my $i_Ntotal =  $i_missing + $i_N_madeCall;
						
									my $i_percOK_madeCall_fullAccuracy = ($i_correct / $i_N_madeCall) if($i_N_madeCall > 0);
									my $i_percOK_total = ($i_correct / $i_Ntotal) if($i_Ntotal > 0);
									my $i_perc_missing = ($i_missing / $i_Ntotal) if($i_Ntotal > 0);						
									
									# die Dumper("Weird - N is $N, but missing is $missing?", [$readLevel, $evaluationLevel, $variety, $methodName]) unless($missing <= $N);
																		
									my $i_N_total_truthDefined = $i_N_madeCall_truthDefined + $i_missing;									
									
									my $i_callRate = $i_N_madeCall / $i_Ntotal if($i_Ntotal > 0);
									my $i_percOK_total_truthDefined = ($i_correct_truthDefined / $i_N_total_truthDefined) if($i_N_total_truthDefined > 0);
									
									my $i_percOK_madeCall = ($i_correct / $i_N_madeCall) if($i_N_madeCall > 0);
																
																
									my $i_percOK_madeCall_truthDefined = ($i_correct_truthDefined / $i_N_madeCall_truthDefined) if($i_N_madeCall_truthDefined > 0);
									
									my $i_percMissing = ($i_missing / $i_Ntotal) if($i_Ntotal > 0);
									
									push(@N_total, $i_Ntotal);
									push(@N_total_truthDefined, $i_N_total_truthDefined);
									push(@N_madeCall, $i_N_madeCall);
									push(@N_madeCall_truthDefined, $i_N_madeCall_truthDefined);
									push(@callRate, $i_callRate);
									push(@percOK_total, $i_percOK_total);
									push(@percOK_madeCall_fullAccuracy, $i_percOK_madeCall_fullAccuracy);
									push(@percOK_total_truthDefined, $i_percOK_total_truthDefined);
									push(@percOK_madeCall, $i_percOK_madeCall);
									push(@percMissing, $i_percMissing);
									push(@percOK_madeCall_truthDefined, $i_percOK_madeCall_truthDefined);
								}
							}		

							push(@header_fields_3_byLevelCorrect, 'Ntotal', 'OKtotal', 'NmadeCall', 'OKmadeCall', 'noCall');

							if(scalar(@callRate))
							{
								die unless(scalar(@percOK_madeCall_fullAccuracy) == scalar(@callRate));
								die unless(scalar(@N_total) == scalar(@callRate));
								die unless(scalar(@N_total_truthDefined) == scalar(@callRate));
								die unless(scalar(@N_madeCall) == scalar(@callRate));
								die unless(scalar(@N_madeCall_truthDefined) == scalar(@callRate));
								die unless(scalar(@percOK_total) == scalar(@callRate));
								die unless(scalar(@percOK_total_truthDefined) == scalar(@callRate));
								die unless(scalar(@percOK_madeCall) == scalar(@callRate));
								die unless(scalar(@callRate) == scalar(@callRate));
								die unless(scalar(@percOK_madeCall_truthDefined) == scalar(@callRate));
								die unless(scalar(@percMissing) == scalar(@callRate));
							}
						 
							my $callRate = (scalar(@callRate)) ? Util::mean(@callRate) : 'NA';
							my $percOK_madeCall_fullAccuracy = (scalar(@callRate)) ? Util::mean(@percOK_madeCall_fullAccuracy) : 'NA';
							
							my $N_total = (scalar(@callRate)) ? Util::mean(@N_total) : 'NA';
							my $N_total_truthDefined = (scalar(@callRate)) ? Util::mean(@N_total_truthDefined) : 'NA';
							
							my $N_madeCall = (scalar(@callRate)) ? Util::mean(@N_madeCall) : 'NA';
							my $N_madeCall_truthDefined = (scalar(@callRate)) ? Util::mean(@N_madeCall_truthDefined) : 'NA';						 	
							
							my $percOK_total = (scalar(@callRate)) ? Util::mean(@percOK_total) : 'NA';
							my $percOK_total_truthDefined = (scalar(@callRate)) ? Util::mean(@percOK_total_truthDefined) : 'NA';
							my $percOK_madeCall = (scalar(@callRate)) ? Util::mean(@percOK_madeCall) : 'NA';
							my $percOK_madeCall_truthDefined = (scalar(@callRate)) ? Util::mean(@percOK_madeCall_truthDefined) : 'NA';
							my $percMissing = (scalar(@callRate)) ? Util::mean(@percMissing) : 'NA';
							
							push(@output_fields_byLevelCorrect,
								scalar(@callRate),
								($N_total ne $N_total_truthDefined) ? join(' / ', $N_total, $N_total_truthDefined) : $N_total,
								($percOK_total ne $percOK_total_truthDefined) ? join(' / ', $percOK_total, $percOK_total_truthDefined) : $percOK_total, 
								($N_madeCall ne $N_madeCall_truthDefined) ? join(' / ', $N_madeCall, $N_madeCall_truthDefined) : $N_madeCall,
								($percOK_madeCall ne $percOK_madeCall_truthDefined) ? join(' / ', $percOK_madeCall, $percOK_madeCall_truthDefined) : $percOK_madeCall,
								$percMissing
							);
							
							print "Generating ", $prefix_for_outputFiles . '_forPlot_barplots_fullDB', " $readLevel $variety $methodName: averaging over ", scalar(@callRate), " iterations.\n";
						 
							if($evaluationLevel ne 'absolute')
							{
								print BARPLOTSFULLDB join("\t",  $readLevel, $variety, $methodName, $evaluationLevel, $callRate, $percOK_madeCall_fullAccuracy, scalar(@N_total)), "\n";						
							}
						}
					}
					print READSCORRECTBYLEVEL join("\t", @output_fields_byLevelCorrect), "\n";
					print READSCORRECTBYLEVEL_ALL join("\t", $suffix, @output_fields_byLevelCorrect), "\n";
				}
			}
			
			close(READSCORRECTBYLEVEL);
		}
		
		close(BARPLOTSFULLDB);

	}
	

	# print frequency evaluation text tables and data for XY plots
	{
		my %_methods;
		my %_readStratification;
		my %_evaluationLevels;
		foreach my $variety (@varieties)
		{
			foreach my $methodName (keys %{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}})
			{
				$_methods{$methodName}++;
				foreach my $level (keys %{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}})
				{
					$_evaluationLevels{$level}++;
				}
			}
		}
		
		my @methods = sort keys %_methods;
		 
		my @evaluationLevels = qw/absolute species genus family/;
		# my @evaluationLevels = sort keys %_evaluationLevels;
		die Dumper("Missing evaluation levels II", \@evaluationLevels, \%_evaluationLevels, \@varieties, \%_methods, $allSimulations_data_href->{freq_byVariety_byLevel}) unless(all {exists $_evaluationLevels{$_}} @evaluationLevels);
			
		my %unclassified_is_shouldBe;
				
		{		
			my @evaluateAccuracyAtLevels = validation::getEvaluationLevels();
			
			my %_getLightning_cache;
			my $getLightning = sub {
				my $taxonID = shift;
				my $directlyMappable_href = shift;
				die unless(defined  $directlyMappable_href);
				
				if(exists $_getLightning_cache{$taxonID})
				{
					return $_getLightning_cache{$taxonID};
				}
				else
				{
					my $lightning = validation::getAllRanksForTaxon_withUnclassified($fullTaxonomy_simulation, $taxonID, $directlyMappable_href);
					$_getLightning_cache{$taxonID} = $lightning;
					return $lightning;
				}
			};
			
			my $get_taxonID_category = sub {
				my $taxonID = shift;
				my $directlyMappable_href = shift;
				die unless(defined  $directlyMappable_href);
								
				my $taxonID_lightning = $getLightning->($taxonID, $directlyMappable_href);
				
				my $shouldBeAssignedTo;
				RANK: foreach my $rank (@evaluateAccuracyAtLevels)
				{
					die unless(defined $taxonID_lightning->{$rank});
					die if($taxonID_lightning->{$rank} eq 'NotLabelledAtLevel');
					if(($taxonID_lightning->{$rank} ne 'Unclassified') and ($taxonID_lightning->{$rank} ne 'NotLabelledAtLevel'))
					{
						die if($taxonID_lightning->{$rank} eq 'Undefined');
						$shouldBeAssignedTo = $rank;
						last RANK;
					}
				}
				$shouldBeAssignedTo = '>' . $evaluateAccuracyAtLevels[$#evaluateAccuracyAtLevels] unless(defined $shouldBeAssignedTo);
				die Dumper("Can't find target assignment for $taxonID", $taxonID_lightning) unless(defined $shouldBeAssignedTo);
				return $shouldBeAssignedTo;
			};
			
	
			open(XYPLOTS, '>', $prefix_for_outputFiles . '_forPlot_frequencies_xy') or die;
			print XYPLOTS join("\t", qw/simulationI variety method level taxonID taxonLabel taxonIDCategory isMappable freqTarget freqIs proportionNovelDirectly proportionNovelTotal/), "\n";
			# die Dumper($allSimulations_data_href->{frequencyComparisons_bySimulation}->[0], $allSimulations_data_href->{realizedN});
			
			for(my $simulationI = 0; $simulationI < $allSimulations_data_href->{realizedN}; $simulationI++)  
			{
				next unless(defined $allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]);
				foreach my $variety (keys %{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]})
				{
					(my $variety_forStore = $variety) =~ s/_\d+$//;		
				
					foreach my $label (keys %{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}})
					{
						# die Dumper([keys %{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}}]);
						foreach my $level (keys %{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}})
						{
							my $unclassified_shouldBe = 0;						
							my $unclassified_is = 0;

							if(exists $allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}{$level}{'Unclassified'}) 
							{
								$unclassified_shouldBe = $allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}{$level}{'Unclassified'}[0];
								$unclassified_is = $allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}{$level}{'Unclassified'}[1];
							}	
							
							push(@{$unclassified_is_shouldBe{$variety_forStore}{$label}{$level}}, [$unclassified_shouldBe, $unclassified_is]);
							
						}
					}
					
					next unless(exists $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety});
					die unless(defined $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety});
					die unless(defined $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{mappable});
					die unless(defined $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{truthReadsNovelTree});
					die unless(defined $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{nReads});
					foreach my $label (keys %{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}})
					{
						foreach my $level (keys %{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}})
						{
							next if($level eq 'absolute');
							foreach my $taxonID (keys %{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}{$level}})
							{
								my $taxonID_label = (($taxonID eq 'Unclassified') or ($taxonID eq 'NotLabelledAtLevel')) ? $taxonID : taxTree::taxon_id_get_name($taxonID, $fullTaxonomy_simulation);
								die unless(defined $allSimulations_data_href->{directlyMappable_bySimulation}->[$simulationI]{$variety});
								my $taxonIDCategory = (($taxonID eq 'Unclassified') or ($taxonID eq 'NotLabelledAtLevel')) ? $taxonID : $get_taxonID_category->($taxonID, $allSimulations_data_href->{directlyMappable_bySimulation}->[$simulationI]{$variety});
								my $isMappable = $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{mappable}{$taxonID} ? 1 : 0;
								
								my $have_truth_freq_data = (exists $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{truthReadsNovelTree}{$taxonID});
								my $controlFrequency = (not $have_truth_freq_data) ? 0 : ($allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{truthReadsNovelTree}{$taxonID}[0] / $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{nReads});
								my $freqNovelNew = (not $have_truth_freq_data) ? 0 : ($allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{truthReadsNovelTree}{$taxonID}[1] / $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{nReads});
								my $freqNovelTotal = (not $have_truth_freq_data) ? 0 : ($allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{truthReadsNovelTree}{$taxonID}[2] / $allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{nReads});
								
								if(($taxonID ne 'Unclassified') and ($taxonID ne 'NotLabelledAtLevel') and ($level ne 'definedAndHypotheticalGenomes') and ($level ne 'definedGenomes'))
								{
									unless (abs($controlFrequency - $allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}{$level}{$taxonID}[0]) <= 1e-4)
									{
										die Dumper("Frequency mismatch -- $controlFrequency vs $allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}{$level}{$taxonID}[0]", $simulationI, $variety, $label, $level, $taxonID) ;
										#warn Dumper($allSimulations_data_href->{frequencyComparisons_details_bySimulation}->[$simulationI]{$variety}{truthReadsNovelTree});
									}
								}
								
								print XYPLOTS join("\t",
									$simulationI,
									$variety,
									$label,
									$level,
									$taxonID,
									$taxonID_label,
									$taxonIDCategory,
									$isMappable,
									@{$allSimulations_data_href->{frequencyComparisons_bySimulation}->[$simulationI]{$variety}{$label}{$level}{$taxonID}},
									$freqNovelNew,
									$freqNovelTotal),
									#$controlFrequency),
									"\n";
							}
						}
					}
				}
			}
			close(XYPLOTS);
		}				
		
		{
				
			open(FREQEVALUATION, '>', $prefix_for_outputFiles . '_frequenciesCorrectByLevel') or die;
			my @header_fields_1_freqCorrect = ('EvaluationLevel');
			my @header_fields_2_freqCorrect = ('');
			my @header_fields_3_freqCorrect = ('');
			
			foreach my $variety (@varieties)
			{
				my $hf2_before = $#header_fields_2_freqCorrect;
				foreach my $method (@methods)
				{
					push(@header_fields_2_freqCorrect, $method, '', '');	
					push(@header_fields_3_freqCorrect, 'nExperiments', 'L1_avg', 'r2_avg');					
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
		

			# unless($assumeHaveHeader_FREQEVALUATION_ALL)
			{
				print FREQEVALUATION_ALL join("\t", 'Experiment', @header_fields_1_freqCorrect), "\n";
				print FREQEVALUATION_ALL join("\t", '', @header_fields_2_freqCorrect), "\n";
				print FREQEVALUATION_ALL join("\t", '', @header_fields_3_freqCorrect), "\n";
			}
			
			

			foreach my $evaluationLevel (@evaluationLevels)
			{		
				my @output_fields_freqCorrect = ($evaluationLevel);	
			 
				foreach my $variety (@varieties)
				{				
					foreach my $methodName (@methods)
					{		
						my $n = 0;
						my $M_freqOK = 'NA';
						my $M_AVGRE = 'NA';
						my $M_RRMSE = 'NA';
						my $M_L1 = 'NA';
						my $M_L2 = 'NA';
						my $M_r2 = 'NA';
						
						if(exists $allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel})
						{
							# die Dumper(keys %{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}});
							$n = scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{freqOK}});
							die Dumper("Count mismatch", scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{freqOK}}), scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{L1}})) unless(scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{freqOK}}) == scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{L1}}));
							
							# $freqOK = $allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{correct}/$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{total};
							$M_freqOK = sum(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{freqOK}}) / scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{freqOK}});
							$M_AVGRE = sum(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{AVGRE}}) / scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{AVGRE}});
							$M_RRMSE = sum(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{RRMSE}}) / scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{RRMSE}});
							$M_L1 = sum(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{L1}}) / scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{L1}});
							$M_L2 = sum(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{L2}}) / scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{L2}});
							$M_r2 = sum(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{r2}}) / scalar(@{$allSimulations_data_href->{freq_byVariety_byLevel}->{$variety}{$methodName}{$evaluationLevel}{r2}});
						}
						
						push(@output_fields_freqCorrect, $n, $M_L1, $M_r2);
					}				
				}	
				
				print FREQEVALUATION join("\t", @output_fields_freqCorrect), "\n";
				print FREQEVALUATION_ALL join("\t", $suffix, @output_fields_freqCorrect), "\n";
			}
			
			close(FREQEVALUATION);
		}
		
		
		
		{
			my $fn_unclassified_frequencues =  $prefix_for_outputFiles . '_unclassifiedSummary_frequencies';
			open(UNCLASSIFIED_FREQ, '>', $fn_unclassified_frequencues) or die;				 
		 
			my @header_fields_1_freqUnknown = ('EvaluationLevel');
			my @header_fields_2_freqUnknown = ('');
			my @header_fields_3_freqUnknown = ('');
			
			foreach my $variety (@varieties)
			{
				my $hf2_before = $#header_fields_2_freqUnknown;
				foreach my $method (@methods)
				{
					push(@header_fields_2_freqUnknown, $method, '', '', '');	
					push(@header_fields_3_freqUnknown, 'averagedOver', 'avgUnknownTarget', 'avgUnknownIs', 'avgUnknownFreqDiff');					
				}
				
				my $hf2_after = $#header_fields_2_freqUnknown;
				my $requiredFields = $hf2_after - $hf2_before;
				die unless($requiredFields > 0);
				my @addToHeader1 = ($variety, (('') x ($requiredFields - 1)));
				die unless(scalar(@addToHeader1) == $requiredFields);
				push(@header_fields_1_freqUnknown, @addToHeader1);
			}
			
			print UNCLASSIFIED_FREQ join("\t", @header_fields_1_freqUnknown), "\n";
			print UNCLASSIFIED_FREQ join("\t", @header_fields_2_freqUnknown), "\n";
			print UNCLASSIFIED_FREQ join("\t", @header_fields_3_freqUnknown), "\n";
		

			# unless($assumeHaveHeader_FREQEVALUATION_ALL)
			{
				print UNCLASSIFIED_FREQ_ALL join("\t", 'Experiment', @header_fields_1_freqUnknown), "\n";
				print UNCLASSIFIED_FREQ_ALL join("\t", '', @header_fields_2_freqUnknown), "\n";
				print UNCLASSIFIED_FREQ_ALL join("\t", '', @header_fields_3_freqUnknown), "\n";
			}
			

			foreach my $evaluationLevel ('definedGenomes', @evaluationLevels)
			{		
				next if($evaluationLevel eq 'definedAndHypotheticalGenomes');
				next if($evaluationLevel eq 'absolute');			
				my @output_fields_freqCorrect = ($evaluationLevel);	
			 
				foreach my $variety (@varieties)
				{				
					foreach my $methodName (@methods)
					{		
						my @all_shouldbe;
						my @all_is;
						my @all_diff;
						
						if(defined $unclassified_is_shouldBe{$variety}{$methodName}{$evaluationLevel})
						{
							my $averagedOver = scalar(@{$unclassified_is_shouldBe{$variety}{$methodName}{$evaluationLevel}});
							
							for(my $i = 0; $i < $averagedOver; $i++)
							{ 
								push(@all_shouldbe, $unclassified_is_shouldBe{$variety}{$methodName}{$evaluationLevel}[$i][0]);
								push(@all_is, $unclassified_is_shouldBe{$variety}{$methodName}{$evaluationLevel}[$i][1]);
								push(@all_diff,  abs($unclassified_is_shouldBe{$variety}{$methodName}{$evaluationLevel}[$i][0] - $unclassified_is_shouldBe{$variety}{$methodName}{$evaluationLevel}[$i][1]));
							}
						}
						
						my $avg_shouldBe = (scalar(@all_shouldbe)) ? Util::mean(@all_shouldbe) : 'NA';
						my $avg_is = (scalar(@all_is)) ? Util::mean(@all_is) : 'NA';
						my $avg_diff =  (scalar(@all_diff)) ? Util::mean(@all_diff) : 'NA';
					
						push(@output_fields_freqCorrect, scalar(@all_shouldbe), $avg_shouldBe, $avg_is, $avg_diff);
					}				
				}	
				
				print UNCLASSIFIED_FREQ join("\t", @output_fields_freqCorrect), "\n";
				print UNCLASSIFIED_FREQ_ALL join("\t", $suffix, @output_fields_freqCorrect), "\n";
			}
			
			close(UNCLASSIFIED_FREQ);
		}
		
	}
	
	{
		
		my %n_reads_correct_byVariety_byLevel_byLength;

		for(my $jobI = 0; $jobI <= $#{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel_bySimulation}}; $jobI++)
		{
			my $n_reads_correct_byVariety_byLevel_local_href = $allSimulations_data_href->{n_reads_correct_byVariety_byLevel_bySimulation}->[$jobI];
			my $n_reads_correct_byVariety_byLevel_byLength_local_href = $allSimulations_data_href->{n_reads_correct_byVariety_byLevel_byLength_bySimulation}->[$jobI];
			die unless(defined $n_reads_correct_byVariety_byLevel_byLength_local_href);
		
			foreach my $variety (keys %$n_reads_correct_byVariety_byLevel_local_href)
			{		
				my $variety_forStore = $variety;
				$variety_forStore =~ s/_\d+$//;							
						
				foreach my $label (keys %{$n_reads_correct_byVariety_byLevel_local_href->{$variety}})
				{		
					foreach my $category (keys %{$n_reads_correct_byVariety_byLevel_local_href->{$variety}{$label}})
					{	
						die unless(defined $n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category});
						foreach my $level (keys %{$n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}})
						{
							foreach my $rL (keys %{$n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}{$level}})
							{						
								my $d = $n_reads_correct_byVariety_byLevel_byLength_local_href->{$variety}{$label}{$category}{$level}{$rL};
								die Dumper('N', $level, $rL, $d) unless(exists $d->{N});
								die Dumper('missing', $level, $rL, $d) unless(exists $d->{missing});
								die unless(exists $d->{correct});				
								
								my $N = $d->{N} + $d->{missing}; die unless($d > 0);
								my $CR = $d->{N} / $N; die unless(($CR >= 0) and ($CR <= 1));
								my $accuracy = ($d->{N} > 0) ? ($d->{correct} / $d->{N}) : -1; die unless(($accuracy >= -1) and ($accuracy <= 1));
								
								
								push(@{$n_reads_correct_byVariety_byLevel_byLength{$variety_forStore}{$label}{$category}{$level}{$rL}{N}}, $d->{N});
								push(@{$n_reads_correct_byVariety_byLevel_byLength{$variety_forStore}{$label}{$category}{$level}{$rL}{totalN}}, $N);
								push(@{$n_reads_correct_byVariety_byLevel_byLength{$variety_forStore}{$label}{$category}{$level}{$rL}{CR}}, $CR);
								push(@{$n_reads_correct_byVariety_byLevel_byLength{$variety_forStore}{$label}{$category}{$level}{$rL}{accuracy}}, $accuracy);
							}				
						}
					}
				}
			}
		}
		
		# print data for read-length plot
		
		open(BYREADLENGTH_FULLDB, '>', $prefix_for_outputFiles . '_forPlot_byReadLength_fullDB') or die;
		print BYREADLENGTH_FULLDB join("\t", qw/variety readCategory evaluationLevel method readLength averagedOver Ntotal callRateAvg Ncalled accuracyAvg/), "\n";
		foreach my $variety (sort keys %n_reads_correct_byVariety_byLevel_byLength)
		{	
			next unless($variety eq 'fullDB');
			foreach my $label (sort keys %{$n_reads_correct_byVariety_byLevel_byLength{$variety}})
			{		
				foreach my $category (sort keys %{$n_reads_correct_byVariety_byLevel_byLength{$variety}{$label}})
				{
					foreach my $level (sort keys %{$n_reads_correct_byVariety_byLevel_byLength{$variety}{$label}{$category}})
					{
						foreach my $rL (sort {$a <=> $b} keys %{$n_reads_correct_byVariety_byLevel_byLength{$variety}{$label}{$category}{$level}})
						{						
							my $d = $n_reads_correct_byVariety_byLevel_byLength{$variety}{$label}{$category}{$level}{$rL};
							die Dumper('totalN', $level, $rL, $d) unless(exists $d->{totalN});
							die Dumper('CR', $level, $rL, $d) unless(exists $d->{CR});
							die Dumper('accuracy', $level, $rL, $d) unless(exists $d->{accuracy});
							
							die unless(scalar(@{$d->{N}}));
							die unless(scalar(@{$d->{totalN}}));
							die unless(scalar(@{$d->{CR}}));
							die unless(scalar(@{$d->{accuracy}}));
							die unless(scalar(@{$d->{N}}) == scalar(@{$d->{accuracy}}));
							
							my $N = Util::mean(@{$d->{N}}); die unless($N >= 0);
							my $totalN = Util::mean(@{$d->{totalN}}); die unless($totalN >= 0);
							my $CR = Util::mean(@{$d->{CR}}); die unless(($CR >= 0) and ($CR <= 1));
							my $accuracy = Util::mean(@{$d->{accuracy}}); die unless(($accuracy >= -1) and ($accuracy <= 1));
							 
							die unless(scalar(@{$d->{totalN}}) == scalar(@{$d->{CR}}));
							die unless(scalar(@{$d->{CR}}) == scalar(@{$d->{accuracy}}));
							
							print BYREADLENGTH_FULLDB join("\t",
								$variety,
								$category,
								$level,
								$label,
								$rL,
								scalar(@{$d->{N}}),
								$totalN,
								$CR,
								$N,
								$accuracy
							), "\n";
						}
					}
				}
			}
		}
		close(BYREADLENGTH_FULLDB);
	}
	
	if(1 == 0)
	{
		foreach my $variety (@varieties)
		{
			print $variety, "\n";
			print "\tReads correctly assigned:\n";
			foreach my $methodName (keys %{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}})
			{
				print "\t\t", $methodName, "\n";
				foreach my $category (keys %{$allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}})
				{
					print "\t\t\t", $category, "\n";
					print "\t\t\t\t", "Missing:  ", $allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}{$category}{missing}, "\n";
					print "\t\t\t\t", "N      :  ", $allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}{$category}{N}, "\n";
					print "\t\t\t\t", "Correct:  ", $allSimulations_data_href->{n_reads_correct_byVariety}->{$variety}{$methodName}{$category}{correct}, "\n";
				}
			}
			
			print "\tReads correctly assigned, by level:\n";
			foreach my $methodName (keys %{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}})
			{
				print "\t\t", $methodName, "\n";
				foreach my $category (keys %{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}})
				{
					print "\t\t\t", $category, "\n";			
					foreach my $level (keys %{$allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category}})
					{
						print "\t\t\t\t", $level, "\n";			
						print "\t\t\t\t\t", "Missing   :  ", $allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category}{$level}{missing}, "\n";
						print "\t\t\t\t\t", "N         :  ", $allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category}{$level}{N}, "\n";
						print "\t\t\t\t\t", "N_td      :  ", $allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category}{$level}{N_truthDefined}, "\n";
						print "\t\t\t\t\t", "Correct   :  ", $allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category}{$level}{correct}, "\n";
						print "\t\t\t\t\t", "Correct_tD:  ", $allSimulations_data_href->{n_reads_correct_byVariety_byLevel}->{$variety}{$methodName}{$category}{$level}{correct_truthDefined}, "\n";
					}
				}
			}		
		}
	}	
}

1;
