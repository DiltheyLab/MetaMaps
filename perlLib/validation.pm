package validation;

use strict;
use Data::Dumper;
use List::Util qw/all/;
use List::MoreUtils qw/mesh/;
use taxTree;

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
					$forReturn{$rank} = 'NotLabelledAtLevel';
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
	die unless(all {exists $masterTaxonomy->{$_}} values %$reads_inferred);

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
				if(($lightning_truth_inDB->{$rank} ne 'Unclassified') and ($lightning_truth_inDB->{$rank} ne 'NotLabelledAtLevel'))
				{
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
		die unless(defined $assignedToRank);
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
	# my %taxonID_across_ranks;
	foreach my $readID (@readIDs)
	{
		my @read_categories = $get_read_categories->($readID);
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

			foreach my $level (@evaluateAccuracyAtLevels)
			{
				unless(defined $n_reads_correct_byLevel{$category}{$level}{N})
				{
					$n_reads_correct_byLevel{$category}{$level}{missing} = 0;
					$n_reads_correct_byLevel{$category}{$level}{N} = 0;
					$n_reads_correct_byLevel{$category}{$level}{N_truthDefined} = 0;
					$n_reads_correct_byLevel{$category}{$level}{correct} = 0;
					$n_reads_correct_byLevel{$category}{$level}{correct_exactlyAtLevel} = 0;
					$n_reads_correct_byLevel{$category}{$level}{correct_truthDefined} = 0;	
				}
				
				unless(defined $n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{N})
				{
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{missing} = 0;
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{N} = 0;
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{correct} = 0;
				}
			}
		}

		if(exists $reads_inferred->{$readID})
		{
			my $trueTaxonID_inUsedDB = $reads_truth_mappingDB->{$readID};		
			my $inferredTaxonID = $reads_inferred->{$readID};
			
			my $lightning_truth = $getLightning->($trueTaxonID_inUsedDB);
			my $lightning_inferred = $getLightning->($inferredTaxonID);
			my $attachedTo_inInference = $get_read_assignedToLevel->($readID);
			
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
			foreach my $level (@evaluateAccuracyAtLevels)
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
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{N}++;
					
					$n_reads_correct_byLevel{$category}{$level}{'attachedTo_' . $attachedTo_inInference}++;
					$n_reads_correct_byLevel{$category}{$level}{'attachedToDirectlyMappable'} += ((exists $mappableTaxonIDs->{$inferredTaxonID}) ? 1 : 0);
					
					if($lightning_truth->{$level} eq $lightning_inferred->{$level})
					{
						$n_reads_correct_byLevel{$category}{$level}{correct}++;
						$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{correct}++;						
						
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
			foreach my $category (@read_categories)
			{
				$n_reads_correct{$category}{missing}++;			
				die unless($n_reads_correct{$category}{missing} <= scalar(@readIDs));
				$n_reads_correct_byLevel_byLength{$category}{'absolute'}{$readLengthBin}{missing}++;
				
				foreach my $level (@evaluateAccuracyAtLevels)
				{
					$n_reads_correct_byLevel{$category}{$level}{missing}++;
					die unless($n_reads_correct_byLevel{$category}{$level}{missing} <= scalar(@readIDs));
					
					$n_reads_correct_byLevel_byLength{$category}{$level}{$readLengthBin}{missing}++;
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
}

sub distributionLevelComparison
{	
	my $masterTaxonomy = shift;
	my $distribution_truth = shift;
	my $distribution_inferred = shift;
	my $label = shift;
	my $external_comparison = shift;
	my $frequencyComparison_href = shift;
	
	foreach my $level ('definedGenomes', 'definedAndHypotheticalGenomes', @evaluateAccuracyAtLevels)
	{
		next unless(defined $distribution_inferred->{$level});
		die unless(defined $distribution_truth->{$level});
		
		my $totalFreq = 0;
		my $totalFreqCorrect = 0;
		foreach my $inferredTaxonID (keys %{$distribution_inferred->{$level}})
		{
			my $isFreq = $distribution_inferred->{$level}{$inferredTaxonID}[1];
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
			my $isFreq = (exists $distribution_inferred->{$level}{$trueTaxonID}) ? $distribution_inferred->{$level}{$trueTaxonID}[1] : 0;
			$S_AVGRE += ( abs($shouldBeFreq - $isFreq) / $shouldBeFreq);
			$S_RRMSE += (($shouldBeFreq - $isFreq) / $shouldBeFreq)**2;
		}
		
		my $L1_sum = 0;
		my $L2_sum = 0;
		my %joint_taxonIDs = map {$_ => 1} ((keys %{$distribution_inferred->{$level}}), (keys %{$distribution_truth->{$level}}));
		foreach my $taxonID (keys %joint_taxonIDs)
		{
			my $isFreq = (exists $distribution_inferred->{$level}{$taxonID}) ? $distribution_inferred->{$level}{$taxonID}[1] : 0;
			my $shouldBeFreq = (exists $distribution_truth->{$level}{$taxonID}) ? $distribution_truth->{$level}{$taxonID} : 0;
			my $L1_diff = abs($isFreq - $shouldBeFreq);
			my $L2_diff = ($isFreq - $shouldBeFreq)**2;
			$L1_sum += $L1_diff;
			$L2_sum += $L2_diff;
			
			if($taxonID eq '1280')
			{
				# print join("\t", "VALUE IN COMPARSION for $label // $level: ", $taxonID, $shouldBeFreq, $isFreq), "\n";
			}
			
			$frequencyComparison_href->{$label}{$level}{$taxonID} = [$shouldBeFreq, $isFreq];
		}
		
		my $L1 = $L1_sum;
		my $L2 = sqrt($L2_sum);
		
		my $AVGRE *= (1 / scalar(keys %{$distribution_inferred->{$level}}));
		my $RRMSE *= (1 / scalar(keys %{$distribution_inferred->{$level}}));
		$RRMSE = sqrt($RRMSE);
		
		print join("\t", $label, $level, $totalFreqCorrect), "\n";
		
		
		if(defined $external_comparison)
		{
			$external_comparison->{$label}{$level}{total} += $totalFreq; 
			$external_comparison->{$label}{$level}{correct} += $totalFreqCorrect; 
			push(@{$external_comparison->{$label}{$level}{AVGRE}}, $AVGRE);
			push(@{$external_comparison->{$label}{$level}{RRMSE}}, $RRMSE);
			push(@{$external_comparison->{$label}{$level}{L1}}, $L1);
			push(@{$external_comparison->{$label}{$level}{L2}}, $L2);
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
				
		die Dumper(\%line) unless(defined $taxonID_master);
		die Dumper("Unknown taxon ID $taxonID_master in file $f", $taxonomy->{$taxonID_master}, $taxonID_nonMaster, \%line) unless(($taxonID_master eq 'Undefined') or ($taxonID_master eq 'Unclassified') or ($taxonID_master eq 'NotLabelledAtLevel') or (defined $taxonomy->{$taxonID_master}));
		die Dumper("Weird line in $f", \%line) unless(defined $line{Absolute});
		die unless(defined $line{PotFrequency});
		$inference{$line{AnalysisLevel}}{$taxonID_master}[0] += $line{Absolute};
		die unless(exists $line{PotFrequency});
		$inference{$line{AnalysisLevel}}{$taxonID_master}[1] += $line{PotFrequency};
		
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
			$taxonID_master = 1;
		}
		
		die "Undefined taxon ID $taxonID_master in file $file $." unless(exists $masterTaxonomy->{$taxonID_master});
		die if(defined $inferred_raw_reads{$readID});
		
		$inferred_raw_reads{$readID} = $taxonID_master;
	}
	close(I);	
	
	return \%inferred_raw_reads;	
}



1;