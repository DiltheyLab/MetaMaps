package taxTree;

use strict;
use Data::Dumper;
use List::Util qw/all/;
use Storable qw/dclone/;

sub readTaxonomy
{
	my $buildDB_taxonomyDir = shift;
		
	my $taxonomy_names_f = $buildDB_taxonomyDir . '/names.dmp';
	my $taxonomy_nodes_f = $buildDB_taxonomyDir . '/nodes.dmp';
	my $merged_nodes_f = $buildDB_taxonomyDir . '/merged.dmp';
	
	unless(-e $taxonomy_names_f)
	{
		die "File $taxonomy_names_f missing, $buildDB_taxonomyDir is not a valid taxonomy directory.";
	}
	unless(-e $taxonomy_nodes_f)
	{
		die "File $taxonomy_nodes_f missing, $buildDB_taxonomyDir is not a valid taxonomy directory.";
	}
	
	my $nodeid_2_names_href = {};
	
	foreach my $preference ('scientific name', 'genbank common name', undef)
	{
		open(NAMES, '<', $taxonomy_names_f) or die "Cannot open $taxonomy_names_f";
		while(<NAMES>)
		{
			# print "\rNAMES $.       " if(($. % 10000) == 0);
			my $line = $_;
			chomp($line);
			$line =~ s/\t?\|$//;
			my @fields = split(/\t\|\t/, $line, -1);
			if($preference and ($fields[3] ne $preference))
			{
				next;
			}
			if(not defined $nodeid_2_names_href->{$fields[0]})
			{	
				$nodeid_2_names_href->{$fields[0]} = [$fields[1], $fields[2]];
			}
		}
		close(NAMES);
	}
	my $tree_href = {};
	open(NODES, '<', $taxonomy_nodes_f) or die "Cannot open $taxonomy_nodes_f";
	while(<NODES>)
	{
		# print "\rNODES $.       " if(($. % 10000) == 0);	
		my $line = $_;
		chomp($line);
		$line =~ s/\t?\|$//;		
		my @fields = split(/\t\|\t/, $line, -1);
		my $node_id = $fields[0];
		my $parent_node = $fields[1];
		my $rank = $fields[2];
		$tree_href->{$node_id} = [$parent_node, $rank];
	}
	close(NODES);
	
	# print "\n";

	my %full_tree = map {$_ => {}} keys %$tree_href; 
	foreach my $node (keys %$tree_href)
	{
		my $parent_node = $tree_href->{$node}[0];
		my $rank_node = $tree_href->{$node}[1];
		die unless(exists $full_tree{$node});
		die unless(exists $full_tree{$parent_node});
		
		die if(exists $full_tree{$node}{rank});
		$full_tree{$node}{rank} = $rank_node;
		
		if($parent_node ne $node)
		{
			$full_tree{$node}{parent} = $parent_node;
			push(@{$full_tree{$parent_node}{children}}, $node);
		}
	}
	
	foreach my $node (keys %full_tree)
	{
		$full_tree{$node}{children} = [] if(not exists $full_tree{$node}{children});
		$full_tree{$node}{parent} = undef if(not exists $full_tree{$node}{parent});
		if(exists $nodeid_2_names_href->{$node})
		{
			$full_tree{$node}{names} = $nodeid_2_names_href->{$node};
		}
		else
		{
			$full_tree{$node}{names} = [];
		}
	}

	taxonomy_checkConsistency(\%full_tree);
	
	return \%full_tree;
	
	# my @compatible_node_ids = search_leave_names(\%full_tree, "coli");
	# print_node_list(\%full_tree, \@compatible_node_ids);
	
	print "\nTree in memory now.\n";
	
	my @inner_nodes = grep {scalar(@{$full_tree{$_}{children}}) > 1} keys %full_tree;
	
	print "Total nodes ", scalar(keys %full_tree), ", inner nodes ", scalar(@inner_nodes), "\n";
	
	for(my $innerNodeI = 0; $innerNodeI <= $#inner_nodes; $innerNodeI++)
	{
		my $innerNodeID = $inner_nodes[$innerNodeI];
		my $nodeName = join('; ', @{$full_tree{$innerNodeID}{names}});
		print "Children for $innerNodeI / $#inner_nodes, $nodeName \n";
		expected_distance_newNode(\%full_tree, $innerNodeID);
	}
	
	
	my $one_E_coli_leaf = 1182710;
	my $E_coli = 562;
	my $Escherichia = 561;
	
	expected_distance_newNode(\%full_tree, $Escherichia);
	
	exit;
	
	print "Leaves:\n";
	print_node_list(\%full_tree, [$Escherichia, descendants_leaves(\%full_tree, $Escherichia)]);
	print "\n\nAll descendants:\n";
	print_node_list(\%full_tree, [$Escherichia, descendants(\%full_tree, $Escherichia)]);
	
	exit;
	
	print_node_list(\%full_tree, [$one_E_coli_leaf, get_ancestors(\%full_tree, $one_E_coli_leaf)]);
	print "Number of nodes in tree: ", scalar(keys %full_tree), "\n";
	print "OK2\n\n";
	sleep 3600;
	exit;
	
	return ($nodeid_2_names_href, $tree_href);
}

sub cloneTaxonomy_integrateX
{
	my $masterTaxonomy_template = shift;
	my $masterTaxonomy_merged = shift;
	my $taxonomyWithX = shift;
	
	die unless(defined $taxonomyWithX);
	
	my @x_nodes = grep {$_ =~ /^x/} keys %$taxonomyWithX;
	
	print "Now integrate ", scalar(@x_nodes), " into clone of master taxonomy.\n";
	
	my $masterTaxonomy = dclone $masterTaxonomy_template;
	
	foreach my $xNode (@x_nodes)
	{
		my $nodeData = $taxonomyWithX->{$xNode};
		die unless(defined $nodeData);
		die if(exists $masterTaxonomy->{$xNode});
		die unless(scalar(@{$nodeData->{children}}) == 0);
		
		my $parent_in_taxonomyWithX = $nodeData->{parent};
		die unless(defined $taxonomyWithX->{$parent_in_taxonomyWithX});
		
		my $newParent = findCurrentNodeID($masterTaxonomy, $masterTaxonomy_merged, $parent_in_taxonomyWithX);
		my $nodeData_new = dclone $nodeData;
		
		$nodeData_new->{parent} = $newParent;
		push(@{$masterTaxonomy->{$newParent}{children}}, $xNode);
		$masterTaxonomy->{$xNode} = $nodeData_new;
	}
	
	taxonomy_checkConsistency($masterTaxonomy);	
		
	return $masterTaxonomy;
}

sub cloneTaxonomy_removeNodes
{
	my $tree_href = shift;
	my $removeNodes_href = shift;
	
	die unless(all {exists $tree_href->{$_}} keys %$removeNodes_href);
	
	my $reducedTree = dclone $tree_href;
	
	foreach my $removeNode (keys %$removeNodes_href)
	{
		delete $reducedTree->{$removeNode};
	}
	
	foreach my $remainingNodeID (keys %$reducedTree)
	{
		my $remainingNode = $reducedTree->{$remainingNodeID};
		die "Parent node should have been deleted" if(exists $removeNodes_href->{$remainingNode->{parent}});
		@{$remainingNode->{children}} = grep {not exists $removeNodes_href->{$_}} @{$remainingNode->{children}};
	}
	
	taxonomy_checkConsistency($reducedTree);
		
	die unless(all {not exists $reducedTree->{$_}} keys %$removeNodes_href);

	return $reducedTree;
}

sub expected_distance_newNode
{
	my $tree_href = shift;
	my $newNodeParent = shift;
	die unless(defined $tree_href->{$newNodeParent});
	
	my @children = @{$tree_href->{$newNodeParent}{children}};
	unless(scalar(@children))
	{
		die "Please provide a non-leaf node";
	}
	unless(scalar(@children) > 1)
	{
		warn "Non-leaf node, but can't estimate distance";
		return -1;
	}
	
	my %descendants_per_node_cache;
	foreach my $nodeID (@children)
	{
		die unless(defined $tree_href->{$nodeID});
		my @desc;
		if(scalar(@{$tree_href->{$nodeID}{children}}) == 0)
		{
			@desc = ($nodeID);
		}
		else
		{
			@desc = descendants_leaves($tree_href, $nodeID); 
		}	
		$descendants_per_node_cache{$nodeID} = \@desc;
	}
	print "Find new-node expected distance for node $newNodeParent // children $#children \n";
	for(my $childI = 0; $childI <= $#children; $childI++)
	{
		my $childI_nodeID = $children[$childI];
		die unless(defined $descendants_per_node_cache{$childI_nodeID});
		my @childI_corresponding_leaves = @{$descendants_per_node_cache{$childI_nodeID}};
		my @combined_other_leaves;
		for(my $childJ = 0; $childJ <= $#children; $childJ++)
		{				
			next if($childI == $childJ);
			my $childJ_nodeID = $children[$childJ];
			die unless(defined $descendants_per_node_cache{$childJ_nodeID});
			push(@combined_other_leaves, @{$descendants_per_node_cache{$childJ_nodeID}});
		}
		
		# print "For child $childI / $#children, evaluate ", scalar(@childI_corresponding_leaves), " against ", scalar(@combined_other_leaves), " genomes.\n";
	}	
		
}

sub get_leave_ids
{
	my $tree_href = shift;
	my @forReturn;
	foreach my $nodeID (keys %$tree_href)
	{
		if(scalar(@{$tree_href->{$nodeID}{children}}) == 0)
		{
			push(@forReturn, $nodeID);
		}
	}
	return @forReturn;
}

sub taxonomy_checkConsistency
{
	my $taxonomy = shift;
	my $root_nodes = 0;
	my $nI = -1;
	my %children_index;
	foreach my $node (keys %$taxonomy)
	{
		foreach my $childID (@{$taxonomy->{$node}{children}})
		{
			$children_index{$node}{$childID}++;
		}
	}
	foreach my $node (keys %$taxonomy)
	{
		$nI++;
		# print "\rCHECK $nI     ";
		
		my $parentID = $taxonomy->{$node}{parent};
		die unless((defined $parentID) or ($node eq '1'));
		if($parentID)
		{
			my $parent_node = $taxonomy->{$parentID};
			# die unless(scalar(grep {$_ eq $node} @{$parent_node->{children}}) == 1);
			die unless($children_index{$parentID}{$node});
			# die unless(scalar(grep {$_ eq $node} @{$parent_node->{children}}) == 1);
		}
		else
		{
			$root_nodes++;
		}
		
		my @children = @{$taxonomy->{$node}{children}};
		foreach my $childID (@children)
		{
			my $childNode = $taxonomy->{$childID};
			die unless(defined $childNode);
			die unless($childNode->{parent} eq $node);
		}
	}
	print "\n";
	
	die unless($root_nodes == 1);
}

sub removeUnmappableParts
{
	my $taxonomy = shift;
	my $taxonID_mappable_href = shift;
	
	my %keepNode;
	
	foreach my $mappableNode (keys %$taxonID_mappable_href)
	{
		if($taxonID_mappable_href->{$mappableNode} and exists $taxonomy->{$mappableNode})
		{			
			$keepNode{$mappableNode} = 1;
			my $runningNode = $mappableNode;
			while($taxonomy->{$runningNode}{parent})
			{
				$runningNode = $taxonomy->{$runningNode}{parent};
				last if($keepNode{$runningNode});
				$keepNode{$runningNode} = 1;
			}
		}
	}

	print "taxTree::removeUnmappableParts(..): Out of ", scalar(keys %$taxonomy), ", keep ", scalar(keys %keepNode), " nodes.\n";
	
	my @previousNodes = keys %$taxonomy;
	
	foreach my $nodeID (@previousNodes)
	{
		if($keepNode{$nodeID})
		{
			my @node_new_children = grep {$keepNode{$_}} @{$taxonomy->{$nodeID}{children}};
			$taxonomy->{$nodeID}{children} = \@node_new_children;
		}
		else
		{
			delete $taxonomy->{$nodeID};
		}
	}
	
	taxonomy_checkConsistency($taxonomy);
}

sub descendants
{
	my $tree_href = shift;
	my $nodeID = shift;
	die "Unknown taxon ID $nodeID" unless($tree_href->{$nodeID});	
	
	my @forReturn;
	my @left_to_visit = @{$tree_href->{$nodeID}{children}};
	
	while(@left_to_visit)
	{
		my $visitNow = pop(@left_to_visit);
		die unless(defined $tree_href->{$visitNow});
		push(@forReturn, $visitNow);
		push(@left_to_visit, @{$tree_href->{$visitNow}{children}});
	}
	
	return @forReturn;
}
	
sub descendants_leaves
{
	my $tree_href = shift;
	my $nodeID = shift;
	die unless($tree_href->{$nodeID});
	my @descendants = descendants($tree_href, $nodeID);
	my @descendants_leaves = grep {scalar(@{$tree_href->{$_}{children}}) == 0} @descendants;
	return @descendants_leaves;
}

sub print_node_list
{
	my $tree_href = shift;
	my $node_ids_aref = shift;
	
	foreach my $nodeID (@$node_ids_aref)
	{
		die unless(defined $tree_href->{$nodeID}{names});
		print $nodeID, "\n", "\t", join('; ', @{$tree_href->{$nodeID}{names}}), "\n";
	}
}

sub search_leave_names
{
	my $tree_href = shift;
	my $search_string = shift;
	die unless(defined $search_string);
	
	my @forReturn;
	
	foreach my $nodeID (keys %$tree_href)
	{
		next unless(scalar(@{$tree_href->{$nodeID}{children}}) == 0);
		my $compatible =  0;
		foreach my $name (@{$tree_href->{$nodeID}{names}})
		{
			if(index($name, $search_string) != -1)
			{
				$compatible = 1;
			}
		}
		push(@forReturn, $nodeID) if ($compatible);
	}
	
	return @forReturn;
}

sub get_ancestors
{
	my $tree_href = shift;
	my $node = shift;
	die unless(defined $node);
	die unless(exists $tree_href->{$node});
	my $current_id = $node;
	my @forReturn;
	while($tree_href->{$current_id}{parent})
	{	
		push(@forReturn, $tree_href->{$current_id}{parent});
		$current_id = $tree_href->{$current_id}{parent};
	}
	return @forReturn;
}

sub get_ancestors_by_rank
{
	my $tree_href = shift;
	my $node = shift;
	my @ancestor_nodeIDs = get_ancestors($tree_href, $node);
	my %byRank;
	foreach my $ancestorID (@ancestor_nodeIDs)
	{
		my $rank = $tree_href->{$ancestorID}{rank};
		die unless(defined $rank);
		next if($rank eq 'no rank');
		die if(defined $byRank{$rank});
		$byRank{$rank} = $ancestorID;
	}
	return \%byRank;
}



sub get_root_node_id
{
	my $tree_href = shift;
	my @_root_nodes = grep {not defined $tree_href->{$_}{parent}} keys %$tree_href;
	die unless(scalar(@_root_nodes) == 1);
	return $_root_nodes[0];
}

sub get_taxon_id_information
{
	my $id = shift;
	my $taxonomy = shift;
	my $taxonomy_fields_aref = shift;
	
	die "Undefined taxon ID $id" unless(defined $taxonomy->{$id});
	
	my %taxon_data;
	my %taxon_data_ids;
	my $current_id = $id;
	while(1)
	{
		die "Missing tree information for ID $current_id" unless(defined $taxonomy->{$current_id});	
		my $rank = $taxonomy->{$current_id}{rank};
		my $name = taxon_id_get_name($current_id, $taxonomy);
		$taxon_data{$rank} = $name;		
		last unless(defined $taxonomy->{$current_id}{parent});
		$current_id = $taxonomy->{$current_id}{parent};
	}
	
	if($taxonomy_fields_aref)
	{
		my %ret_taxon_data;
		foreach my $field (@$taxonomy_fields_aref)
		{
			if($taxon_data{$field})
			{
				$ret_taxon_data{$field} = $taxon_data{$field};
			}
			else
			{
				$ret_taxon_data{$field} = 'NA';
			}		
		}
		return \%ret_taxon_data;
	}
	else
	{
		return \%taxon_data;
	}
}
	
sub taxon_id_get_name
{
	my $node_id = shift;
	my $taxonomy = shift;
	die "ID $node_id undefined" unless(exists $taxonomy->{$node_id});
	
	my $useName;
	if($taxonomy->{$node_id}{names}->[1])
	{
		$useName = $taxonomy->{$node_id}{names}->[1];
	}
	else
	{
		$useName = $taxonomy->{$node_id}{names}->[0];		
	}
	
	return $useName;
}	


sub readMerged
{
	my $buildDB_taxonomyDir = shift;
	my $merged_nodes_f = $buildDB_taxonomyDir . '/merged.dmp';
	unless(-e $merged_nodes_f)
	{
		die "File $merged_nodes_f missing, but want to read merged nodes";
	}
	
	my %merged;
	
	open(MERGED, '<', $merged_nodes_f) or die "Cannot open $merged_nodes_f";
	while(<MERGED>)
	{
		my $line = $_;
		chomp($line);
		$line =~ s/\t?\|$//;
		my @fields = split(/\t\|\t/, $line, -1);
		die unless(scalar(@fields) == 2);
		
		die if (exists $merged{$fields[0]});
		$merged{$fields[0]} = $fields[1];
	}
	close(MERGED);	
	
	return \%merged;
}

sub findCurrentNodeID
{
	my $taxonomy = shift;
	my $merged = shift;
	my $originalID = shift;
	
	if(($originalID == 0) or ($originalID !~ /^\d+$/))
	{
		return $originalID;
	}
	
	if(exists $taxonomy->{$originalID})
	{
		return $originalID;
	}
	else
	{
		my $runningID = $originalID;
		while(exists $merged->{$runningID})
		{
			$runningID = $merged->{$runningID};
		}
		if(exists $taxonomy->{$runningID})
		{
			return $runningID;
		}
		else
		{
			die "Cannot transate ID $originalID (running $runningID)";
		}
	}
}

sub getNodesForPotentialAttachmentOfNovelSpecies
{
	my $taxonomy_href = shift;
	my @nodeIDs = keys %$taxonomy_href;
	my @nodes_direct_attachment = grep {my $rank = $taxonomy_href->{$_}{rank}; (($rank eq 'species') or ($rank eq 'genus') or ($rank eq 'family'))} @nodeIDs;
	my @full_potential_node_list = map {taxTree::descendants($taxonomy_href, $_)} @nodes_direct_attachment;
	@full_potential_node_list = grep {scalar(@{$taxonomy_href->{$_}{children}}) > 1} @full_potential_node_list;
	my %_u = map {$_ => 1} @full_potential_node_list;
	
	my @forReturn = keys %_u;
	
	my %rankStats;
	my $nodes_multi_rank_children = 0;
	foreach my $nodeID (@forReturn)
	{
		$rankStats{$taxonomy_href->{$nodeID}{rank}}++;
		
		my %ranks_children = map {$taxonomy_href->{$_}{rank} => 1} @{$taxonomy_href->{$nodeID}{children}};
		if(scalar(keys %ranks_children) > 1)
		{
			$nodes_multi_rank_children++;
		}
	}
	
	print "Total nodes considered for attachment (>1 child; in rank <= (species, genus or family); from taxonomy): ", scalar(@forReturn), "\n";
	print "\tOf these, $nodes_multi_rank_children have multi-rank child sets.\n";
	print "\tNode rank stats:\n";
	foreach my $rank (keys %rankStats)
	{
		print "\t\t", $rank, ": ", $rankStats{$rank}, "\n";
	}
	
	return @forReturn;
}

sub getSubComputationsForAttachment
{
	my $taxonomy_href = shift;
	my $nodeID = shift;
	my $mappable_href = shift;
	die unless(defined $mappable_href);
	
	my @node_children = @{$taxonomy_href->{$nodeID}{children}};
	die unless(scalar(@node_children) > 0);
	
	my %mappable_descendants_per_childNode;
	foreach my $childNodeID (@node_children)
	{
		my @nodes_to_consider = ($childNodeID, descendants($taxonomy_href, $childNodeID));
		my @nodes_mappable = grep {exists $mappable_href->{$_}} @nodes_to_consider;
		die Dumper("Node $childNodeID doesn't seem to have mappable children (coming from $nodeID)", \@node_children, \@nodes_to_consider, \@nodes_mappable) unless(scalar(@nodes_mappable) > 0);
		$mappable_descendants_per_childNode{$childNodeID} = \@nodes_mappable;
	}

	my @forReturn;
	foreach my $childNodeID (@node_children)
	{
		my @mappable_descendants = @{$mappable_descendants_per_childNode{$childNodeID}};
		my @otherChildren_mappable_descendants;
		foreach my $childNodeID2 (@node_children)
		{
			next if($childNodeID eq $childNodeID2);
			push(@otherChildren_mappable_descendants, @{$mappable_descendants_per_childNode{$childNodeID2}});
		}
		die unless(scalar(@otherChildren_mappable_descendants));
		
		foreach my $mappable_descendant (@mappable_descendants)
		{
			push(@forReturn, [$childNodeID, $mappable_descendant, \@otherChildren_mappable_descendants]);
			
			# sanity check I`
			die unless(exists $mappable_href->{$mappable_descendant});
		}
			
		# sanity check II
		die unless(all { exists $mappable_href->{$_} } @otherChildren_mappable_descendants);		
	}
	
	return @forReturn;
	
	
}

sub getTaxonomyFileNames
{
	return qw/delnodes.dmp merged.dmp names.dmp nodes.dmp/;
}


sub getRelevantRanks
{
	return qw/species genus family order phylum superkingdom/;
}
	
1;
