package simulation;

use strict;
use Data::Dumper;
use List::Util qw/all/;
use List::MoreUtils qw/mesh/;
use taxTree;

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
		if($taxonID eq 0)
		{
			foreach my $l ('EqualCoverageUnit', keys %relevantLevels)
			{
				$thisTaxon_levels{$l} = 0;		
			}
		}
		else
		{
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

1;