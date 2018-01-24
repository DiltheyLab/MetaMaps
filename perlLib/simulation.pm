package simulation;

use strict;
use Data::Dumper;
use List::Util qw/all/;
use List::MoreUtils qw/mesh/;
use taxTree;

# this assumes that the taxonomy is complete, i.e. there are no unclassified
# inputs
sub truthReadFrequenciesFromReadCounts
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
				$thisTaxon_levels{$l} = 'Unclassified';
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
			my $name;
			my $taxonIDForPrint = $taxonID;
			if($taxonID eq 'Unclassified')
			{
				$name = $taxonID;
				$taxonIDForPrint = 0;
			}
			else
			{ 
				$name = taxTree::taxon_id_get_name($taxonID, $taxonomy_href) 
			}
			print O join("\t", $level, $taxonIDForPrint, $name, $nReads, $f), "\n";
		}
		
	}
	close(O);
	
}


sub truthGenomeFrequenciesFromReadCounts
{
	my $outputFn = shift;
	my $taxonID_2_bases_href = shift;
	my $readCounts_href = shift;
	my $taxa_genome_lengths_href = shift;
	my $taxonomy_href = shift;
	
	open(O, '>', $outputFn) or die "Cannot open $outputFn";
	print O join("\t", qw/taxonID Name Bases nReads Genomes genomesProportion/), "\n";
	
	my $sum_genomes = 0;
	foreach my $taxonID (keys %$taxonID_2_bases_href)
	{
		$sum_genomes += ($taxonID_2_bases_href->{$taxonID} / $taxa_genome_lengths_href->{$taxonID});
	}
	foreach my $taxonID (sort keys %$taxonID_2_bases_href)
	{
		die unless(defined $taxa_genome_lengths_href->{$taxonID});
		die unless($taxa_genome_lengths_href->{$taxonID});

		my $n_genomes = $taxonID_2_bases_href->{$taxonID} / $taxa_genome_lengths_href->{$taxonID};
		print O join("\t",
			$taxonID, 
			taxTree::taxon_id_get_name($taxonID, $taxonomy_href),
			$taxonID_2_bases_href->{$taxonID},
			$readCounts_href->{$taxonID},
			$n_genomes,
			$n_genomes / $sum_genomes,
		), "\n";		
	}
	
	close(O);
}
1;