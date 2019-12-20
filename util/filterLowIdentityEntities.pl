use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../perlLib"; 
use File::Basename;
use Getopt::Long;

$| = 1;

use Util;
use taxTree;
use simulation;
use validation;

$| = 1;

my $DB = 'miniSeq+H';
my $mappings = 'tmp/CAMI_metamaps';
my $identityThreshold = 0.8;
GetOptions (
	'DB:s' => \$DB, 
	'mappings:s' => \$mappings, 
	'identityThreshold:s' => \$identityThreshold, 
);

unless($DB and $mappings)
{
	print_help();
}

my $taxonomy_dir = 'databases/' . $DB . '/taxonomy';
unless(-d $taxonomy_dir)
{
	die "Directory $taxonomy_dir not found - did you specify a correct --DB and are you calling me from the MetaMaps main directory?";
}

unless(-e $mappings)
{
	die "File $mappings not found";
}

unless(($identityThreshold >= 0) and ($identityThreshold <= 1))
{
	die "Please provide --identityThreshold between 0 and 1";
}
$identityThreshold *= 100;

my $WIMP_file = $mappings . '.EM.WIMP';
unless(-e $WIMP_file)
{
	die "Expected file $WIMP_file not found - have you carried out classification?";
}

my $read2Taxon_file = $mappings . '.EM.reads2Taxon';
unless(-e $read2Taxon_file)
{
	die "Expected file $read2Taxon_file not found - have you carried out classification?";
}

my $n_reads_unmapped = getUnmappedReads($WIMP_file);

my @identities = readIdentityDistribution($mappings);
@identities = sort @identities;
die if(scalar(@identities) == 0);

my $median_index = int(scalar(@identities)/2);
die unless(($median_index >= 0) and ($median_index <= $#identities));
my $median_identity = $identities[$median_index];

my $outputFn_extractedIdentities = $mappings . '.extractedIdentities';
open(IDENTITIES, '>', $outputFn_extractedIdentities) or die "Can't open $outputFn_extractedIdentities for writing";
print IDENTITIES join("\n", @identities), "\n";
close(IDENTITIES);

print "\nMedian identity $median_identity - full distribution printed to $outputFn_extractedIdentities\n";

my $n_filter_howMany = scalar(grep {$_ <= $identityThreshold} @identities);

print "\nWith --identityThreshold $identityThreshold -- $n_filter_howMany of ", scalar(@identities), " *individual* best mappings are below that threshold.\n";

my $DB_taxonomy = taxTree::readTaxonomy($taxonomy_dir);
my $mappingUnit_identities = getIdentitiesForMappingUnits($mappings . '.EM', $DB_taxonomy);

my %mappingUnit_2_remove;
print "Hit database genomes:\n";
foreach my $mappingUnit (keys %$mappingUnit_identities)
{
	my @mappingUnit_identities = sort @{$mappingUnit_identities->{$mappingUnit}};
	@mappingUnit_identities = sort @mappingUnit_identities;
	die if(scalar(@mappingUnit_identities) == 0);	
		
	my $mappingUnit_median_index = int(scalar(@mappingUnit_identities)/2);
	die unless(($mappingUnit_median_index >= 0) and ($mappingUnit_median_index <= $#mappingUnit_identities));
	my $mappingUnit_median_identity = $mappingUnit_identities[$mappingUnit_median_index];
	die unless($mappingUnit_median_identity);
	
	print " - ", $mappingUnit, "  ", $mappingUnit_median_identity, "\n";
	if($mappingUnit_median_identity < $identityThreshold)
	{
		$mappingUnit_2_remove{$mappingUnit}++;
	}
}
print "Out of ", scalar(keys %$mappingUnit_identities), ", now remove ", scalar(keys %mappingUnit_2_remove), " because their median identity is below $identityThreshold \n";

(my $mappingUnits_href, my $setReads20) = modifyMainMappingFileAndReturnCounts($mappings . '.EM', \%mappingUnit_2_remove);

my $outputFN_WIMP = $mappings . '.EM-filtered.WIMP';

open(OUTPUTWIMP, '>', $outputFN_WIMP) or die "Cannot open $outputFN_WIMP for writing";
print OUTPUTWIMP join("\t", qw/AnalysisLevel taxonID Name Absolute EMFrequency PotFrequency/), "\n";
foreach my $level (qw/definedGenomes species genus family/)
{
	my $totalReads = $n_reads_unmapped;
	my %outputDist;
	$outputDist{0} = $n_reads_unmapped;
	foreach my $taxonID (keys %$mappingUnits_href)
	{
		my $taxonID_at_level;
		if(($level eq 'definedGenomes') or ($taxonID eq '0'))
		{
			$taxonID_at_level = $taxonID;
		}
		else
		{
			my $ancestors_href = taxTree::get_ancestors_by_rank($DB_taxonomy, $taxonID);
			if(exists $ancestors_href->{$level})
			{
				$taxonID_at_level = $ancestors_href->{$level}
			}
			else
			{
				$taxonID_at_level = 0;
			}
		}
		$outputDist{$taxonID_at_level} += $mappingUnits_href->{$taxonID};
		$totalReads += $mappingUnits_href->{$taxonID};
	}
	die if($totalReads == 0);
	
	foreach my $outputTaxonID (keys %outputDist)
	{
		my $outputName = ($outputTaxonID eq '0') ? 'Unclassified' : $DB_taxonomy->{$outputTaxonID}{names}[0];
		print OUTPUTWIMP join("\t", $level, $outputTaxonID, $outputName, $outputDist{$outputTaxonID}, 'NA', $outputDist{$outputTaxonID}/$totalReads), "\n";
	}
}
close(OUTPUTWIMP);

print "\nProduced filtered WIMP file $outputFN_WIMP\n";

my $fn_read2Taxon_output = $mappings . '.EM-filtered.reads2Taxon';
open(R2T_IN, '<', $read2Taxon_file) or die "Cannot open $read2Taxon_file for reading";
open(R2T_OUT, '>', $fn_read2Taxon_output) or die "Cannot open $fn_read2Taxon_output for writing";
my $sanity_check_set2Zero = 0;
while(<R2T_IN>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line);
	die unless(scalar(@line_fields) == 2);
	if($setReads20->{$line_fields[0]})
	{
		print R2T_OUT join("\t", $line_fields[0], 0), "\n";
		$sanity_check_set2Zero++;
	}
	else
	{
		print R2T_OUT $line, "\n";
	}
}
close(R2T_IN);
close(R2T_OUT);

print "\nProduced filtered read2Taxon file $fn_read2Taxon_output [sanity check: $sanity_check_set2Zero entries modified]\n";

print "\n";

sub getUnmappedReads
{
	my $WIMP = shift;
	my $n_unmapped;
	open(WIMP, '<', $WIMP) or die "Cannot open $WIMP";
	while(<WIMP>)
	{
		if($_ = /definedGenomes\t-3\treadsLongEnough_unmapped\t(\d+)\t/)
		{
			die if($n_unmapped);
			$n_unmapped = $1;
		}
	}	
	close(WIMP);
	return $n_unmapped;
}

sub getIdentitiesForMappingUnits
{
	my $mappings = shift;
	my $taxonomy = shift;
	
	my $currentReadID;
	my $currentBestMappingQuality;
	my $bestMapping_identity;
	my $bestMapping_readLine;
		
	my $n_reads_nonFiltered = 0;
	my %processedReadID;
	my %identities_mappingUnits;
	my $flush = sub {
		my @line_fields = split(/ /, $bestMapping_readLine);
		my $readID = $line_fields[0];	
		die if($processedReadID{$readID});
		$processedReadID{$readID}++;
					
		die unless($readID eq $currentReadID);

		my $mappingUnit = $line_fields[5];
		my $taxonID = Util::extractTaxonID($mappingUnit, $mappings, -1);
		unless(exists $taxonomy->{$taxonID})
		{
			die "Taxon ID $taxonID not in taxonomy, have you specified the correct DB?";
		}
		push(@{$identities_mappingUnits{$taxonID}}, $bestMapping_identity);
	};
	
	open(MAPPINGS, '<', $mappings) or die "Cannot open $mappings";
	while(<MAPPINGS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		
		my @line_fields = split(/ /, $line);
		my $readID = $line_fields[0];
		my $identity = $line_fields[12];
		my $mQ = $line_fields[13];
		die unless($mQ);

		if(defined $currentReadID and ($currentReadID ne $readID))
		{
			$flush->();
			$currentBestMappingQuality = undef;			
		}
		
		$currentReadID = $readID;
		if((not defined $currentBestMappingQuality) or ($currentBestMappingQuality < $mQ))
		{
			$currentBestMappingQuality = $mQ;
			$bestMapping_readLine = $line;
			$bestMapping_identity = $identity;
		}
	}
	close(MAPPINGS);
	$flush->();
		
	return \%identities_mappingUnits;
}

sub modifyMainMappingFileAndReturnCounts
{
	my $mappings = shift;
	my $mappingUnits2Remove_href = shift;
	
	my $outputFn_filtered = $mappings;
	die unless($outputFn_filtered =~ /\.EM/);
	$outputFn_filtered =~ s/\.EM/.EM-filtered/;
	
	
	my $currentReadID;
	my $currentBestMappingQuality;
	my $bestMapping_identity;
	my $bestMapping_readLine;
		
	my $n_reads_nonFiltered = 0;
	my %processedReadID;
	my %setReads2Zero;
	my %distribution_definedGenomes;
	my $flush = sub {
		my @line_fields = split(/ /, $bestMapping_readLine);
		my $readID = $line_fields[0];	
		die if($processedReadID{$readID});
		$processedReadID{$readID}++;
		die unless($readID eq $currentReadID);
		
		my $mappingUnit = $line_fields[5];
		my $taxonID = Util::extractTaxonID($mappingUnit, $mappings, -1);
		if(not exists $mappingUnits2Remove_href->{$taxonID})
		{
			print MAPPINGSOUT $bestMapping_readLine, "\n";	
			$distribution_definedGenomes{$taxonID}++;
			$n_reads_nonFiltered++;
		}
		else
		{
			$distribution_definedGenomes{0}++;
			$setReads2Zero{$readID}++;		
		}
	};
	
	open(MAPPINGSOUT, '>', $outputFn_filtered) or die "Cannot open $outputFn_filtered for writing";
	open(MAPPINGS, '<', $mappings) or die "Cannot open $mappings";
	while(<MAPPINGS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		
		my @line_fields = split(/ /, $line);
		my $readID = $line_fields[0];
		my $identity = $line_fields[12];
		my $mQ = $line_fields[13];
		die unless($mQ);

		if(defined $currentReadID and ($currentReadID ne $readID))
		{
			$flush->();
			$currentBestMappingQuality = undef;			
		}
		
		$currentReadID = $readID;
		if((not defined $currentBestMappingQuality) or ($currentBestMappingQuality < $mQ))
		{
			$currentBestMappingQuality = $mQ;
			$bestMapping_readLine = $line;
			$bestMapping_identity = $identity;
		}
	}
	close(MAPPINGS);
	$flush->();
	
	print "\nProduced filtered file $outputFn_filtered (", scalar(keys %setReads2Zero), " of $n_reads_nonFiltered alignments were filtered out).\n";
	
	return (\%distribution_definedGenomes, \%setReads2Zero);
}

sub readIdentityDistribution
{
	my $mappings = shift;
	
	my $currentReadID;
	my $currentBestIdentity;
	
	my @identities;
	
	my %processedRead;
	open(MAPPINGS, '<', $mappings) or die "Cannot open $mappings";
	while(<MAPPINGS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/ /, $line);
		my $readID = $line_fields[0];
		my $identity = $line_fields[12];

		if(defined $currentReadID and ($currentReadID ne $readID))
		{
			die if($processedRead{$currentReadID});
			$processedRead{$currentReadID}++;
			die unless(defined $currentBestIdentity);
			push(@identities, $currentBestIdentity);
			$currentBestIdentity = undef;
		}
		$currentReadID = $readID;
		if((not defined $currentBestIdentity) or ($currentBestIdentity < $identity))
		{
			$currentBestIdentity = $identity;
		}
	}
	close(MAPPINGS);
	if(defined $currentBestIdentity)
	{
		push(@identities, $currentBestIdentity);
	}
	
	return @identities;
}


sub print_help
{
	print qq(
filterLowIdentityEntities.pl

  Filter out WIMP results that have low average identities.
  
Usage:

  perl filterLowIdentityEntities.pl --DB dbNAME --mappings mappingsFILE --identityThreshold FLOAT
  
Example:

  perl filterLowIdentityEntities.pl --DB miniSeq+H --mappings myMappings --identityThreshold 0.8
  
);
exit;
}