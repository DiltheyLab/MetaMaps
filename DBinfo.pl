#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;
use File::Copy;

$| = 1;

use taxTree;
use Util;

unless(scalar(@ARGV) == 2)
{
	print_help();
}

my $DB = $ARGV[0];
my $reportLevel = $ARGV[1];

my %taxonID_2_contigs;
my %contigLength;
Util::read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLength);

my $taxonomyDir = $DB . '/taxonomy';
my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

my %report;

foreach my $contigID (keys %contigLength)
{
	my $taxonID = Util::extractTaxonID($contigID, '?', '?');
	my @ancestors = ($taxonID, taxTree::get_ancestors($taxonomy, $taxonID));
	
	my %taxonID_ranks;
	foreach my $nodeID (@ancestors)
	{
		my $rank = $taxonomy->{$nodeID}{rank};
		die unless(defined $rank);
		$taxonID_ranks{$rank} = taxTree::taxon_id_get_name($nodeID, $taxonomy);
	}
	
	$taxonID_ranks{$reportLevel} = 'Undefined' unless(defined $taxonID_ranks{$reportLevel});
	
	$report{$taxonID_ranks{$reportLevel}}[0]++;
	$report{$taxonID_ranks{$reportLevel}}[1] += $contigLength{$contigID};
}

print "\nDB statistics at level '$reportLevel':\n";
foreach my $v (sort keys %report)
{
	print "\t - ${v}: $report{$v}[0] genomes, ", sprintf("%.2f", $report{$v}[1] / (1024**2)), "mb.\n";
}
print "\n";

sub print_help
{
	print qq(
DBinfo.pl

  Print some statistics on database composition (# genomes, megabytes).
  
Usage:

  perl DBinfo.pl dbNAME taxonomyLevel
  
Example:

  perl DBinfo.pl databases/miniSeq superfamily
  
  
);
exit;
}