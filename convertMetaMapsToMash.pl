#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;
use File::Copy;

$| = 1;

use SimulationsKraken;
use Util;

unless(scalar(@ARGV) == 1)
{
	print_help();
}

my $mash_output_version = `mash`;
die "Mash correctly installed?" unless($mash_output_version =~ /Mash version/);

my $database = $ARGV[0];

# test that database is valid
# todo reinstate
#my %taxonID_2_contigs;
#my %contigLength;
#Util::read_taxonIDs_and_contigs($database, \%taxonID_2_contigs, \%contigLength);

my $DB_fa = $database . '/DB.fa';
die "File $DB_fa not present" unless(-e $DB_fa);

my $outputDir_mash = $database . '/mash';
unless(-d $outputDir_mash)
{
	mkdir($outputDir_mash) or die "Cannot mkdir $outputDir_mash";
}

my $currentTaxonID;
my %sequences_per_taxonID;
open(DB, '<', $DB_fa) or die "Cannot open $DB_fa";
while(<DB>)
{
	my $line = $_;
	next unless($line);
	if(substr($line, 0, 1) eq '>')
	{
		$currentTaxonID = Util::extractTaxonID($line);
	}
	die unless(defined $currentTaxonID);
	$sequences_per_taxonID{$currentTaxonID} .= $line; 
}
close(DB);

my @fastas;
foreach my $taxonID (keys %sequences_per_taxonID)
{
	my $output_fn = $outputDir_mash . '/' . $taxonID . '.fa';
	push(@fastas, $output_fn);
	open(F, '>', $output_fn) or die "Cannot open $output_fn";
	print F $sequences_per_taxonID{$taxonID}, "\n";
	close(F);
}

my $fn_filelist = $outputDir_mash . '/fileList';
open(FILELIST, '>', $fn_filelist) or die "Cannot open $fn_filelist";
print FILELIST join("\n", @fastas), "\n";
close(FILELIST);

my $prefix_sketch = $outputDir_mash . '/DB';
my $cmd_mash = "mash sketch -l $fn_filelist -o $prefix_sketch";
system($cmd_mash) and die "Command failed: $cmd_mash";

foreach my $file (@fastas, $fn_filelist)
{
	unlink($file) or die "Cannot delete $file";
}

print "\n\nProduced ${prefix_sketch}.msh\n\n";

sub print_help
{
	print qq(
convertMetaMapsToMash.pl

  Convert a MetaMaps DB to a Mash sketch (in memory, requires RAM)!
  
Usage:

  perl convertMetaMapToMash.pl dbNAME
  
Example:

  perl convertMetaMapToMash.pl databases/miniSeq
  
  
);
exit;
}