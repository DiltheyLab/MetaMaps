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

my $kraken_binPrefix = SimulationsKraken::getKrakenBinPrefix();
my $Bracken_dir = SimulationsKraken::getBrackenDir();
my $krakenDBTemplate = SimulationsKraken::getKrakenDBTemplate();

unless(scalar(@ARGV) == 1)
{
	print_help();
}

my $database = $ARGV[0];

# test that database is valid
# todo reinstate
#my %taxonID_2_contigs;
#my %contigLength;
#Util::read_taxonIDs_and_contigs($database, \%taxonID_2_contigs, \%contigLength);

my $outputDir_kraken = $database . '/kraken';

SimulationsKraken::translateMetaMapToKraken (
	$outputDir_kraken,
	$database,
	$krakenDBTemplate,
	$kraken_binPrefix,
	$Bracken_dir
);

sub print_help
{
	print qq(
convertMetaMapToKraken.pl

  Convert a MetaMap DB to Kraken/Bracken.
  
Usage:

  perl convertMetaMapToKraken.pl dbNAME
  
Example:

  perl convertMetaMapToKraken.pl databases/miniSeq
  
  
);
exit;
}