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

my $centrifugeDir = SimulationsKraken::getCentrifugeDir();

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

my $outputDir_centrifuge = $database . '/centrifuge';

SimulationsKraken::translateMetaMapToCentrifuge (
	$outputDir_centrifuge,
	$database, 
	$centrifugeDir,
);


sub print_help
{
	print qq(
convertMetaMapsToKraken.pl

  Convert a MetaMaps DB to Kraken/Bracken/Kraken2.
  
Usage:

  perl convertMetaMapsToKraken.pl dbNAME
  
Example:

  perl convertMetaMapsToKraken.pl databases/miniSeq
  
  
);
exit;
}