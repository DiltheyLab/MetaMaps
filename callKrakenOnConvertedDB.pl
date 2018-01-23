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

unless(scalar(@ARGV) == 3)
{
	print_help();
}

my $database = $ARGV[0];
my $FASTQ = $ARGV[1];
my $outputDir = $ARGV[2];

# test that database is valid
# todo reinstate
#my %taxonID_2_contigs;
#my %contigLength;
#Util::read_taxonIDs_and_contigs($database, \%taxonID_2_contigs, \%contigLength);

unless(-e $outputDir)
{
	mkdir($outputDir) or die "Cannot mkdir $outputDir";
}


my %taxonID_original_2_contigs;
my %contigLength;
Util::read_taxonIDs_and_contigs($database, \%taxonID_original_2_contigs, \%contigLength);

	
my $krakenDir = $database . '/kraken';

SimulationsKraken::doKrakenOnExistingDB (
	$krakenDir,
	$FASTQ,
	$outputDir,
	$kraken_binPrefix,
	$Bracken_dir,
	\%taxonID_original_2_contigs
);

sub print_help
{
	print qq(
callKrakenOnConvertedDB.pl

  Call Kraken/Bracken with a FASTQ file on a converted DB.
  
Usage:

  perl callKrakenOnConvertedDB.pl dbNAME FASTQFILE outputDir
  
Example:

  perl callKrakenOnConvertedDB.pl databases/miniSeq test.fastq testKrakenResults
  
  
);
exit;
}