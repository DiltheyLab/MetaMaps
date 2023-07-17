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

my $kraken2_binPrefix = SimulationsKraken::getKraken2BinPrefix();

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
$outputDir = abs_path($outputDir);


my %taxonID_original_2_contigs;
my %contigLength;
Util::read_taxonIDs_and_contigs($database, \%taxonID_original_2_contigs, \%contigLength);

	
my $kraken2Dir = $database . '/kraken2';

SimulationsKraken::doKraken2OnExistingDB (
	$kraken2Dir,
	$FASTQ,
	$outputDir,
	$kraken2_binPrefix,
	\%taxonID_original_2_contigs
);

sub print_help
{
	print qq(
callKraken2OnConvertedDB.pl

  Call Kraken2 with a FASTQ file on a converted DB.
  
Usage:

  perl callKraken2OnConvertedDB.pl dbNAME FASTQFILE outputDir
  
Example:

  perl callKraken2OnConvertedDB.pl databases/miniSeq test.fastq testKrakenResults
  
  
);
exit;
}