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

my $centrifugeBinDir = SimulationsKraken::getCentrifugeDir();


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

	
my $centrifugeDBDir = $database . '/centrifuge';

SimulationsKraken::doCentrifugeOnExistingDB (
	$centrifugeDBDir,
	$FASTQ,
	$outputDir,
	$centrifugeBinDir,
	\%taxonID_original_2_contigs,
	$database
);

sub print_help
{
	print qq(
callCentrifugeOnConvertedDB.pl

  Call Centrifuge with a FASTQ file on a converted DB.
  
Usage:

  perl callCentrifugeOnConvertedDB.pl dbNAME FASTQFILE outputDir
  
Example:

  perl callCentrifugeOnConvertedDB.pl databases/miniSeq test.fastq testCentrifugeResults
  
  
);
exit;
}