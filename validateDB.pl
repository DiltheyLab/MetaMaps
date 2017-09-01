use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/getcwd abs_path/;

$| = 1;

use taxTree;
use Util;
unless(scalar(@ARGV) == 1)
{
	print_help();
}

my $DB = $ARGV[0];


my %taxonID_2_contigs;
my %contigLength;
Util::read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLength);

my $taxonomyDir = $DB . '/taxonomy';
my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

my $DB_fa = $DB . '/DB.fa';

my %contigLengths_DB;
my $runningContigID;
open(DB, '<', $DB_fa) or die "Cannot open $DB_fa";
while(<DB>)
{
	my $l = $_;
	chomp($l);
	if(substr($l, 0, 1) eq '>')
	{
		$runningContigID = substr($l, 1);
		$contigLengths_DB{$runningContigID} = 0;
		my $taxonID = Util::extractTaxonID($runningContigID, $DB_fa, $.);
		die unless(exists $taxonomy->{$taxonID});
		die unless(defined $contigLength{$runningContigID});
	}
	else
	{
		$contigLengths_DB{$runningContigID} += length($l);
	}
}
close(DB);

my $n_contigs = 0;
foreach my $taxonID (keys %taxonID_2_contigs)
{
	foreach my $contigID (keys %{$taxonID_2_contigs{$taxonID}})
	{
		die unless(defined $contigLength{$contigID});
		die Dumper("Contig length discrepancy", $contigID, $contigLength{$contigID}, $contigLengths_DB{$contigID}) unless($contigLength{$contigID} == $contigLengths_DB{$contigID});
	}
}

print "\n\n$DB validated\n\n";

sub print_help
{
	print qq(
validateDB.pl

  Validate database consistency.
  
Usage:

  perl validateDB.pl dbNAME
  
Example:

  perl validateDB.pl databases/miniSeq
  
  
);
exit;
}