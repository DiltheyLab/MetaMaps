use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../perlLib";
use taxTree;

$| = 1;

my $DB;
my $r2t;
my $FASTQ;
my $target;

GetOptions (
	'DB:s' => \$DB, 
	'r2t:s' => \$r2t, 
	'FASTQ:s' => \$FASTQ,
	'target:s' => \$target,
);

unless($DB and $r2t and $FASTQ and (defined $target))
{
	print_help();
}

unless(-e $r2t)
{
	die "Reads-to-taxon file $r2t not existing";
}

my $taxonomyDir = $DB . '/taxonomy';
my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

unless(($target eq '0') or (exists $taxonomy->{$target}))
{
	die "Specified --target $target is no known node of the $DB database taxonomy";
}

my @_targetNodes = ($target, ($target ne '0') ? taxTree::descendants($taxonomy, $target) : ());
my %targetNodes = map {$_ => 1} @_targetNodes;

print STDERR "\nHave ", scalar(keys %targetNodes), " target nodes.\n";

my %targetReads;
open(R2T, '<', $r2t) or die "Cannot open --r2t $r2t";
while(<R2T>)
{
	chomp;
	next unless($_);
	my @fields = split(/\t/, $_);
	die "Line $. in $r2t doesn't have 2 fields" unless(scalar(@fields) == 2);
	die "Taxonomy mismatch - ID $fields[1] neither 0 nor part of DB taxonomy" unless(($fields[1] eq '0') or (exists $taxonomy->{$fields[1]}));
	if($targetNodes{$fields[1]})
	{
		$targetReads{$fields[0]}++;
	}
}
close(R2T);

print STDERR "Targeting ", scalar(keys %targetReads), " reads for extraction.\n";

my %printedReads;
open(FASTQ, '<', $FASTQ) or die "Cannot open $FASTQ";
while(<FASTQ>)
{
	chomp;
	next unless($_);
	die "Corrupted FASTQ: Line $. of $FASTQ should be a read ID line, but isn't" unless(substr($_, 0, 1) eq '@');
	my $readID_line = $_;
	my $readID = substr($_, 1);
	$readID =~ s/\s.*//;
	my $sequence = <FASTQ>;
	my $plus = <FASTQ>; die unless(substr($plus, 0, 1) eq '+');
	my $qualities = <FASTQ>;
	die unless(length($sequence) == length($qualities));
	if($targetReads{$readID})
	{
		print $readID_line, "\n", $sequence, $plus, $qualities;
		$printedReads{$readID}++;
	}
}
close(FASTQ);

print STDERR "Done. Got ", scalar(keys %printedReads), " of ", scalar(keys %targetReads), " reads.\n";

sub print_help
{
	print qq(
extractReads.pl

  Extract reads mapping to a particular node of the taxonomy
  (and its ancestors).
  
Usage:

  perl extractReads.pl --DB dbDIR --r2t R2TFILE --FASTQ FASTQFILE --target taxonomyID|0
  
Parameters:

  DB
      Path to database directory.
  
  r2t
  
      Read-to-taxon file (e.g. /path/to/analysis.EM.reads2Taxon).

  FASTQ

      Path to the original FASTQ file.

   target

      Target node in the taxonomy for extraction
      (or '0' for reads that are long enough but unmapped).	    
);
exit;
}
