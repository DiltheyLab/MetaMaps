use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;
use FindBin;
use lib "$FindBin::Bin/../perlLib";
use Util;

$| = 1;

my $DB = 'databases/miniSeq+H';
my @files = ('/data/projects/phillippy/projects/MetaMap/databases/miniSeq+H/simulations_p25_logNormal/0/truth_genomeFrequencies.txt', '/data/projects/phillippy/projects/MetaMap/databases/miniSeq+H/simulations_i100_specifiedFrequencies/0/truth_genomeFrequencies.txt');

my %taxonID_2_contigs;
my %contigLength;
Util::read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLength);


my $taxonomyDir = $DB . '/taxonomy';
my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

foreach my $file (@files)
{
	my $f_out = $file . '.withNC';
	open(F, '<', $file) or die "Cannot open $file";
	open(FOUT, '>', $f_out) or die "Cannot open $f_out for writing";
	my $headerLine = <F>;
	chomp($headerLine);
	my @headerFields = split(/\t/, $headerLine);
	print FOUT join("\t", @headerFields, 'NCs'), "\n";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @lineFields = split(/\t/, $line);
		die unless($#lineFields == $#headerFields);
		my %line = (mesh @headerFields, @lineFields);
		my $taxonID = $line{taxonID};
		die Dumper("Unknown taxon ID in $file: $taxonID", \%line) unless($taxonomy->{$taxonID});
		die Dumper("No contigs for $taxonID", \%line) unless ($taxonID_2_contigs{$taxonID});
		my @contigs = keys %{$taxonID_2_contigs{$taxonID}};
		my %NCs = map {$_ => 1} map {die unless($_ =~ /\|(NC_.+)/); $1} @contigs;
		my @NCs = (keys %NCs);
		print FOUT join("\t", @lineFields, join(',', @NCs)), "\n";
	}
	close(F);
	close(FOUT);
	print "Produced $f_out\n";
}	