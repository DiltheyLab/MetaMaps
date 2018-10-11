use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../perlLib"; 
use File::Basename;

$| = 1;

use Util;
use taxTree;
use simulation;
use validation;
	
my $prefix_out = '../tmp/truthCAMI';

my $masterTaxonomy_dir = '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';
my $master_taxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);



my $fn_reads = '/data/projects/phillippy/software/camiClient/19122017_mousegut_pacbio_scaffolds/2018.02.13_14.02.01_sample_0/reads/anonymous_reads.fq';
my $fn_reads_truth = '/data/projects/phillippy/software/camiClient/19122017_mousegut_pacbio_scaffolds/2018.02.13_14.02.01_sample_0/reads/reads_mapping.tsv';

my $fn_out_reads = $prefix_out . '.perRead';
my $fn_out_distribution = $prefix_out . '.distribution';
my $fn_out_distribution_genomeFreqs = $prefix_out . '.distribution_genomes';
my $fn_out_origin = $prefix_out . '.genomesOfOrigin';

my %readID_2_length;
open(READS, '<', $fn_reads) or die "Cannot open $fn_reads";
while(<READS>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	die Dumper("Weird line $line", $., $fn_reads) unless(substr($line, 0, 1) eq '@');
	my $readID = substr($line, 1);
	$readID =~ s/\s.+//; 
	my $seq = <READS>;
	my $plus = <READS>;
	my $qual = <READS>;
	$readID_2_length{$readID} = length($seq);		
}
close(READS);

my %read_2_taxonID;

my %origin_genome_ids;
my %origin_to_taxonID;
my %taxonID_to_origins;
open(TRUTH, '<', $fn_reads_truth) or die "Cannot open $fn_reads_truth";
my $truth_header = <TRUTH>;
chomp($truth_header);
my @header_fields = split(/\t/, $truth_header);
die unless($header_fields[0] eq '#anonymous_read_id');
die unless($header_fields[1] eq 'genome_id');
die unless($header_fields[2] eq 'tax_id');
die unless($header_fields[3] eq 'read_id');
while(<TRUTH>)
{
	chomp;
	my @line_fields = split(/\t/, $_);
	my $readID = $line_fields[0];
	my $taxonID = $line_fields[2];
	my $origin = $line_fields[3];
	die "Weird origin $origin" unless($origin =~ /^(\w+\.\d)/);
	$origin = $1;
	$origin_genome_ids{$origin}++;
	die "Taxon ID $taxonID not part of master taxonomy" unless($master_taxonomy->{$taxonID});
	die unless(defined $readID_2_length{$readID});
	$read_2_taxonID{$readID} = $taxonID;
	
	die unless((not defined $origin_to_taxonID{$origin}) or ($origin_to_taxonID{$origin} eq $taxonID));
	$origin_to_taxonID{$origin} = $taxonID;
	$taxonID_to_origins{$taxonID}{$origin}++;
}
close(TRUTH);

foreach my $readID (keys %readID_2_length)
{
	die unless(defined $read_2_taxonID{$readID});
}



open(ORIGIN, '>', $fn_out_origin) or die "Cannot open $fn_out_origin";
print ORIGIN join("\n", keys %origin_genome_ids), "\n";
close(ORIGIN);

my $fn_in_origins_fasta = $fn_out_distribution . '.fasta';
my $genomes_href = readFASTA($fn_in_origins_fasta);

my %taxa_genome_lengths;
foreach my $taxonID (keys %taxonID_to_origins)
{
	foreach my $origin (keys %{$taxonID_to_origins{$taxonID}})
	{
		die "Origin seqeunce $origin not in file $fn_in_origins_fasta" unless($genomes_href->{$origin});
		$taxa_genome_lengths{$taxonID} += length($genomes_href->{$origin});
	}
}

my %taxonID_read_counts;
my %taxonID_2_bases;
open(OUT_PERREAD, '>', $fn_out_reads) or die "Cannot open file $fn_out_reads";
foreach my $readID (keys %read_2_taxonID)
{
	print OUT_PERREAD join("\t", $readID, $read_2_taxonID{$readID}), "\n";
	$taxonID_read_counts{$read_2_taxonID{$readID}}++;
	my $length = $readID_2_length{$readID};
	die unless(defined $length);
	$taxonID_2_bases{$read_2_taxonID{$readID}} += $length;
}	 
close(OUT_PERREAD);



simulation::truthReadFrequenciesFromReadCounts($fn_out_distribution, \%taxonID_read_counts, $master_taxonomy);
simulation::truthGenomeFrequenciesFromReadCounts($fn_out_distribution_genomeFreqs, \%taxonID_2_bases, \%taxonID_read_counts, \%taxa_genome_lengths, $master_taxonomy);

sub readFASTA
{
	my $file = shift;	
	my $cut_sequence_ID_after_whitespace = shift;
	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			if($cut_sequence_ID_after_whitespace)
			{
				$line =~ s/\s+.+//;
			}
			$currentSequence = substr($line, 1);
			$R{$currentSequence} = '';
		}
		else
		{
			die "Weird input in $file" unless (defined $currentSequence);
			$R{$currentSequence} .= uc($line);
		}
	}	
	close(F);
		
	return \%R;
}



