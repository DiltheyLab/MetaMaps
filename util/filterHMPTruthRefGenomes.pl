use strict;
use warnings;
use Data::Dumper;

my %sawTaxa;
my %taxa_genome_lengths;
open(GL, '<', '../tmp/HMP_genome_lengths.txt') or die "Cannot open ../tmp/HMP_genome_lengths.txt";
while(<GL>)
{
	my $l = $_;
	chomp($l);
	next unless($l);
	my @f = split(/\t/, $l);
	die unless(scalar(@f) == 2);
	$taxa_genome_lengths{$f[0]} = $f[1];
}
close(GL);


my %gi_2_taxon;
open(GI2TAXON, '<', '/data/projects/phillippy/projects/mashsim/db/gi_taxid_nucl.dmp.HMP') or die;
while(<GI2TAXON>)
{
	my $line = $_; 
	chomp($line);
	my @f = split(/\s+/, $line);
	die unless($#f == 1);	
	$gi_2_taxon{$f[0]} = $f[1];
}
close(GI2TAXON);

$gi_2_taxon{126640115} = '400667';
$gi_2_taxon{126640097} = '400667';
$gi_2_taxon{126640109} = '400667';
$gi_2_taxon{161510924} = '451516';
$gi_2_taxon{32470532} = '176280';

foreach my $genome ('/data/projects/phillippy/projects/mash_map/Jobs/blasr/hmp/targetAll/mock.all.genome.fa')
{
	my $genome_href = readFASTA($genome);
	my $genome_filtered_href = {};
	foreach my $contigID (keys %$genome_href)
	{
		die "Invalid contig ID - no GI! $contigID" unless($contigID =~ /gi\|(\d+)\|/);
		my $gi = $1;
		next unless(exists $gi_2_taxon{$gi});
		my $taxonID = $gi_2_taxon{$gi};
		if(exists $taxa_genome_lengths{$taxonID})
		{
			$genome_filtered_href->{$contigID} = $genome_href->{$contigID};
			$sawTaxa{$taxonID}++;
		}
					
	}
	writeFASTA($genome . '.filteredToHMP', $genome_filtered_href);
	
	print "Generated file: ", $genome . '.filteredToHMP', "\n";
}

foreach my $taxonID (keys %taxa_genome_lengths)
{
	unless($sawTaxa{$taxonID})
	{
		warn "Taxon $taxonID not in genomes file";
	}
}


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

sub writeFASTA
{
	my $file = shift;
	# print "Writing $file\n";
	my $href = shift;
	open(F, '>', $file) or die "Cannot open $file";
	foreach my $key (sort keys %$href)
	{
		my $seq = $href->{$key};
		print F '>', $key, "\n";
		# print "\t", $key, "\t", length($seq), "\n";
		while($seq)
		{
			my $toPrint;
			if(length($seq) > 50)
			{
				$toPrint = substr($seq, 0, 50);
				substr($seq, 0, 50) = '';
			}
			else
			{
				$toPrint = $seq;
				$seq = '';
			}	
			print F $toPrint, "\n";
		}
	}
	close(F);	
}



