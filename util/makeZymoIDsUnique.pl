use strict;
use File::Basename;
use Data::Dumper;

my $referenceDir = '/data/projects/phillippy/projects/MetaMap/loman/ZymoBIOMICS.STD.refseq.v2/Genomes/nonUniqueIDs';
my $targetDir = '/data/projects/phillippy/projects/MetaMap/loman/ZymoBIOMICS.STD.refseq.v2/Genomes/';

my %file_2_taxonID = (
	'Lactobacillus_fermentum_complete_genome.fasta' => '1613',
	'Bacillus_subtilis_complete_genome.fasta' => '1423',
	'Saccharomyces_cerevisiae_draft_genome.fa' => '4932',
	'Staphylococcus_aureus_complete_genome.fasta' => '1280',
	'Salmonella_enterica_complete_genome.fasta' => '28901',
	'Pseudomonas_aeruginosa_complete_genome.fasta' => '287',
	'Listeria_monocytogenes_complete_genome.fasta' => '1639',
	'Escherichia_coli_complete_genome.fasta' => '562',
	'Enterococcus_faecalis_complete_genome.fasta' => '1351',
	'Cryptococcus_neoformans_draft_genome.fasta' => '5207',
);

my @genomes = (glob($referenceDir . '/*.fasta'), glob($referenceDir . '/*.fa'));

my $combined_fn = $targetDir . '/combined.fa';
open(COMBINED, '>', $combined_fn) or die "Cannot open $combined_fn";

my %contigIDs;
foreach my $file (@genomes)
{
	my $basename = fileparse($file);
	my $taxonID = $file_2_taxonID{$basename};
	die unless(defined $taxonID);
	my $targetFile = $targetDir . '/' . $basename;
	open(GENOME, '<', $file) or die "Cannot open $file";
	open(SINGLE, '>', $targetFile) or die "Cannot open $targetFile";
	while(<GENOME>)
	{
		next unless($_);
		if(substr($_, 0, 1) eq '>')
		{
			my $contigID = substr($_, 1);
			chomp($contigID);
			my $newContigID = 'tx' . $taxonID . '|' . $contigID;
			# print $contigID, " -> ", $newContigID, "\n";
			die if($contigIDs{$newContigID});
			$contigIDs{$newContigID}++;
			print SINGLE '>', $newContigID, "\n";
			print COMBINED '>', $newContigID, "\n";
		}
		else
		{
			print SINGLE $_;
			print COMBINED $_;
		}
	}		
	close(SINGLE);
	close(GENOME);
}
print "\n\nDone - created genomes in $targetDir\n\nNow combine, index $combined_fn and remap!\n\n";