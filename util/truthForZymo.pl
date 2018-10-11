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

my $prefix_out = '../tmp/truthZymp';
# my $DB = '../databases/miniSeq+H';

	
# my $metaMap_taxonomy_dir = $DB . '/taxonomy';
# my $MetaMap_taxonomy = taxTree::readTaxonomy($metaMap_taxonomy_dir);

# my %taxonID_2_contigs;
# my %contigLengths;
# Util::read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLengths);

my $masterTaxonomy_dir = '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';
my $master_taxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);



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

my %contig_from_genome;
my %taxa_genome_lengths;
my %contig_lengths;
my $referenceDir = '/data/projects/phillippy/projects/MetaMap/loman/ZymoBIOMICS.STD.refseq.v2/Genomes/';
my @genomes = (glob($referenceDir . '/*.fasta'), glob($referenceDir . '/*.fa'));
foreach my $genome (@genomes)
{
	next if($genome =~ /combined/);
	my $genome_basename = fileparse($genome);
	die "Genome $genome_basename undefined" unless($file_2_taxonID{$genome_basename});
	my $currentContigID;
	open(GENOME, '<', $genome) or die;
	while(<GENOME>)
	{	
		if(substr($_, 0, 1) eq '>')
		{
			chomp($_);
			substr($_, 0, 1) = '';
			my $contigID = $_;
			$contigID =~ s/\s.+//;			
			die if($contig_from_genome{$contigID});
			$contig_from_genome{$contigID} = $genome_basename;
			$currentContigID = $contigID;
		}		
		else
		{
			chomp($_);
			die unless(defined $currentContigID);
			$contig_lengths{$currentContigID} += length($_);
				
			die unless($file_2_taxonID{$genome_basename});
			$taxa_genome_lengths{$file_2_taxonID{$genome_basename}} += length($_);
		}	
	}
	close(GENOME);
}

my $referenceGenome = $referenceDir . '/combined.fasta';
my %reference_contigs_2_taxon;
my %reference_contigs_unique;
open(REF, '<', $referenceGenome) or die;
while(<REF>)
{
	next unless($_);
	if(substr($_, 0, 1) eq '>')
	{
		chomp($_);
		substr($_, 0, 1) = '';
		my $contigID = $_;
		$contigID =~ s/\s.+//;
		die unless($contig_from_genome{$contigID});
		die unless($file_2_taxonID{$contig_from_genome{$contigID}});
		my $taxonID = $file_2_taxonID{$contig_from_genome{$contigID}};
		die unless(defined $master_taxonomy->{$taxonID});
		$reference_contigs_2_taxon{$contigID} = $taxonID;
	}	
}
close(REF);

foreach my $config (['bwa_nanopore', '/data/projects/phillippy/projects/MetaMap/loman/ZymoBIOMICS.STD.refseq.v2/Genomes/GridION-Zymo_CS_MPZBB_LSK109.all.fq.bam', '/data/projects/phillippy/projects/MetaMap/loman/GridION-Zymo_CS_MPZBB_LSK109.all.fq'])
{
	print $config->[0], "\n";
	 
	my $Zymo_fastQ = $config->[2];
	my $Zymo_readIDs_href = getReadIDs($Zymo_fastQ);
	
	my $fn_out_reads = $prefix_out . '_' . $config->[0] . '.perRead';
	my $fn_out_distribution = $prefix_out . '_' . $config->[0] . '.distribution';
	my $fn_out_distribution_genomeFreqs = $prefix_out . '_' . $config->[0] . '.distribution_genomes';

	my %alignments_per_readID;
	my %alignments_per_longReadID;
	my %NC_2_gi;
	my %haveReadInFastQ;
	my %noReadInFastQ;
	my %readID_2_length;
	my %read_2_taxonID;

	die unless($config->[1] =~ /\.bam/);
	{
		my $n_reads = 0;
					
		open(BAM, '-|', "samtools view -F 0x800 -F 0x100 -F 0x4 $config->[1]") or die "Cannot pipe-open $config->[1]";

		while(<BAM>)
		{
			my $line = $_;
			chomp($line);
			my @fields = split(/\s+/, $line);
			my $longReadID = $fields[0];
			if(exists $Zymo_readIDs_href->{$longReadID})
			{
				#next;
				#die "Read ID $readID not in HMP FASTQ $Zymo_fastQ";
				die if($haveReadInFastQ{$longReadID});
				$haveReadInFastQ{$longReadID}++;
			}
			else
			{
				die Dumper("Read ID not in FASTQ?", $longReadID);
				$noReadInFastQ{$longReadID}++;
				next;
			}

			my $contigID = $fields[2];
			my $mapQ = $fields[4];
			my $seq = $fields[9];
			
			die unless($mapQ =~ /^\d+$/);
			
			$n_reads++;
			$alignments_per_readID{$longReadID}++;
			$alignments_per_longReadID{$longReadID}++;
		
			if(exists $readID_2_length{$longReadID})
			{
				die unless($readID_2_length{$longReadID} = length($seq));
			}
			
			$readID_2_length{$longReadID} = length($seq);		

			die "Undefined contig ID $contigID" unless(defined $reference_contigs_2_taxon{$contigID});
			$read_2_taxonID{$longReadID} = $reference_contigs_2_taxon{$contigID};
			
		}
		close(BAM);	
	}
	
	print "Reads - no alignments   - also in FASTQ:", scalar(grep {not exists $alignments_per_longReadID{$_}} keys %$Zymo_readIDs_href), "\n";
	print "Reads - with alignments - also in FASTQ:", scalar(keys %haveReadInFastQ), "\n";
	print "Reads - with alignments - not  in FASTQ: ", scalar(keys %noReadInFastQ), "\n";
	
	# statistics

	my %histogram_n_alignments;
	foreach my $readID (keys %alignments_per_readID)
	{
		my $n_alignments = $alignments_per_readID{$readID};
		die "This no primary?" unless($n_alignments == 1);
		$histogram_n_alignments{$n_alignments}++;
	}

	print "Number of reads: ", scalar(keys %alignments_per_readID), "\n";
	print "Number-of-alignments histogram:\n";
	foreach my $n_alignment (sort keys %histogram_n_alignments)
	{
		next if($histogram_n_alignments{$n_alignment} < 100);
		print "\t", $n_alignment, "\t", $histogram_n_alignments{$n_alignment}, "\n";
	}
	
	open(OUT_PERREAD, '>', $fn_out_reads) or die "Cannot open file $fn_out_reads";
	my %taxonID_read_counts;
	my %taxonID_2_bases;
	foreach my $readID (keys %read_2_taxonID)
	{
		print OUT_PERREAD join("\t", $readID, $read_2_taxonID{$readID}), "\n";
		$taxonID_read_counts{$read_2_taxonID{$readID}}++;
		my $length = $readID_2_length{$readID};
		die unless(defined $length);
		
		$taxonID_2_bases{$read_2_taxonID{$readID}} += $length;
	}	 

	foreach my $readID (keys %$Zymo_readIDs_href)
	{
		next if(defined $read_2_taxonID{$readID});
		$taxonID_read_counts{0}++;
	}

	close(OUT_PERREAD);

	foreach my $taxonID (keys %taxonID_2_bases)
	{
		unless(exists $taxa_genome_lengths{$taxonID})
		{
			die "Missing length information for taxon $taxonID";
		}
	}
	
	simulation::truthReadFrequenciesFromReadCounts($fn_out_distribution, \%taxonID_read_counts, $master_taxonomy);
	simulation::truthGenomeFrequenciesFromReadCounts($fn_out_distribution_genomeFreqs, \%taxonID_2_bases, \%taxonID_read_counts, \%taxa_genome_lengths, $master_taxonomy);

	my $fn_out_FASTQ = $Zymo_fastQ . '.mappable';
	open(FASTQ_OUT, '>', $fn_out_FASTQ) or die "Cannot open $fn_out_FASTQ";
	open(FASTQ_IN, '<', $Zymo_fastQ) or die "Cannot open $Zymo_fastQ";
	while(<FASTQ_IN>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		die Dumper("Weird line $line", $., $Zymo_fastQ) unless(substr($line, 0, 1) eq '@');
		my $readID = substr($line, 1);
		$readID =~ s/\s.+//; 
		my $seq = <FASTQ_IN>;
		my $plus = <FASTQ_IN>;
		my $qual = <FASTQ_IN>;
		if(defined $read_2_taxonID{$readID})
		{
			die unless($read_2_taxonID{$readID});
			print FASTQ_OUT $line, "\n", $seq, $plus, $qual;
		}
	}
	close(FASTQ_IN);
	close(FASTQ_OUT);
	print "\n\nDone $config->[0] . Produced files:\n";
	print "\t - $fn_out_reads \n";
	print "\t - $fn_out_distribution \n";
	print "\t - $fn_out_distribution_genomeFreqs \n";
	print "\t - $fn_out_FASTQ \n"; 
	print "\n";
}

sub getReadIDs
{
	my $fn = shift;
	
	my %forReturn;
	open(F, '<', $fn) or die;
	my $isFirstLine = 1;
	my $FASTA = 0;
	while(<F>)
	{
		chomp;
		next unless($_);
		if($isFirstLine)
		{
			$FASTA = 1 if(substr($_, 0, 1) eq '>');
			$isFirstLine = 0;
		}
		
		if($FASTA)
		{
			my $readID = $_;
			die unless(substr($readID, 0, 1) eq '>');
			substr($readID, 0, 1) = '';
			<F>;
			$readID =~ s/\s.+//; 
			die if($forReturn{$readID});
			$forReturn{$readID}++;		
		}
		else
		{
			my $readID = $_;
			die unless(substr($readID, 0, 1) eq '@');
			substr($readID, 0, 1) = '';
			$readID =~ s/\s.+//; 
			<F>;
			my $plus = <F>;
			die unless(substr($plus, 0, 1) eq '+');
			<F>;
			$forReturn{$readID}++;
		}
	}
	close(F);
	return \%forReturn;
}