use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Set::IntervalTree;
use List::MoreUtils qw/mesh/;


$| = 1;

my $DB;
my $mappings;
GetOptions (
	'DB:s' => \$DB, 
	'mappings:s' => \$mappings, 
);

unless($DB and $mappings)
{
	print_help();
}

my $EM_file = $mappings . '.EM';
die "Mappings file $EM_file not existing -- have you run MetaMaps in 'classify' mode?" unless(-e $EM_file);

my $annotations_file = $DB . '/' . 'DB_annotations.txt';
die "Please specify a gene-annotated database (file $annotations_file not found)." unless(-e $annotations_file);

my $getBestMapping = sub {
	my $readID = shift;
	my $lines_aref = shift;
	
	my $bestMapQ;
	my $bestMapQ_location;
	foreach my $line (@$lines_aref)
	{
		my @line_fields = split(/ /, $line);
		die unless($line_fields[0] eq $readID);
		my $mapQ = $line_fields[13];
		die unless(defined $mapQ);
		die unless($mapQ >= 0);
		die unless($mapQ <= 1);
		
		my $contigID = $line_fields[5];
		my $contigID_start = $line_fields[7];
		my $contigID_stop = $line_fields[8];
		die unless($contigID_start <= $contigID_stop);
		my $identity = $line_fields[9] / 100;
		
		if((not defined $bestMapQ) or ($bestMapQ < $mapQ))
		{
			$bestMapQ = $mapQ;
			$bestMapQ_location = [$contigID, $contigID_start, $contigID_stop, $mapQ, $identity];	
		}
	}
	die unless(defined $bestMapQ);
	return $bestMapQ_location;
};

my %relevantContigIDs;
my $readFunc_getRelevantContigIDs = sub {
	my $readID = shift;
	my $read_lines_aref = shift;
	my $bestMapping_aref = $getBestMapping->($readID, $read_lines_aref);
	$relevantContigIDs{$bestMapping_aref->[0]}++;
};
processAllReads($EM_file, $readFunc_getRelevantContigIDs);

print "Found ", scalar(keys %relevantContigIDs),  " relevant contig IDs.\n";

my %foundAnnotations_perContig;
my %gene_2_protein_and_product;
open(ANNOTATIONS, '<', $annotations_file) or die "Cannot open $annotations_file";
my $annotations_header_line = <ANNOTATIONS>;
chomp($annotations_header_line);
my @annotations_header_fields = split(/\t/, $annotations_header_line);
die unless($annotations_header_fields[0] eq 'ContigId');
die unless($annotations_header_fields[$#annotations_header_fields - 1] eq 'CDSProduct');

while(<ANNOTATIONS>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @line_fields = split(/\t/, $line, -1);
	if(exists $relevantContigIDs{$line_fields[0]})
	{
		my %line = (mesh @annotations_header_fields, @line_fields);
		die unless(defined $line{Start});
		die unless(defined $line{Stop});
		
		unless(defined $foundAnnotations_perContig{$line{ContigId}})
		{
			$foundAnnotations_perContig{$line{ContigId}} = Set::IntervalTree->new;
		}
		
		my $geneId = $line{GeneName} . '//' . $line{GeneLocusTag};
		$foundAnnotations_perContig{$line{ContigId}}->insert($geneId, $line{Start}, $line{Stop}+1);
		die Dumper("Multi-defined gene - $line{GeneName} - $annotations_file $.") if(defined $gene_2_protein_and_product{$line{GeneName}});
		$gene_2_protein_and_product{$geneId} = [$line{GeneName}, $line{GeneLocusTag}, $line{CDSProduct}, $line{CDSProteinId}]; # CDSProteinId is actually product, and vice versa
	}
}
close(ANNOTATIONS);



print "\t... of which we have annotations for ", scalar(keys %foundAnnotations_perContig), ".\n";

my $nReads_mapped_toContigWithAnnotations = 0;
my $nReads_mapped_toContigWithoutAnnotations = 0;
my %summary_per_gene_name;
my $readFunc_getOverlappingGenes = sub {
	my $readID = shift;
	my $read_lines_aref = shift;
	my $bestMapping_aref = $getBestMapping->($readID, $read_lines_aref);
	if(exists $foundAnnotations_perContig{$bestMapping_aref->[0]})
	{
		$nReads_mapped_toContigWithAnnotations++;
		my $overlapping_genes = $foundAnnotations_perContig{$bestMapping_aref->[0]}->fetch($bestMapping_aref->[1], $bestMapping_aref->[2]);
		foreach my $gene_name (@$overlapping_genes)
		{
			$summary_per_gene_name{$gene_name}[0]++;
			push(@{$summary_per_gene_name{$gene_name}[1]}, $bestMapping_aref->[4]);
		}
	}
	else
	{
		$nReads_mapped_toContigWithoutAnnotations++;
	}
};
processAllReads($EM_file, $readFunc_getOverlappingGenes);


print "Of ", ($nReads_mapped_toContigWithAnnotations + $nReads_mapped_toContigWithoutAnnotations), " mapped reads, $nReads_mapped_toContigWithAnnotations go to contigs with annotations and $nReads_mapped_toContigWithoutAnnotations to contigs without.\n";

my $output_file = $EM_file . '.geneLevelAnalysis';
open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
print OUTPUT join("\t", "GeneName", "GeneLocusTag", "ProteinId", "Product", "nReads", "medianIdentity"), "\n";
foreach my $geneName (keys %summary_per_gene_name)
{
	my @identities = @{$summary_per_gene_name{$geneName}[1]};
	my $medianIdentity = getMedian(\@identities);
	print OUTPUT join("\t", $gene_2_protein_and_product{$geneName}[0], $gene_2_protein_and_product{$geneName}[1], $gene_2_protein_and_product{$geneName}[2], $gene_2_protein_and_product{$geneName}[3], $summary_per_gene_name{$geneName}[0], $medianIdentity), "\n";
	
}
close(OUTPUT);
print "\n\nDone. Produced $output_file\n\n";

sub getMedian
{
	my $aref = shift;
	my @l = sort @$aref;
	if(scalar(@l) == 1)
	{
		return $l[0];
	}
	else
	{
		my $medianI = int(scalar(@l) / 2 + 0.5) - 1;
		die unless($medianI >= 0);
		die unless($medianI <= $#l);
		return $l[$medianI];
	}
}
sub processAllReads
{
	my $EM_file = shift;
	my $F = shift;
	
	my @runningLines;
	my $runningReadID = '';
	open(EM, '<', $EM_file) or die "Cannot open $EM_file";
	while(<EM>)
	{
		chomp;
		next unless($_);
		my $line = $_;
		die "Weird line $. in $EM_file:\n$line" unless($line =~ /^(\S+) /);
		my $readID = $1;
		if($runningReadID ne $readID)
		{
			if(@runningLines)
			{
				$F->($runningReadID, \@runningLines);
			}	
			$runningReadID = $readID;
			@runningLines = ();
		}
		
		push(@runningLines, $line);
	}
	close(EM);
	if(@runningLines)
	{
		$F->($runningReadID, \@runningLines);
	}
}

sub print_help
{
	print qq(
geneLevelAnalysis.pl

  Varry out a gene-level analysis. Only works for databases with
  annotation information (i.e. with a file DB_annotations.txt in the DB dir).
  
Usage:

  perl geneLevelAnalysis.pl --DB dbPATH --mappings PREFIX
  
Parameters:

  DB
      Path to database directory.
  
  mappings
  
      Path to the mappings file. Requires that these have
	  been analyzed with MetaMaps in 'classify' mode. 
);
exit;
}