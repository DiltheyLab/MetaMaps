use strict;
use warnings;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../perlLib"; 
$| = 1;

use Util;
use taxTree;
use simulation;
use validation;

my $prefix_out = '../tmp/truthHMP7';
my $targetDB = '../databases/miniSeq';

my $DB = '../databases/miniSeq+H';

	
my $metaMap_taxonomy_dir = $DB . '/taxonomy';
my $MetaMap_taxonomy = taxTree::readTaxonomy($metaMap_taxonomy_dir);




	
my %taxonID_2_contigs;
my %contigLengths;
Util::read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLengths);

my %NC_2_taxon;
my %NC_noPoint_2_taxon;
foreach my $taxonID (keys %taxonID_2_contigs)
{
	foreach my $contigID (keys %{$taxonID_2_contigs{$taxonID}})
	{
		my @contigID_parts = split(/\|/, $contigID);
		my $NC = $contigID_parts[$#contigID_parts];
		if(defined $NC_2_taxon{$NC})
		{
			die unless($NC_2_taxon{$NC} eq $taxonID);
		}
		$NC_2_taxon{$NC} = $taxonID;
		(my $NC_noPoint = $NC) =~ s/\.\d+$//;
		if(defined $NC_noPoint_2_taxon{$NC_noPoint})
		{
			die unless($NC_noPoint_2_taxon{$NC_noPoint} eq $taxonID);
		}		
		$NC_noPoint_2_taxon{$NC_noPoint} = $taxonID;
		
	}
}

my $masterTaxonomy_dir = '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';
(my $master_taxonomy, my $master_taxonomy_merged) = validation::prepare_masterTaxonomy_withX($masterTaxonomy_dir, $MetaMap_taxonomy);

# my $master_taxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
# my $master_taxonomy_merged = taxTree::readMerged($masterTaxonomy_dir);

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


# foreach my $config (['bwa_nanopore', '/data/projects/phillippy/projects/mash_map/Jobs/blasr/hmp/targetAll/mock.all.genome.fa.nanopore.sorted.bam', '/scratch/tmp/hmp-nanopore.fasta'], ['blasr_pacbio', '/data/projects/phillippy/projects/mash_map/Jobs/blasr/hmp/target/all.m4', '/scratch/tmp/hmp_set7_combined.fastq'], ['bwa_pacbio', '/data/projects/phillippy/projects/mash_map/Jobs/blasr/hmp/targetAll/mock.all.genome.fa.pacbioReads.bam', '/scratch/tmp/hmp_set7_combined.fastq'])
foreach my $config (['bwa_pacbio', '/data/projects/phillippy/projects/mash_map/Jobs/blasr/hmp/targetAll/mock.all.genome.fa.pacbioReads.bam', '/scratch/tmp/hmp_set7_combined.fastq'], ['bwa_nanopore', '/data/projects/phillippy/projects/mash_map/Jobs/blasr/hmp/targetAll/mock.all.genome.fa.nanopore.sorted.bam', '/scratch/tmp/hmp-nanopore.fasta'])
{
	print $config->[0], "\n";
	
	my $HMP_fastQ = $config->[2];
	my $HMP_readIDs_href = getReadIDs($HMP_fastQ);

	my $fn_out_reads = $prefix_out . '_' . $config->[0] . '.perRead';
	my $fn_out_distribution = $prefix_out . '_' . $config->[0] . '.distribution';
	my $fn_out_distribution_genomeFreqs = $prefix_out . '_' . $config->[0] . '.distribution_genomes';

	my %alignments_per_readID;
	my %alignments_per_longReadID;
	my %read_2_gis;
	my %read_2_NCs;
	my %NC_2_gi;
	my %haveReadInFastQ;
	my %noReadInFastQ;
	my %readID_2_length;
	if($config->[1] =~ /\.m4/)
	{
		my $n_read_blasr = 0;
		open(BLASRTRUTH, '<', $config->[1]) or die "Cannot open $config->[1]";

		while(<BLASRTRUTH>)
		{
			my $line = $_;
			chomp($line);
			my @fields = split(/\s+/, $line);
			my $longReadID = $fields[0];
			die "Can't parse read ID $longReadID from $config->[1]" unless($longReadID =~ /(^.+)\/\d+_\d+$/);
			my $readID = $1;
			if(exists $HMP_readIDs_href->{$readID})
			{
				#next;
				#die "Read ID $readID not in HMP FASTQ $HMP_fastQ";
				$haveReadInFastQ{$readID}++;
			}
			else
			{
				$noReadInFastQ{$readID}++;
				next;
			}
			# print join("\t", $longReadID, $readID), "\n";
			my $contigID = $fields[1];
			my $identity = $fields[3];
			die unless($identity >= 2); die unless($identity <= 100);
			
			my $alignment_read_start = $fields[5];
			my $alignment_read_stop = $fields[6];
			die unless($alignment_read_start < $alignment_read_stop);
			my $alignment_read_length = $alignment_read_stop - $alignment_read_start + 1;
			
			my $read_length = $fields[7];
			die Dumper($read_length, $alignment_read_stop) unless($read_length >= $alignment_read_stop);
			
			my $alignment_cover = $alignment_read_length/$read_length;
			#next unless($alignment_cover >= 0.7);
			$alignments_per_readID{$readID}++;
			$alignments_per_longReadID{$longReadID}++;
			
			die "Invalid contig ID - no GI! $contigID" unless($contigID =~ /gi\|(\d+)\|/);
			my $gi = $1;
			#next if($gi eq '148642060');
			
			push(@{$read_2_gis{$readID}}, [$gi, $alignment_read_length * ($identity/100)]);
			if(exists $readID_2_length{$readID})
			{
				die unless($readID_2_length{$readID} == $read_length);
			}
			$readID_2_length{$readID} = $read_length;
			
			$n_read_blasr++;

		}
		close(BLASRTRUTH);
	}
	elsif($config->[1] =~ /\.bam/)
	{
		my $n_reads = 0;
					
		open(BAM, '-|', "samtools view -F 0x800 -F 0x100 -F 0x4 $config->[1]") or die "Cannot pipe-open $config->[1]";

		while(<BAM>)
		{
			my $line = $_;
			chomp($line);
			my @fields = split(/\s+/, $line);
			my $longReadID = $fields[0];
			
			my $readID;
			if($config->[1] =~ /nanopore/i)
			{
				$readID = $longReadID;

			}
			else
			{
				die "Can't parse read ID $longReadID from $config->[1]" unless($longReadID =~ /(^.+)\/\d+_\d+$/);
				$readID = $1;	
				$readID = $longReadID;
			}

			if(exists $HMP_readIDs_href->{$longReadID})
			{
				#next;
				#die "Read ID $readID not in HMP FASTQ $HMP_fastQ";
				die if($haveReadInFastQ{$longReadID});
				$haveReadInFastQ{$longReadID}++;
			}
			else
			{
				die Dumper("Read ID not in FASTQ?", $longReadID, $readID);
				$noReadInFastQ{$longReadID}++;
				next;
			}
			# print join("\t", $longReadID, $readID), "\n";
			my $contigID = $fields[2];
			my $mapQ = $fields[4];
			my $seq = $fields[9];
			
			die unless($mapQ =~ /^\d+$/);
			die "Invalid contig ID - no GI! $contigID" unless($contigID =~ /gi\|(\d+)\|/);
			my $gi = $1;

			die "Invalid contig ID" unless(substr($contigID, length($contigID) - 1, 1) eq '|');
			my @contigID_parts = split(/\|/, substr($contigID, 0, length($contigID) - 1));
			my $NC = $contigID_parts[$#contigID_parts];
			# die Dumper($contigID, $NC);
			#die "Invalid contig ID - no NC! $contigID" unless($contigID =~ /(((NC)|(AE)|(DS)).+?)\|/);
			#my $NC = $1;
					
			#next if($gi eq '148642060');
					
			die "Duplicate short read ID $readID from $longReadID in file $config->[1]" if(exists $read_2_gis{$readID});
			
			# these are ID substitutions - we checked that the represented sequences are identical
			if($NC eq 'AE017194.1')
			{
				$NC = 'NC_003909.8';
			}
			if($NC eq 'ABPI01000001.1')
			{
				$NC = 'NC_017316.1';
			}
			
			if(defined $NC_2_gi{$NC})
			{
				die unless($NC_2_gi{$NC} eq $gi);
			}
			
			$NC_2_gi{$NC} = $gi;
			
			
			push(@{$read_2_gis{$readID}}, [$gi, $mapQ]);
			push(@{$read_2_NCs{$readID}}, [$NC, $mapQ]);
			
			
			$n_reads++;
			$alignments_per_readID{$readID}++;
			$alignments_per_longReadID{$longReadID}++;
		
			if(exists $readID_2_length{$readID})
			{
				die unless($readID_2_length{$readID} = length($seq));
			}
			
			$readID_2_length{$readID} = length($seq);		

			
		}
		close(BAM);	
	}
	else
	{
		die "Not sure how to deal with config $config->[1]";
	}

	
	print "Reads - no alignments   - also in FASTQ:", scalar(grep {not exists $alignments_per_longReadID{$_}} keys %$HMP_readIDs_href), "\n";
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

	my %gis_present;
	my %NCs_present;
	foreach my $readID (keys %read_2_gis)
	{
		{
			my @alignments = @{$read_2_gis{$readID}};
			my $sortAlignments = sub {
				my $a = shift;
				my $b = shift;
				if($a->[1] == $b->[1])
				{
					return ($a->[0] cmp $b->[0]);
				}
				else
				{
					return ($a->[1] <=> $b->[1]);
				}
			};	
			if(scalar(@alignments) > 1)
			{
				@alignments = sort {$sortAlignments->($a, $b)} @alignments;
				@alignments = reverse @alignments;
				die unless($alignments[0][1] >= $alignments[1][1]);
			}
			
			$read_2_gis{$readID} = $alignments[0][0];
			$gis_present{$alignments[0][0]}++;
		}
		{
			my @alignments = @{$read_2_NCs{$readID}};
			my $sortAlignments = sub {
				my $a = shift;
				my $b = shift;
				if($a->[1] == $b->[1])
				{
					return ($a->[0] cmp $b->[0]);
				}
				else
				{
					return ($a->[1] <=> $b->[1]);
				}
			};	
			if(scalar(@alignments) > 1)
			{
				@alignments = sort {$sortAlignments->($a, $b)} @alignments;
				@alignments = reverse @alignments;
				die unless($alignments[0][1] >= $alignments[1][1]);
			}
			
			$read_2_NCs{$readID} = $alignments[0][0];
			$NCs_present{$alignments[0][0]}++;
		
		}
	}

	print "Reading gi-2-taxon...\n";
	unless(-e '/data/projects/phillippy/projects/mashsim/db/gi_taxid_nucl.dmp.HMP')
	{
		open(GI2TAXON, '<', '/data/projects/phillippy/projects/mashsim/db/gi_taxid_nucl.dmp') or die;
		open(GI2TAXONOUT, '>', '/data/projects/phillippy/projects/mashsim/db/gi_taxid_nucl.dmp.HMP') or die;
		while(<GI2TAXON>)
		{
			my $line = $_; 
			chomp($line);
			my @f = split(/\s+/, $line);
			die unless($#f == 1);	
			if($gis_present{$f[0]})
			{
				print GI2TAXONOUT $line, "\n";
			}
			if(($. % 100000) == 0)
			{
				print "\rGI line $. ...";
			}		
		}
		close(GI2TAXON);
		print "\n";
	}

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
	$gi_2_taxon{148642060} = '420247';	
	
	my %NC_2_taxon_final;
	foreach my $NC (keys %NCs_present)
	{
		(my $NC_noPoint = $NC) =~ s/\.\d+$//;
		if(defined $NC_2_taxon{$NC})
		{
			$NC_2_taxon_final{$NC} =  $NC_2_taxon{$NC};
		}
		elsif(defined $NC_noPoint_2_taxon{$NC_noPoint})
		{
			$NC_2_taxon_final{$NC} =  $NC_noPoint_2_taxon{$NC_noPoint};			
		}
		else
		{
			my $gi = $NC_2_gi{$NC};
			die unless(defined $gi);
			my $taxonID_fromGI = $gi_2_taxon{$gi};		
			die "No taxon info for gi $gi" unless(defined $taxonID_fromGI);
			die "Unexpected translation - $taxonID_fromGI" unless(($taxonID_fromGI eq '400667') or ($taxonID_fromGI eq '272943') or ($taxonID_fromGI eq '411466'));
			$NC_2_taxon_final{$NC} = $taxonID_fromGI;
		}
	}	 
	# print "\nGIs present: ", scalar(keys %gis_present), "\n";

	# gi 2 taxon ID

	# also see /data/projects/phillippy/projects/MetaMap/util/annotateHMPTruthTablesWithNCs.pl, where these files are being used

				
	open(OUT_PERREAD, '>', $fn_out_reads) or die "Cannot open file $fn_out_reads";
	my %read_2_taxonID;
	my %taxonID_read_counts;
	my %taxonID_2_bases;
	foreach my $readID (keys %read_2_gis)
	{
		my $gi = $read_2_gis{$readID};
		my $NC = $read_2_NCs{$readID};
		my $taxonID_original = $NC_2_taxon_final{$NC};
		die "No translation for NC number $NC" unless(defined $taxonID_original);
		
		my $taxonID_current = taxTree::findCurrentNodeID($master_taxonomy, $master_taxonomy_merged, $taxonID_original);
		print OUT_PERREAD join("\t", $readID, $taxonID_current), "\n";
		$read_2_taxonID{$readID} = $taxonID_current;
		$taxonID_read_counts{$taxonID_current}++;
		my $length = $readID_2_length{$readID};
		die unless(defined $length);
		
		my $taxonID_genomeLength_original = $gi_2_taxon{$gi};
		die "No translation for GI number $gi" unless(defined $taxonID_genomeLength_original);
		my $taxonID_genomeLength_current = taxTree::findCurrentNodeID($master_taxonomy, $master_taxonomy_merged, $taxonID_genomeLength_original);
		$taxonID_2_bases{$taxonID_genomeLength_current} += $length;
	}	

	foreach my $readID (keys %$HMP_readIDs_href)
	{
		next if(defined $read_2_taxonID{$readID});
		print OUT_PERREAD join("\t", $readID, 0), "\n";
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

	print "\n\nDone $config->[0] . Produced files:\n";
	print "\t - $fn_out_reads \n";
	print "\t - $fn_out_distribution \n";
	print "\t - $fn_out_distribution_genomeFreqs \n";
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