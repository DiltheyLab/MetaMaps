#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use File::Find;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/abs_path/;

use taxTree;

my $refSeqDirectory;
my $taxonomyInDirectory;
my $taxonomyOutDirectory;
GetOptions (
	'refSeqDirectory:s' => \$refSeqDirectory, 
	'taxonomyInDirectory:s' => \$taxonomyInDirectory, 
	'taxonomyOutDirectory:s' => \$taxonomyOutDirectory, 
);

unless($refSeqDirectory and $taxonomyInDirectory and $taxonomyOutDirectory)
{
	print_help();
}

$refSeqDirectory = abs_path($refSeqDirectory);

die "Please specify valid directory for --refSeqDirectory" unless(-d $refSeqDirectory);
die "Please specify valid directory for --taxonomyInDirectory" unless(-d $taxonomyInDirectory);

print "Reading taxonomy from $taxonomyInDirectory ..\n";
my $taxonomy_href = taxTree::readTaxonomy($taxonomyInDirectory);
my $merged_href = taxTree::readMerged($taxonomyInDirectory);
print "\tdone.\n\n";

my @expected_taxonomy_files = taxTree::getTaxonomyFileNames();
unless(-d $taxonomyOutDirectory)
{
	mkdir($taxonomyOutDirectory) or die "Can't mkdir $taxonomyOutDirectory (specified via --taxonomyOutDirectory)";
}
foreach my $f (@expected_taxonomy_files)
{
	my $cmd_cp_f = qq(cp ${taxonomyInDirectory}/${f} ${taxonomyOutDirectory}/);
	system($cmd_cp_f) and die "Cannot execute cp command:\n$cmd_cp_f";
}

my %file_2_taxonID;
my %file_2_organismName;
my %taxonID_2_files;

my %refSeq_categories;
my %assemblyLevels;
my %category_level_combined;
my $included_genome = 0;
my $missed_representative = 0;
my $missed_referenceGenome = 0;
my $missed_reference_and_representative = 0;
print "Scanning $refSeqDirectory for *_assembly_report.txt\n";
find(\&wanted, $refSeqDirectory);
sub wanted {
	my $file = $File::Find::name; 
	
	{	
		no warnings;
		$file = $File::Find::name; 
	}
	
	# warn Dumper("Expected file nt existing", $file) unless(-e $file);
	
	if((-e $file) and (not -d $file) and ($file =~ /_assembly_report\.txt$/))
	{	
		my $gzFile = get_gz_for_assemblyReport($file);
		if(-e $gzFile)
		{
			remove_existing_fnaFile_for_assemblyReportFile($file);
			
			my $taxonID;
			my $organismName = 'NA';
			my $refSeq_category = '';
			my $assembly_level = '';
			open(REPORT, '<', $file) or die "Cannot open assembly report $file";
			while(<REPORT>)
			{
				if($_ =~ /^# Taxid:\s*(\d+)/)
				{
					$taxonID = $1;
				}
				
				if($_ =~ /^# Organism name:\s*(.+)/)
				{
					$organismName = $1;
					$organismName =~ s/[\n\r]//g;
				}				
				
				if($_ =~ /^# RefSeq category:\s*(.+)/)
				{
					$refSeq_category = $1;
					$refSeq_category =~ s/[\n\r]//g;
				}				
				
				if($_ =~ /^# Assembly level:\s*(.+)/)
				{
					$assembly_level = $1;
					$assembly_level =~ s/[\n\r]//g;
				}				
								
			}
			close(REPORT);
			unless(defined $taxonID)
			{
				die "Assembly report $file has not taxon ID";
			}
			
			$refSeq_categories{$refSeq_category}++;
			$assemblyLevels{$assembly_level}++;
			if(not $refSeq_category)
			{
				$refSeq_category = 'Undefined';
				# warn "No RefSeq category info in $file";
			}
			$category_level_combined{join('_', $refSeq_category, $assembly_level)}++;

			next unless($assembly_level eq 'Complete Genome');

			$included_genome++;
			
			unless(exists $taxonomy_href->{$taxonID})
			{
				warn "Taxon ID $taxonID not defined in tree in $taxonomyInDirectory - try recovering from merged nodes.";
				$taxonID = taxTree::findCurrentNodeID($taxonomy_href, $merged_href, $taxonID);
			}
			if(exists $taxonomy_href->{$taxonID})
			{
				$file_2_taxonID{$file} = $taxonID;
				$file_2_organismName{$file} = $organismName;
				push(@{$taxonID_2_files{$taxonID}}, $file);			
			}
			else
			{
				die "Taxon ID $taxonID not defined in tree in $taxonomyInDirectory -- update your taxonomy directory?";
			}
				
			# if(1) #  or $assembly_level eq 'Complete Genome'
			# {

			# }
			# else
			# {	
				# die Dumper("Weird assembly_level", $assembly_level) unless(($assembly_level eq 'Contig') or ($assembly_level eq 'Scaffold') or ($assembly_level eq 'Chromosome'));
				# if(($refSeq_category eq 'Representative Genome') and ($refSeq_category eq 'Reference Genome'))
				# {
					# $missed_reference_and_representative++;
				# }
				# elsif($refSeq_category eq 'Representative Genome')
				# {
					# $missed_representative++;
				# }
				# elsif($refSeq_category eq 'Reference Genome')
				# {
					# $missed_referenceGenome++;
				# }
				# else
				# {
					# die Dumper("Weird refSeq_category", $refSeq_category) unless(($refSeq_category eq '') or ($refSeq_category eq 'Representative Genome'));
				# }
			# }
		}
		else
		{
			warn "No assembly data file - $gzFile";
		}
	}
}

print "Summary input data:\n";
print "\tRefseq categories:\n";
foreach my $key (keys %refSeq_categories)
{
	print "\t\t", $key, ": ", $refSeq_categories{$key}, "\n";
}
print "\tAssembly levels:\n";
foreach my $key (keys %assemblyLevels)
{
	print "\t\t", $key, ": ", $assemblyLevels{$key}, "\n";
}
print "\t{Category} X {Assembly level}:\n";
foreach my $key (keys %category_level_combined)
{
	print "\t\t", $key, ": ", $category_level_combined{$key}, "\n";
}

print "Total genomes: $included_genome \n";

# print "\tMissed reference and representative: $missed_reference_and_representative \n";
# print "\tMissed representative: $missed_representative \n";
# print "\tMissed reference: $missed_referenceGenome \n";

my $running_new_ids = 0;

my %new_tree_relationships;
my %new_tree_names;
my %need_taxonID = (1 => 1);
foreach my $taxonID (keys %taxonID_2_files)
{
	unless(exists $taxonomy_href->{$taxonID})
	{
		die "Taxon ID $taxonID not defined in tree in $taxonomyInDirectory -- update your taxonomy directory?";
	}
	
	my $current_id = $taxonID;
	do {
		$need_taxonID{$current_id} = 1;
		die unless(defined $taxonomy_href->{$current_id});	
		$current_id = $taxonomy_href->{$current_id}{parent};
	} while($current_id != 1);	

	if(scalar(@{$taxonID_2_files{$taxonID}}) > 1)
	{
		my $thisNode_rank = $taxonomy_href->{$taxonID}{rank};
		die unless(defined $thisNode_rank);
		unless(($thisNode_rank eq 'species') or ($thisNode_rank eq 'no rank') or ($thisNode_rank eq 'subspecies') or ($thisNode_rank eq 'varietas'))
		{
			die Dumper("Unexpected rank", $thisNode_rank, $taxonID, $taxonID_2_files{$taxonID});
		}	
		my @associated_files = @{$taxonID_2_files{$taxonID}};
		foreach my $f (@associated_files)
		{
			$running_new_ids++; 
			my $newID = 'x' . $running_new_ids;
			$file_2_taxonID{$f} = $newID;
			$new_tree_relationships{$newID} = $taxonID;
			$new_tree_names{$newID} = $file_2_organismName{$f};
			$taxonID_2_files{$newID} = [$f];
		}
		
		delete $taxonID_2_files{$taxonID};
	}
}

print "Introduced $running_new_ids new taxonomic IDs\n";

my $contigCounter = 0;
foreach my $assemblyReportFile (keys %file_2_taxonID)
{
	my $taxonID = $file_2_taxonID{$assemblyReportFile};
	my $fresh_fnaFile = get_fresh_fnaFile_for_assemblyReportFile($assemblyReportFile);
	my $f2 = $fresh_fnaFile . '_';
	open(F, '<', $fresh_fnaFile) or die "Cannot open $fresh_fnaFile";
	open(F2, '>', $f2) or die "Cannot open $f2";
	while(<F>)
	{
		my $line = $_;
		if(substr($line, 0, 1) eq '>')
		{
			substr($line, 0, 1) = '';
			die "File $fresh_fnaFile already conains kraken segment?" if ($line =~ /kraken:taxid\|(\d+)/);
			$contigCounter++;
			my $line_for_inclusion = $line;
			$line_for_inclusion =~ s/\s.*//;
			$line_for_inclusion =~ s/[\r\n]//g;			
			my $newID = 'C' . $contigCounter . '|kraken:taxid|' . $taxonID . '|' . $line_for_inclusion;
			die "Invalid ID: $newID" if ($newID =~ /\s/);
			print F2 '>', $newID, "\n";
		}
		else
		{
			print F2 $line;
		}
	}
	close(F2);
	close(F);
	
	my $mv_command = qq(mv $f2 $fresh_fnaFile);
	die "Command $mv_command failed" if (system($mv_command));
}

print "Annotated $contigCounter contigs\n";
my $taxonomy_names_f_out = $taxonomyOutDirectory . '/names.dmp';
my $taxonomy_nodes_f_out = $taxonomyOutDirectory . '/nodes.dmp';

die unless(-e $taxonomy_names_f_out);
open(NAMESOUT, '>>', $taxonomy_names_f_out) or die "Cannot open $taxonomy_names_f_out";
foreach my $newID (keys %new_tree_relationships)
{
	print NAMESOUT join("\t|\t", $newID, $new_tree_names{$newID}, '', 'scientific name', ''), "\n";
}
close(NAMESOUT);

die unless(-e $taxonomy_nodes_f_out);
open(NODESOUT, '>>', $taxonomy_nodes_f_out) or die "Cannot open $taxonomy_nodes_f_out";
foreach my $newID (keys %new_tree_relationships)
{
	print NODESOUT join("\t|\t", $newID,  $new_tree_relationships{$newID}, 'pseudospecies', ''), "\n";
}	
close(NODESOUT);

print "\nOutput new taxonomy into $taxonomyOutDirectory\n\n";

sub get_gz_for_assemblyReport
{
	my $assemblyReport = shift;
	die "Weird apparent assembly report file $assemblyReport" unless($assemblyReport =~ /_assembly_report\.txt$/);
	
	my $gzFile = $assemblyReport;
	$gzFile =~ s/_assembly_report\.txt/_genomic.fna.gz/;
	
	return $gzFile;
}

sub remove_existing_fnaFile_for_assemblyReportFile
{
	my $assemblyReport = shift;
	die "Weird apparent assembly report file $assemblyReport" unless($assemblyReport =~ /_assembly_report\.txt$/);
	
	my $gzFile = get_gz_for_assemblyReport($assemblyReport);
	
	unless(-e $gzFile)
	{
		warn "Expected file $gzFile  not existing";
		return 0;
	}
	
	my $fna_file = $gzFile;
	$fna_file =~ s/\.gz$//;
	if(-e $fna_file)
	{
		unlink($fna_file) or die "Cannot unlink $fna_file";
	}			
	
	return 1;
}

sub get_fresh_fnaFile_for_assemblyReportFile
{
	my $assemblyReport = shift;
	die "Weird apparent assembly report file $assemblyReport" unless($assemblyReport =~ /_assembly_report\.txt$/);
	
	my $gzFile = $assemblyReport;
	$gzFile =~ s/_assembly_report\.txt/_genomic.fna.gz/;
	die "Expected file $gzFile  not existing" unless(-e $gzFile);
	
	my $fna_file = $gzFile;
	$fna_file =~ s/\.gz$//;
	if(-e $fna_file)
	{
		unlink($fna_file) or die "Cannot unlink $fna_file";
	}
	my $cmd = qq(gunzip -c $gzFile > $fna_file);
	if(system($cmd))
	{
		unlink($fna_file);
		die "Command $cmd failed";
	}

	die unless(-e $fna_file);
	
	return $fna_file;			
}



sub print_help
{
	print qq(
annotateRefSeqSequencesWithUniqueTaxonIDs.pl

  Prepares a RefSeq/GenBank download for MetaMaps.
  
Usage:

  perl annotateRefSeqSequencesWithUniqueTaxonIDs.pl --refSeqDirectory DIR --taxonomyInDirectory DIR --taxonomyOutDirectory DIR
  
Details:

  Reads a directory structure as downloaded from RefSeq / Genbank servers,
  unzips assemblies, gets information from *_assembly_report files, and writes
  FASTA files with taxonID information in FASTA IDs.
  
  Each sub-directory is treated as a separate 'mappping unit' that will receive a
  unique mapping ID - the script will create pseudo taxonomic IDs, prefixed with an
  'x', if necessary.
  
  The specified --refSeqDirectory is manipulated in-place, and the amended taxonomy
  is written into --taxonomyOutDirectory\n\n
  
);
 
exit; 
}