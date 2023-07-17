#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use File::Find;
use lib "$FindBin::Bin/perlLib";
use Cwd qw/abs_path/;

use taxTree;

my $inputFileList;
my $outputFile;
my $taxonomyInDirectory;
my $taxonomyOutDirectory;
GetOptions (
	'inputFileList:s' => \$inputFileList, 
	'outputFile:s' => \$outputFile, 
	'taxonomyInDirectory:s' => \$taxonomyInDirectory, 
	'taxonomyOutDirectory:s' => \$taxonomyOutDirectory, 
);

unless($inputFileList and $outputFile and $taxonomyInDirectory and $taxonomyOutDirectory)
{
	print_help();
}

die "Please specify valid file for --inputFileList" unless(-e $inputFileList);
die "Please specify valid directory for --taxonomyInDirectory" unless(-d $taxonomyInDirectory);

# get list of input files and rudimentary checks

my %taxonID_2_files;
my $n_files = 0;
open(INPUT, '<', $inputFileList) or die "Cannot open $inputFileList";
while(<INPUT>)
{
	chomp;
	$_ =~ s/[\r\n]//g;
	next unless($_);
	my @fields = split(/\s/, $_);
	die "Input file $inputFileList, line $.: Invalid format - I expect two whitespace-separated files, the first specifying a valid taxon ID, the second specifying a path to a FASTA file" unless(scalar(@fields) == 2);
	die "Input file $inputFileList, line $.: The specified path ($fields[1]) does not exist" unless(-e $fields[1]);
	open(F, '<', $fields[1]) or die "Cannot open $fields[1]";
	my $firstLine = <F>;
	close(F);
	unless(substr($firstLine, 0, 1) eq '>')
	{
		die "The specified file $fields[1] does not seem to be a valid FASTA file (criterion: the first character of the first line is a '>')";
	}
	push(@{$taxonID_2_files{$fields[0]}}, $fields[1]);
	$n_files++;
}
close(INPUT);

print "Input: $n_files genomes and ", scalar(keys %taxonID_2_files), " taxon IDs.\n";

# Read taxonomy and check the specified IDs are valid

print "Reading taxonomy from $taxonomyInDirectory ..\n";
my $taxonomy_href = taxTree::readTaxonomy($taxonomyInDirectory);
my $merged_href = taxTree::readMerged($taxonomyInDirectory);
print "\tdone.\n\n";

my %invalid_taxon_IDs;
foreach my $taxonID (keys %taxonID_2_files)
{
	if(not exists $taxonomy_href->{$taxonID})
	{
		$invalid_taxon_IDs{$taxonID}++;
	}
}
if(scalar(keys %invalid_taxon_IDs))
{
	die "Several taxon IDs specified in $inputFileList are invalid, i.e. they don't appear in the taxonomy specified in $taxonomyInDirectory: " . join("\n", keys %invalid_taxon_IDs);
}

# All rudimentary checks passed - copy taxonomy

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

	my @associated_files = @{$taxonID_2_files{$taxonID}};	
	if(scalar(@associated_files) > 1)
	{
		my $thisNode_rank = $taxonomy_href->{$taxonID}{rank};
		die unless(defined $thisNode_rank);
		unless(($thisNode_rank eq 'species') or ($thisNode_rank eq 'no rank') or ($thisNode_rank eq 'subspecies') or ($thisNode_rank eq 'varietas'))
		{
			die Dumper("Unexpected rank", $thisNode_rank, $taxonID, $taxonID_2_files{$taxonID});
		}	
		foreach my $f (@associated_files)
		{
			$running_new_ids++; 
			my $newID = 'x' . $running_new_ids;
			$file_2_taxonID{$f} = $newID;
			$new_tree_relationships{$newID} = $taxonID;
			$new_tree_names{$newID} = $taxonomy_href->{$taxonID}{names}[0];
			$taxonID_2_files{$newID} = [$f];
		}
		
		delete $taxonID_2_files{$taxonID};
	}
	else
	{
		my $f = $associated_files[0];
		$file_2_taxonID{$f} = $taxonID;	
	}
}

print "Introduced $running_new_ids new taxonomic IDs\n";

open(OUTPUT, '>', $outputFile) or die "Cannot open $outputFile for writing";
my $contigCounter = 0;
foreach my $FASTAfile (keys %file_2_taxonID)
{
	my $taxonID = $file_2_taxonID{$FASTAfile};
	die unless(defined $taxonID);
	open(F, '<', $FASTAfile) or die "Cannot open $FASTAfile";
	while(<F>)
	{
		my $line = $_;
		if(substr($line, 0, 1) eq '>')
		{
			substr($line, 0, 1) = '';
			my $line_for_inclusion = $line;
			
			if($line_for_inclusion =~ /kraken:taxid\|(\d+)/)
			{
				warn "File $FASTAfile already conains kraken segment -- this information will be ignored.";
				$line_for_inclusion =~ s/kraken:taxid\|(\d+)//g;
			}
			$line_for_inclusion =~ s/\s.+//;

			$contigCounter++;
			
			my $newID = 'C' . $contigCounter . '|kraken:taxid|' . $taxonID . '|' . $line_for_inclusion;
			print OUTPUT '>', $newID;
		}
		else
		{
			print OUTPUT $line;
		}
	}
	close(F);
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

sub print_help
{
	print qq(
combineAndAnnotateReferences.pl

  Prepares multiple reference genomes in FASTA format with specified taxon
  IDs for use with MetaMaps.
  
Usage:

  perl combineAndAnnotateReferences.pl --inputFileList FILE --outputFile FILE --taxonomyInDirectory DIR --taxonomyOutDirectory DIR
  
Details:

  Creates a single, taxon-annotated output FASTA that contains the genomes specified
  via --inputFileList. The produced FASTA can be passed to buildDB.pl.
  
  Taxon ID information is taken from the list of input files; any existing taxon ID
  information in the input FASTA files is overwritten or ignored.
  
  The input file list follows a simple space-separated format:
  
  taxon_id_1 path_to_genome_1_FASTA
  taxon_id_2 path_to_genome_2_FASTA
  ...
  
  Each input file is treated as a SEPARATE genome; combine multi-FASTA reference
  genomes into single files for use with this script. If multiple input genomes
  carry the same taxon ID, these will be treated as strains and unique MetaMaps-
  specific taxon IDs (prefixed with an 'x') will be created as necessary.
  
  All specified taxon IDs need to be present in the taxonomy (e.g. downloaded
  from NCBI) specified with --taxonomyInDirectory. An output taxonomy (including
  all utilizied and created taxon IDs) is written into  --taxonomyOutDirectory.
  Do not use the same directory for input and output taxonomies.
  
  The FASTA combination process will generate unique contig identifiers; not all
  contig header information from the input files will be preserved in the combined
  output file.
);
 
exit; 
}