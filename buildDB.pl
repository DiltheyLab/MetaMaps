#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Find;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use List::Util qw/shuffle/;
use List::MoreUtils qw/all/;
use Cwd qw/abs_path getcwd/;

use taxTree;
use Util;

my $DB = 'refseq';
my $FASTAs;
my $taxonomyDir;
my $maxSpecies;
# my $Utest;
my $updateTaxonomy;
my $oldTaxonomy;
my $includeDirPattern;
GetOptions (
	'DB:s' => \$DB, 
	'FASTAs:s' => \$FASTAs, 
	'taxonomy:s' => \$taxonomyDir, 
	'maxSpecies:s' => \$maxSpecies, 
	'updateTaxonomy:s' => \$updateTaxonomy, 
	'oldTaxonomy:s' => \$oldTaxonomy, 
	'includeDirPattern:s' => \$includeDirPattern,
#	'Utest:s' => \$Utest, 
);

# Get input arguments

unless($DB and $FASTAs and $taxonomyDir)
{
	print_help();
}

die "Please specify a taxonomy directory (parameter --taxonomy)" unless(-d $taxonomyDir);
die "Taxonomy directory (--taxonomy) does not contain expected files for taxonomy in NCBI format (names.dmp, nodes.dmp..)." unless(all {-e $_} (map {$taxonomyDir . '/' . $_} taxTree::getTaxonomyFileNames()));

my @includeDirPatterns = split(/,/, $includeDirPattern);

unless(-d $DB)
{
	mkdir($DB) or die "Cannot mkdir $DB";
}

my @files_in_DB = glob($DB.'/*');
#die "Database output directory (parameter --DB) is not empty" unless(scalar(@files_in_DB) == 0); # todo reactivate

# copy taxonomy

my $dir_copy_taxonomy = $DB . '/taxonomy';
unless(-e $dir_copy_taxonomy)
{
	mkdir($dir_copy_taxonomy) or die "Cannot mkdir $dir_copy_taxonomy";
}

foreach my $f (taxTree::getTaxonomyFileNames())
{
	my $fP = $taxonomyDir . '/' . $f;
	my $cmd_cp = qq(cp $fP $dir_copy_taxonomy);
	system($cmd_cp) and die "Copy command $cmd_cp failed";
}

# Find input files

my @FASTAfiles;
foreach my $FASTAComponent (split(/,/, $FASTAs))
{
	if(-d $FASTAComponent)
	{
		find(\&wanted, $FASTAComponent);
		sub wanted {
			my $file = $File::Find::name;
			if($file =~ /(.fa$)|(\.fna$)/)
			{
				my $patternOK = 1;
				if(scalar(@includeDirPatterns))
				{
					$patternOK = 0;
					PATTERNLOOP: foreach my $pattern (@includeDirPatterns)
					{
						my $pattern_dir = '/' . $pattern . '/';
						if(index($file, $pattern_dir) != -1)
						{
							$patternOK = 1;
							last PATTERNLOOP;
						}
					}
				}
				# print join("\t", $file, $patternOK), "\n";
				if($patternOK)
				{
					push(@FASTAfiles, $file);
				}
			}
		}
	}
	else
	{
		die "Specified FASTA file (via --FASTAs) doesn't exist: $FASTAComponent" unless(-e $FASTAComponent);
		push(@FASTAfiles, $FASTAComponent);
	}
}
print "\nNumber of found FASTA input files: ", scalar(@FASTAfiles), "\n";

# read taxonomy

print "\nReading taxonomy from $taxonomyDir ..\n";
print "\tdone.\n\n";
my $taxonomy_href;

my $new_taxonomy_href;
my $new_taxonomy_href_merged;
if($updateTaxonomy)
{	
	die "Please specify --oldTaxonomy" 
	unless(defined $oldTaxonomy);
	
	$taxonomy_href = taxTree::readTaxonomy($oldTaxonomy);
	
	my $new_taxonomy_href_nonX = taxTree::readTaxonomy($taxonomyDir);
	$new_taxonomy_href_merged = taxTree::readMerged($taxonomyDir);
	
	if(scalar(grep {$_ =~ /^x/} keys %$new_taxonomy_href_nonX) > 0)
	{
		die "Your new --taxonomy has x nodes in it, abort";
	}	
	
	$new_taxonomy_href = taxTree::cloneTaxonomy_integrateX($new_taxonomy_href_nonX, $new_taxonomy_href_merged, $taxonomy_href);	
	
	taxTree::storeXInDir($new_taxonomy_href, $dir_copy_taxonomy);
}
else
{
	 $taxonomy_href = taxTree::readTaxonomy($taxonomyDir);		
}

# get list of contigs IDs and their positions			
			
my %taxonID_2_contig;	
my %contig_2_length;		
my %taxonIDs;
my @contigs;
for(my $fileI = 0; $fileI <= $#FASTAfiles; $fileI++)
{
	my $currentContigID;
	my $currentContigRunningLength;
		
	my $fileN = $FASTAfiles[$fileI];
	open(F, '<', $fileN) or die "Cannot open file $fileN";
	while(<F>)
	{
		if(substr($_, 0, 1) eq '>')
		{
			my $position_contig_start = tell(F) - length($_);			
			chomp;
			my $contigID = $_;
			
			push(@contigs, [$fileI, $position_contig_start, $contigID]);
			
			my $taxonID = Util::extractTaxonID($contigID, $fileN, $.);

			if($updateTaxonomy)
			{
				taxTree::findCurrentNodeID($new_taxonomy_href, $new_taxonomy_href_merged, $taxonID);
			}
			else
			{
				unless(exists $taxonomy_href->{$taxonID})
				{
					die "Taxon ID $taxonID - from file $fileN, line $. - not found in taxonomy -- consider updating your taxonomy directory.";
				}			
			}

			push(@{$taxonID_2_contig{$taxonID}}, substr($contigID, 1));
			
			if(defined $currentContigID)
			{
				die unless(defined $currentContigRunningLength);
				die "Duplicate contig ID in $currentContigID" if(defined $contig_2_length{$currentContigID});
				$contig_2_length{$currentContigID} = $currentContigRunningLength;
			}
			$currentContigID = substr($contigID, 1);
			$currentContigRunningLength = 0;
			
		}
		else
		{
			chomp;
			$currentContigRunningLength += length($_);
		}
	}
	close(F);
	if(defined $currentContigID)
	{
		die unless(defined $currentContigRunningLength);
		die "Duplicate contig ID in $currentContigID" if(defined $contig_2_length{$currentContigID});		
		$contig_2_length{$currentContigID} = $currentContigRunningLength;
	}	
}
@contigs = shuffle @contigs;

my $getContigSequence = sub {
	my $contigInfoAref = shift;
	
	my $fn = $FASTAfiles[$contigInfoAref->[0]];
	my $pos = $contigInfoAref->[1];
	my $expectedName = $contigInfoAref->[2];
	
	open(F, '<', $fn) or die "Cannot open file $fn";
	seek(F, $pos, 0);
	my $firstLine = <F>;
	chomp($firstLine);
	die "Couldn't find contig - got $firstLine, wanted $expectedName, in $fn" unless($firstLine eq $expectedName);
	
	if($updateTaxonomy)
	{
		$firstLine = update_contigID_newTaxon_ID($firstLine);
	}
	
	my $contigSequence = $firstLine . "\n";
	
	while(<F>)
	{
		if(substr($_, 0, 1) eq '>')
		{
			last;
		}
		else
		{
			$contigSequence .= $_;
		}
	}
	close(F);
	
	return $contigSequence;
};

my @useTaxonIDs = keys %taxonID_2_contig;	
if((defined $maxSpecies) and (scalar(@useTaxonIDs) > $maxSpecies))
{
	die unless($maxSpecies > 0);
	@useTaxonIDs = @useTaxonIDs[0 .. ($maxSpecies-1)];
}

# none of that stuff works with updated taxonomies
# if($Utest)
# {
	# my $count_ancestors_with_genomes = sub {
		# my $nodeID = shift;
		
		# die unless(defined $taxonomy_href->{$nodeID}{children});
		
		# my @children_nodes = $taxonomy_href->{$nodeID}{children};
		
		# my @children_nodes_mappable = grep {
			# my $child_node_ID = $_;
			
			# my @descendants_of_child = taxTree::descendants($taxonomy_href, $child_node_ID);
			# my @child_and_descendants = ($child_node_ID, @descendants_of_child);
			
			# my @descendants_mappable = grep {exists $taxonID_2_contig{$_}} @child_and_descendants;
			
			# (scalar(@descendants_mappable) > 0);
		# } @children_nodes;
		
		# return scalar(@children_nodes);
	# };
	
	# foreach my $level (qw/species genus family/)
	# {
		# my @nodes_thisLevel = grep {$taxonomy_href->{$_}{rank} eq $level} keys %$taxonomy_href;
		# my @nodes_thisLevel_atLeast3 = grep {$count_ancestors_with_genomes->($_) >= 3} @nodes_thisLevel;
		# print "Mappable nodes $level: ", scalar(@nodes_thisLevel_atLeast3), "\n";
	# }
	
	# exit;
# }


my %_useTaxonID = map {$_ => 1} @useTaxonIDs;
my %_useTaxonIDs_afterUpdate;
my $outputTaxonsAndContigs = $DB . '/' . 'taxonInfo.txt';
open(DB, '>', $outputTaxonsAndContigs) or die "Cannot open $outputTaxonsAndContigs - check whether I have write permissions";
if($updateTaxonomy)
{
	my %new_taxonID_2_contigs;
	my %new_contig_2_length;
	my $made_updates_taxonID = 0;
	my $made_updates_contigID = 0;
	
	foreach my $taxonID (@useTaxonIDs)
	{
		my $newTaxonID = taxTree::findCurrentNodeID($new_taxonomy_href, $new_taxonomy_href_merged, $taxonID);
		$made_updates_taxonID++ if ($taxonID ne $newTaxonID);		
		my @contigs = @{$taxonID_2_contig{$taxonID}};
		my @new_contigs;
		foreach my $contigID (@contigs)
		{
			die unless(Util::extractTaxonID($contigID, '?', '?') eq $taxonID);
			my $newContigID = update_contigID_newTaxon_ID($contigID);
			$made_updates_contigID++ if ($contigID ne $newContigID);
			die if(defined $new_contig_2_length{$newContigID});
			$new_contig_2_length{$newContigID} = $contig_2_length{$contigID};
			push(@new_contigs, $newContigID);
		}
		push(@{$new_taxonID_2_contigs{$newTaxonID}}, @new_contigs);
	}
	
	my @new_useTaxonIDs = keys %new_taxonID_2_contigs;
	foreach my $taxonID (@new_useTaxonIDs)
	{
		my @new_contigs = @{$new_taxonID_2_contigs{$taxonID}};
		print DB $taxonID, " ", join(';', map {die unless(defined $new_contig_2_length{$_}); $_ . '=' . $new_contig_2_length{$_}} @new_contigs), "\n";
	}
	
	%_useTaxonIDs_afterUpdate = map {$_ => 1} @new_useTaxonIDs;
	
	print "Updated taxonomic information -- $made_updates_taxonID taxons, $made_updates_contigID contig IDs.\n";
}
else
{
	foreach my $taxonID (@useTaxonIDs)
	{
		my @contigs = @{$taxonID_2_contig{$taxonID}};
		print DB $taxonID, " ", join(';', map {die unless(defined $contig_2_length{$_}); $_ . '=' . $contig_2_length{$_}} @contigs), "\n";
	}
	
	%_useTaxonIDs_afterUpdate = map {$_ => 1} @useTaxonIDs;
}
close(DB);

my $countNs_windowSize = 1000;
my $outputFN_contigWindows = $DB . '/' . 'contigNstats_windowSize_' . $countNs_windowSize . '.txt';

my $outputFN = $DB . '/' . 'DB.fa';
open(DB, '>', $outputFN) or die "Cannot open $outputFN - check whether I have write permissions";
open(NSTATS, '>', $outputFN_contigWindows) or die "Cannot open $outputFN_contigWindows - check whether I have write permissions";
for(my $contigI = 0; $contigI <= $#contigs; $contigI++)
{
	my $contigInfo = $contigs[$contigI];
	my $contigID = $contigInfo->[2];
	my $taxonID = Util::extractTaxonID($contigID, '?', '?');
	if($_useTaxonID{$taxonID})
	{
		my $contigSequence = $getContigSequence->($contigInfo);
		print DB $contigSequence;
		
		my $contigID_noTrailing = substr($contigID, 1);
		
		die unless(substr($contigSequence, 0, 1) eq '>');
		(my $contigSequence_contiguous = $contigSequence) =~ s/^(.+)//;
		my $translatedContigID = $1;
		die unless(substr($translatedContigID, 0, 1) eq '>');
		$translatedContigID = substr($translatedContigID, 1);
		$contigSequence_contiguous =~ s/[\n\r]//g;
		
		my $taxonID_fromTranslatedContigID = Util::extractTaxonID($translatedContigID, '?', '?');

	
		die Dumper("Length mismatch", $contigID, length($contigSequence_contiguous), $contig_2_length{$contigID_noTrailing}, substr($contigSequence_contiguous, 0, 200))  unless(length($contigSequence_contiguous) == $contig_2_length{$contigID_noTrailing});
		
		my @windows_Ns;
		for(my $startWindow = 0; $startWindow < length($contigSequence_contiguous); $startWindow += $countNs_windowSize)
		{
			my $stopWindow = $startWindow + $countNs_windowSize  - 1;
			$stopWindow = (length($contigSequence_contiguous) - 1) if ($stopWindow >= length($contigSequence_contiguous));
			die unless($stopWindow >= $startWindow);
			my $windowSequence = substr($contigSequence_contiguous, $startWindow, $stopWindow - $startWindow + 1);
			
			my $count_n = ($windowSequence =~ tr/n//);			
			my $count_N = ($windowSequence =~ tr/N//);			
			
			push(@windows_Ns, $count_n + $count_N);
		}
		
		print NSTATS join("\t", $taxonID_fromTranslatedContigID, $translatedContigID, join(";", @windows_Ns)), "\n";
	}
}
close(DB);
close(NSTATS);

taxTree::trimTaxonomyInDir($dir_copy_taxonomy, \%_useTaxonIDs_afterUpdate);

print "\nProduced randomized-order database sequence file $outputFN and\n\ttaxon info file $outputTaxonsAndContigs (" . scalar(@useTaxonIDs) . " taxa).\n\n";


sub update_contigID_newTaxon_ID
{
	my $contigID = shift;
	
	my $oldContigID = $contigID;
	
	my $existingTaxonID = Util::extractTaxonID($contigID, '?', '?');
	my $newTaxonID = taxTree::findCurrentNodeID($new_taxonomy_href, $new_taxonomy_href_merged, $existingTaxonID);
	
	$contigID =~ s/kraken:taxid\|(x?\d+)/kraken:taxid\|$newTaxonID/;
	if($existingTaxonID ne $newTaxonID)
	{
		die if(index($contigID, $existingTaxonID) != -1);
	}
	
	return $contigID;
}

sub print_help
{
	print qq(
buildDB.pl

  Assemble a MetaMap database.
  
Usage:

  perl buildDB.pl --DB dbNAME --FASTAs DIR|FILELIST --taxonomy DIR
  
Parameters:

  DB
      Name (output directory) of DB to be constructed.
  
  FASTAs
  
      Directory with or comma-separated list of FASTA input
	  files. Directories are scanned recursively (*.fa, *.fna).
      Each contig has to be annotated with a taxon ID
      satisfying the regular expression 'kraken:taxid\\|(x?\\d+)'.
      
  taxonomy
      
      Path to taxonomy directory (in NCBI format, containing
      files nodes.dmp, names.dmp etc.) - all taxon IDs present
      in the FASTA must be present in the taxonomy, including
      pseudo IDs starting with an 'x'.
  
);
exit;
}
