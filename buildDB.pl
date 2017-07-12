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

my $DB = 'refseq';
my $FASTAs;
my $taxonomyDir;

GetOptions (
	'DB:s' => \$DB, 
	'FASTAs:s' => \$FASTAs, 
	'taxonomy:s' => \$taxonomyDir, 
);

# Get input arguments

unless($DB and $FASTAs and $taxonomyDir)
{
	print_help();
}


die "Please specify a taxonomy directory (parameter --taxonomy)" unless(-d $taxonomyDir);
die "Taxonomy directory (--taxonomy) does not contain expected files for taxonomy in NCBI format (names.dmp, nodes.dmp..)." unless(all {-e $_} (map {$taxonomyDir . '/' . $_} taxTree::getTaxonomyFileNames()));

unless(-d $DB)
{
	mkdir($DB) or die "Cannot mkdir $DB";
}

my @files_in_DB = glob($DB.'/*');
die "Database output directory (parameter --DB) is not empty" unless(scalar(@files_in_DB) == 0);

# copy taxonomy

my $dir_copy_taxonomy = $DB . '/taxonomy';
mkdir($dir_copy_taxonomy) or die "Cannot mkdir $dir_copy_taxonomy";

foreach my $f (taxTree::getTaxonomyFileNames())
{
	my $fP = $taxonomyDir . '/' . $f;
	my $cmd_cp = qq(cp $fP $dir_copy_taxonomy);
	system($cmd_cp) and die "Copy command $cmd_cp failed";
}
# Find input files

my @FASTAfiles;
if(-d $FASTAs)
{
	find(\&wanted, $FASTAs);
	sub wanted {
		my $file = $File::Find::name;
		if($file =~ /(.fa$)|(\.fna$)/)
		{
			push(@FASTAfiles, $file);
		}
	}
}
else
{
	@FASTAfiles = split(/,/, $FASTAs);
	foreach my $f (@FASTAfiles)
	{
		die "Specified FASTA file (via --FASTAs) doesn't exist: $f" unless(-e $f);
	}
}

print "\nNumber of found FASTA input files: ", scalar(@FASTAfiles), "\n";

# read taxonomy

print "\nReading taxonomy from $taxonomyDir ..\n";
my $taxonomy_href = taxTree::readTaxonomy($taxonomyDir);
print "\tdone.\n\n";

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
			
			unless($contigID =~ /kraken:taxid\|(x?\d+)/)
			{
				die "Expect taxon ID in contig identifier - file $fileN - line $.";
			}			
			my $taxonID = $1;
			unless(exists $taxonomy_href->{$taxonID})
			{
				die "Taxon ID $taxonID - from file $fileN, line $. - not found in taxonomy -- consider updating your taxonomy directory.";
			}
			
			push(@{$taxonID_2_contig{$taxonID}}, substr($contigID, 1));
			
			if(defined $currentContigID)
			{
				die unless(defined $currentContigRunningLength);
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

my $outputTaxonsAndContigs = $DB . '/' . 'taxonInfo.txt';
open(DB, '>', $outputTaxonsAndContigs) or die "Cannot open $outputTaxonsAndContigs - check whether I have write permissions";
foreach my $taxonID (keys %taxonID_2_contig)
{
	my @contigs = @{$taxonID_2_contig{$taxonID}};
	
	print DB $taxonID, " ", join(';', map {die unless(defined $contig_2_length{$_}); $_ . '=' . $contig_2_length{$_}} @contigs), "\n";
}
close(DB);

my $outputFN = $DB . '/' . 'DB.fa';
open(DB, '>', $outputFN) or die "Cannot open $outputFN - check whether I have write permissions";
for(my $contigI = 0; $contigI <= $#contigs; $contigI++)
{
	my $contigInfo = $contigs[$contigI];
	my $contigSequence = $getContigSequence->($contigInfo);
	print DB $contigSequence;
}
close(DB);

print "\nProduced randomized-order database sequence file $outputFN and\n\ttaxon info file $outputTaxonsAndContigs\n\n";

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
