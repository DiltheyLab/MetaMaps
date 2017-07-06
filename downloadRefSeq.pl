#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Net::FTP;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use List::MoreUtils qw/all/;
use Cwd qw/abs_path getcwd/;

my $DB = 'refseq';
my $seqencesOutDirectory = '';
my $taxonomyOutDirectory;

GetOptions (
	'DB:s' => \$DB, 
	'seqencesOutDirectory:s' => \$seqencesOutDirectory, 
	'taxonomyOutDirectory:s' => \$taxonomyOutDirectory, 
);

# Get input arguments

unless($seqencesOutDirectory and $taxonomyOutDirectory)
{
	print_help();
}

die "Please specify --DB as either refseq or genbank" unless(($DB eq 'refseq') or ($DB eq 'genbank'));

# Create output directories

if(not -d $seqencesOutDirectory)
{
	mkdir($seqencesOutDirectory) or die "Cannot mkdir $seqencesOutDirectory";
}

if(not -d $taxonomyOutDirectory)
{
	mkdir($taxonomyOutDirectory) or die "Cannot mkdir $taxonomyOutDirectory";
}

# initial cwd
my $cwd = getcwd();

# Establish FTP

my $ftp_server = 'ftp.ncbi.nlm.nih.gov';
my $ftp = Net::FTP->new($ftp_server, Debug => 0) or die "Cannot connect to some.host.name: $@";
$ftp->login("anonymous",'-anonymous@') or die "Cannot login ", $ftp->message;
$ftp->binary();

# Download taxonomy
use taxTree;
my @expected_taxonomy_files = taxTree::getTaxonomyFileNames();
chdir($taxonomyOutDirectory) or die "Cannot chdir into $taxonomyOutDirectory";
foreach my $f (@expected_taxonomy_files)
{
	if(-e $f)
	{
		unlink($f) or die "Cannot delete $f -- check I have write access to $taxonomyOutDirectory";
	}
}
my $ftp_root_taxonomy = '/pub/taxonomy/';
$ftp->cwd($ftp_root_taxonomy) or die "Cannot change working directory ", $ftp->message;
my @taxonomy_content = $ftp->ls();
die "No file taxdump.tar.gz in $ftp_root_taxonomy on $ftp_server?" unless(scalar(grep {$_ eq 'taxdump.tar.gz'} @taxonomy_content) == 1);
$ftp->get('taxdump.tar.gz') or die "Cannot transfer file", $ftp->message;		
my $cmd_extract = qq(tar -xvzf taxdump.tar.gz);
system($cmd_extract) and die "Could not extract taxonomy file taxdump.tar.gz";
die "Expected taxonomy files missing" unless(all {-e $_} @expected_taxonomy_files);

print "\nTaxonomy downloaded and extracted into $taxonomyOutDirectory\n\n";

# Download species

chdir($cwd) or die "Cannot chdir into $cwd";

my @target_subdirs = qw/archaea bacteria fungi protozoa viral invertebrate plant unknown vertebrate_other vertebrate_mammalian/;
#my @target_subdirs = qw/bacteria fungi protozoa viral invertebrate plant unknown vertebrate_other vertebrate_mammalian/;

my $ftp_root_genomes = '/genomes'; die unless(substr($ftp_root_genomes, 0, 1) eq '/');

my $DB_downloaded_assemblies = 0;

$ftp->cwd($ftp_root_genomes . '/' . $DB) or die "Cannot change working directory ", $ftp->message;

foreach my $subDir (@target_subdirs)
{
	my $subDir_local = $seqencesOutDirectory . '/' . $subDir;
	mkdir($subDir_local);
	die "Directory $subDir_local not existing, but it should (I just tried to mkdir it - do I have write permissions?"  unless(-d $subDir_local);
	
	$ftp->cwd($ftp_root_genomes . '/' . $DB . '/'. $subDir . '/') or die "Cannot change working directory ", $ftp->message;		
	
	my @subDir_content = $ftp->ls();
	
	print "Now download genomes for ", scalar(@subDir_content), " $subDir species (", $DB, ").\n";
	my $downloaded_species = 0;
	my $skipped_species = 0;
	my $downloaded_assemblies = 0;
	
	SPECIES: foreach my $speciesDir (@subDir_content)
	{
		next if(($speciesDir eq '.') or ($speciesDir eq '..'));
		my $speciesDir_local = $subDir_local. '/' . $speciesDir;
		mkdir($speciesDir_local);
		die "Directory $speciesDir_local not existing, but it should (I just tried to mkdir it - do I have write permissions?" unless(-d $speciesDir_local);	

		my $speciesDir_with_latest = $speciesDir . '/latest_assembly_versions';
		
		unless($ftp->cwd(join('/', $ftp_root_genomes,  $DB, $subDir, $speciesDir_with_latest)))
		{
			$skipped_species++;
			next SPECIES;
		}	

		$downloaded_species++;
		my @speciesDir_latest_content = $ftp->ls();

		foreach my $assembly_version (@speciesDir_latest_content)
		{
			next if(($assembly_version eq '.') or ($assembly_version eq '..'));	

			$ftp->cwd(join('/', $ftp_root_genomes,  $DB, $subDir, $speciesDir_with_latest)) or die "Cannot change working directory ", $ftp->message;		
			$ftp->cwd($assembly_version) or die "Cannot change working directory (assembly_version) ", $ftp->message;		
			
			my @assembly_dir_contents = $ftp->ls();
	  
			my @genomic_fna_files = grep {$_ =~ /_genomic.fna.gz$/} grep {$_ !~ /(_cds_from_)|(_rna_from_g)/} @assembly_dir_contents;
			my @assembly_report_files = grep {$_ =~ /_assembly_report.txt$/} @assembly_dir_contents;
			die Dumper("Problem identifying files for download", \@genomic_fna_files, \@assembly_report_files, join('/', $ftp_root_genomes,  $DB, $subDir, $speciesDir_with_latest), $assembly_version) unless((scalar(@assembly_report_files) == 1) and (scalar(@genomic_fna_files) == 1));
			
			my $assemblyVersion_local = $speciesDir_local . '/'. $assembly_version;
			mkdir($assemblyVersion_local);
			die unless(-d $assemblyVersion_local);	
			
			my $cwd_before = getcwd();
			
			chdir($assemblyVersion_local) or die "Cannot chdir to $assemblyVersion_local";
			
			foreach my $file_to_transfer (@genomic_fna_files, @assembly_report_files)
			{
				$ftp->get($file_to_transfer) or die "Cannot transfer file", $ftp->message;		
			}		
			
			# die Dumper($assemblyVersion_local, \@genomic_fna_files);
			# print "Downloaded $assembly_version for $speciesDir\n";
			
			$downloaded_assemblies++;
			
			chdir($cwd_before) or die "Cannot chdir into $cwd_before";
		}
	}
	
	print "Summary for $subDir:\n";
	print "\tDownloaded species: $downloaded_species \n";
	print "\tSkipped species (most likely because there is no 'latest_assembly_versions' link in species directory): $skipped_species \n";
	print "\tDownloaded assemblies: $downloaded_assemblies\n";
	$DB_downloaded_assemblies += $downloaded_assemblies;
}	

print "\nDownload for $DB complete. Have $DB_downloaded_assemblies assemblies.\n";

# Print success message

my $suggested_command = 'NA';
if(-e abs_path($FindBin::Bin . '/annotateRefSeqSequencesWithUniqueTaxonIDs.pl'))
{
	$suggested_command = 'perl ' . abs_path($FindBin::Bin . '/annotateRefSeqSequencesWithUniqueTaxonIDs.pl') . qq( --refSeqDirectory $seqencesOutDirectory --taxonomyInDirectory $taxonomyOutDirectory --taxonomyOutDirectory DIR);
}

print qq(
Download successful - output directories:
- (sequences)  $seqencesOutDirectory
- (taxonomy)   $taxonomyOutDirectory

Suggested command for next step:

$suggested_command

);


sub print_help
{
	print qq(
downloadRefSeq.pl

  Downloads RefSeq/GenBank for MetaMap.
  
Usage:

  perl downloadRefSeq.pl --DB [refseq|genbank] --seqencesOutDirectory DIR --taxonomyOutDirectory DIR
  
Parameters:

  DB
      Name of NCBI database to be downloaded. Default: refseq.
  
  sequencesOutDirectory
  
      Output directory for sequences
	  
  taxonomyOutDirectory
      
	  Output directory for taxonomy
  
);
exit;
}
