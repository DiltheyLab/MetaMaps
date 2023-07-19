#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Net::FTP;
use File::Copy;
use File::Slurp;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use List::MoreUtils qw/all mesh/;
use Cwd qw/abs_path getcwd/;
$| = 1;

my $DB = 'refseq';
my $seqencesOutDirectory = '';
my $taxonomyOutDirectory;
my $targetBranches;
my $skipIncompleteGenomes = 0;
my $timeout = 360;

GetOptions (
	'DB:s' => \$DB, 
	'seqencesOutDirectory:s' => \$seqencesOutDirectory, 
	'taxonomyOutDirectory:s' => \$taxonomyOutDirectory, 
	'targetBranches:s' => \$targetBranches,
	'skipIncompleteGenomes:s' => \$skipIncompleteGenomes,
	'timeout:s' => \$timeout,
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
my $ftp;
initFTP(0);

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

#my @target_subdirs = qw/archaea bacteria fungi protozoa viral invertebrate plant unknown vertebrate_other vertebrate_mammalian/;
my @target_subdirs = qw/archaea bacteria fungi protozoa viral/;
if($targetBranches)
{
	$targetBranches =~ s/\s//g;
	@target_subdirs = split(/,/, $targetBranches);
}
#my @target_subdirs = qw/bacteria fungi protozoa viral invertebrate plant unknown vertebrate_other vertebrate_mammalian/;

my $ftp_root_genomes = '/genomes'; die unless(substr($ftp_root_genomes, 0, 1) eq '/');

my $DB_downloaded_assemblies = 0;

$ftp->cwd($ftp_root_genomes . '/' . $DB) or die "Cannot change working directory ", $ftp->message;

my $report_filename = $seqencesOutDirectory . '/report.txt';
open(my $report_fh, '>', $report_filename) or die "Could not open file '$report_filename' $!";

SUBDIR: foreach my $subDir (@target_subdirs)
{
    my $processed_entries = 1;
	my $subDir_local = $seqencesOutDirectory . '/' . $subDir;
	mkdir($subDir_local);
	die "Directory $subDir_local not existing, but it should (I just tried to mkdir it - do I have write permissions?"  unless(-d $subDir_local);
	
	$ftp->cwd($ftp_root_genomes . '/' . $DB . '/'. $subDir . '/') or die "Cannot change working directory ", $ftp->message;
	
	my @subDir_content = $ftp->ls();
	my $num_subdirs = scalar @subDir_content;
	print "Processing " . $num_subdirs . " entries";
	my @assembly_summary_files = grep {$_ =~ /assembly_summary.txt/} @subDir_content;

	unless(scalar(@assembly_summary_files) == 1)
	{
		say $report_fh Dumper("Could not identify assembly summary file in " . $ftp_root_genomes . '/' . $DB . '/'. $subDir . '/', \@assembly_summary_files);
		next;
	}
	
	my $assembly_summary_file = $assembly_summary_files[0];
	my $assembly_summary_file_local =  $subDir_local . '/' . $assembly_summary_file;

	ATTEMPT: for(my $attempt = 0; $attempt <= 100; $attempt++)
	{
		if (-e $assembly_summary_file_local)
		{
			print "\nAssembly summary file exists - skipping.\n";
			last ATTEMPT;
		}

		if($ftp->get($assembly_summary_file))
		{
			last ATTEMPT;
		}
		else
		{
			warn "Cannot transfer file " . $ftp->message;	
			if($attempt <= 1)
			{
				initFTP(0);
			}
			else
			{
				initFTP(1);
			}
			if($attempt >= 2)
			{
				say $report_fh "Attempt $attempt to get assembly_summary_file $assembly_summary_file failed";
				next SUBDIR;
			}
		}
	}

	if (! -e $assembly_summary_file_local) {
		move($assembly_summary_file, $assembly_summary_file_local) or die "Cannot move $assembly_summary_file to $assembly_summary_file_local - current cwd: " . getcwd();
	}

	my $assemblies_to_download = 0;
	my %assembly_dirs_by_species;
	die "File §assembly_summary_file_local not existing" unless(-e $assembly_summary_file_local);
	open(ASMSUM, '<', $assembly_summary_file_local) or die "Cannot open $assembly_summary_file_local";
	<ASMSUM>;
	my $asmsum_header = <ASMSUM>;
	chomp($asmsum_header);
	my @asmsum_header_fields = split(/\t/, $asmsum_header);
	while(<ASMSUM>)
	{
		chomp;
		next unless($_);
		my @line_fields = split(/\t/, $_, -1);
		die "Field number mismatch in file $assembly_summary_file_local - $#asmsum_header_fields / $#line_fields" unless(scalar(@asmsum_header_fields) == scalar(@line_fields));
		my %line = (mesh @asmsum_header_fields, @line_fields);
		my $ftp_path = $line{'ftp_path'};
		my $assembly_level = $line{'assembly_level'};
		my $organism_name = $line{'organism_name'};
		if (defined $ftp_path) {
			die unless(defined $assembly_level);
			next if($skipIncompleteGenomes and ($assembly_level ne 'Complete Genome'));
			(my $organism_name_safe = $organism_name) =~ s/\W/_/g;
			push(@{$assembly_dirs_by_species{$organism_name_safe}}, $ftp_path);
			$assemblies_to_download++;
		}
	}
	close(ASMSUM);
	my $total_dirs = 1;
	foreach my $check_species_name (keys %assembly_dirs_by_species) {
		my $check_subdir_species = $subDir_local. '/' . $check_species_name;
		if (-e $check_subdir_species) {
			my @check_dirs = read_dir($check_subdir_species);
			my $num_check_dirs = scalar @check_dirs;
			my $processed_species_subdirs = 0;
			my $all_ok = 0;
			my $checked_dirs = 1;
			foreach my $dir (@check_dirs) {
				my @check_files = read_dir($check_subdir_species . '/' . $dir);
				print "[".$total_dirs."/".(scalar(keys %assembly_dirs_by_species))."] : Checking [" . $checked_dirs . "/" . $num_check_dirs . "]" . $check_subdir_species . '/' . $dir . "\n";
				my @check_genomic_fna_files = grep {($_ =~ /_genomic.fna.gz$/) or ($_ =~ /_genomic.gff.gz$/) or ($_ =~ /_protein.faa.gz$/)} grep {$_ !~ /(_cds_from_)||(_rna_from_g)/} @check_files;
				my @check_assembly_report_files = grep {$_ =~ /_assembly_report.txt$/} @check_files;
				if (@check_genomic_fna_files and @check_assembly_report_files)
				{
				}
				else
				{
					$all_ok = 1;
				}
				$checked_dirs++;
			}

			if ($all_ok) {
				print "[".$total_dirs."/".(scalar(keys %assembly_dirs_by_species))."] : All local files OK\n";
				delete $assembly_dirs_by_species{$check_species_name};
				$assemblies_to_download--;
			}
		}
		$total_dirs++;
	}

	print "Now download genomes for ", scalar(keys %assembly_dirs_by_species), " $subDir species ($assemblies_to_download genomes - ", $DB, " - skip incomplete genomes: $skipIncompleteGenomes).\n";
		
	my $downloaded_species = 0;
	my $skipped_species = 0;
	my $downloaded_assemblies = 0;
	
	my $speciesI = 0;
	
	# init connection if checking local files takes ages (often get CLOSE_WAIT)
	initFTP(0);

	my $total_fileI = 0;
	SPECIES: foreach my $speciesName (keys %assembly_dirs_by_species)
	{
		$speciesI++;
				
		my $speciesDir_local = $subDir_local. '/' . $speciesName;
		mkdir($speciesDir_local);
		die "Directory $speciesDir_local not existing, but it should (I just tried to mkdir it - do I have write permissions?" unless(-d $speciesDir_local);	

		my $fileI = 0;
		foreach my $assembly_path_fullURL (@{$assembly_dirs_by_species{$speciesName}})
		{
			$fileI++;
			$total_fileI++;
			
			# last SPECIES  if($downloaded_assemblies > 100);
			(my $assembly_path_FTP = $assembly_path_fullURL) =~ s/ftp:\/\/ftp.ncbi.nlm.nih.gov//g;
			$assembly_path_FTP =~ s/https:\/\/ftp.ncbi.nlm.nih.gov//g;
			
			ATTEMPT: for(my $attempt = 0; $attempt <= 3; $attempt++)
			{
				if($attempt == 3)
				{
						close $report_fh;
						die "Cannot change working directory into assembly path (gave up - $assembly_path_fullURL) $assembly_path_FTP ", $ftp->message;
				}
				
				$ftp->cwd($assembly_path_FTP) or do {
						say $report_fh "Cannot change working directory into assembly path (attempt $attempt - $assembly_path_fullURL) $assembly_path_FTP " . $ftp->message;
						initFTP(0);
						next ATTEMPT;
				};
				last ATTEMPT;
			}

			my @assembly_dir_contents = $ftp->ls();
	  
			my @genomic_fna_files = grep {($_ =~ /_genomic.fna.gz$/) or ($_ =~ /_genomic.gff.gz$/) or ($_ =~ /_protein.faa.gz$/)} grep {$_ !~ /(_cds_from_)|(_rna_from_g)/} @assembly_dir_contents;
			my @assembly_report_files = grep {$_ =~ /_assembly_report.txt$/} @assembly_dir_contents;
						
			unless($assembly_path_fullURL =~ /.+\/(.+?)$/) {
				say $report_fh $speciesName . ": Can't parse assembly version from assembly path: $assembly_path_fullURL";
				next SPECIES;
			}
			my $assembly_version = $1;
			
			unless((scalar(@assembly_report_files) == 1) and (scalar(@genomic_fna_files) >= 1) and (scalar(@genomic_fna_files) <= 3)) {
				say $report_fh Dumper($speciesName . " [".$assembly_version."] Problem identifying files for download", \@genomic_fna_files, \@assembly_report_files, $assembly_path_fullURL, $assembly_version);
				next SPECIES;
			}
			
			my $assemblyVersion_local = $speciesDir_local . '/'. $assembly_version;
			mkdir($assemblyVersion_local);
			die unless(-d $assemblyVersion_local);
			
			my $cwd_before = getcwd();
			
			chdir($assemblyVersion_local) or die "Cannot chdir to $assemblyVersion_local";
			
			FTT: foreach my $file_to_transfer (@genomic_fna_files, @assembly_report_files)
			{
				print "\r\t Genome $total_fileI / $assemblies_to_download ; species $speciesI / ", scalar(keys %assembly_dirs_by_species), " $subDir ($speciesName) -- version $fileI / ", scalar(@{$assembly_dirs_by_species{$speciesName}}), ": GET $file_to_transfer                        ";
				ATTEMPT: for(my $attempt = 0; $attempt <= 100; $attempt++)
				{
					if (-e $file_to_transfer) {
						print "\nFile exists: [FTP] ", $ftp->size($file_to_transfer), " <=> ", (-s $file_to_transfer), " [LOCAL]\n";
						next FTT if($ftp->size($file_to_transfer) == (-s $file_to_transfer));
					}

					if($ftp->get($file_to_transfer))
					{
						last ATTEMPT;
					}	
					else
					{
						warn "Cannot transfer file " . $ftp->message;	
						if($attempt <= 1)
						{
							initFTP(0);
						}
						else
						{
							initFTP(1);
						}
						if($attempt >= 2)
						{
							say $report_fh "Attempt $attempt to get $file_to_transfer failed";
                                                        next FTT;
						}
					}
				}
			}
			
			# die Dumper($assemblyVersion_local, \@genomic_fna_files); 
			# print "Downloaded $assembly_version for $speciesDir\n";
			
			$downloaded_assemblies++;
			
			chdir($cwd_before) or die "Cannot chdir into $cwd_before";
		}
		
		$downloaded_species++;
	}
	
	
	print "\n";
	print "Summary for $subDir:\n";
	print "\tDownloaded species: $downloaded_species \n";
	# print "\tSkipped species (most likely because there is no 'latest_assembly_versions' link in species directory): $skipped_species \n";
	print "\tDownloaded assemblies: $downloaded_assemblies\n";
	$DB_downloaded_assemblies += $downloaded_assemblies;
}

# close reporting file
close $report_fh;

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
	  
  targetBranches  
  
	  Specification of target branches (comma-separated), e.g. archaea,bacteria,fungi
	  
  skipIncompleteGenomes
  
	  Skip incomplete genomes. Default: 0.
   
);
exit;
}

sub initFTP
{
	my $debug = shift;
	$ftp = Net::FTP->new($ftp_server, Debug => $debug, Timeout => $timeout) or die "Cannot connect to some.host.name: $@";
	$ftp->login("anonymous",'-anonymous@') or die "Cannot login ", $ftp->message;
	$ftp->binary();	
}
