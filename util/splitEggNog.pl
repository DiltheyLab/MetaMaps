use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Find;
use FindBin;
use lib "$FindBin::Bin/perlLib";
use List::Util qw/shuffle/;
use List::MoreUtils qw/all mesh/;
use Cwd qw/abs_path getcwd/;

my $input;
my $output;
my $action = '';
my $targetPerSplit = 100e6;
GetOptions (
	'input:s' => \$input, 
	'output:s' => \$output, 
	'action:s' => \$action, 
);

die "Please specify --input" unless($input);
die "Please specify --output" unless($output);
die "Please specify --action as split,submit,collect" unless(($action eq 'split') or ($action eq 'submit') or ($action eq 'collect'));

$input = abs_path($input);
$output = abs_path($output);

my $prefix_split = $output . '.split';
my $flag_split = $prefix_split . '.done';
if($action eq 'split')
{
	die "Input file split already? (flag file $flag_split present)" if(-e $flag_split);
	
	my @existing_split_input_files = glob($prefix_split . '.i.*');
	die Dumper(\@existing_split_input_files) if(@existing_split_input_files);
	foreach my $f (@existing_split_input_files)
	{
		unlink($f);
	}
	
	my $runningCharacters = 0;
	my $splitI = 0;
	my $runningOutputFh;
	
	my $openFh = sub {
		$splitI++;
		my $fn = $prefix_split . '.i.' . $splitI;
		if($runningOutputFh)
		{
			close($runningOutputFh) or die "Weird - couldn't close file";		
		}
		open($runningOutputFh, '>', $fn) or die "Cannot open $fn";
		$runningCharacters = 0;
	};
	
	$openFh->();
	open(I, '<', $input) or die "Cannot open $input";
	while(<I>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		if(substr($line, 0, 1) eq '>')
		{
			if($runningCharacters >= $targetPerSplit)
			{
				$openFh->();
			}	
			print {$runningOutputFh} $line, "\n";
		}		
		else
		{
			print {$runningOutputFh} $line, "\n";
			$runningCharacters += length($line);
		}
	}
	close(I);
	
	close($runningOutputFh) or die "Weird - couldn't close file";
	
	open(FLAG, '>', $flag_split) or die;
	print FLAG 1;
	close(FLAG);
	
	print "\nDone. Produced $splitI files.\n\n";
}
elsif($action eq 'submit')
{
	die "Input file not split yet? (missing flag file $flag_split)" unless(-e $flag_split);
	
	my @existing_split_input_files = glob($prefix_split . '.i.*');

	print "Found ", scalar(@existing_split_input_files), " files.\n";
	
	my $fn_submitall = $prefix_split . '.submit.all';
	open(SUBMITALL, '>', $fn_submitall) or die "Cannot open $fn_submitall";
	foreach my $splitFile (@existing_split_input_files)
	{
		die unless($splitFile =~ /\.i\.(\d+)$/);
		my $jobN = $1;
		(my $submitFile = $splitFile) =~ s/\.i\./.submit./;
		(my $outputFile = $splitFile) =~ s/\.i\./.o./;
		my $outputFile_OKflag = $outputFile . '.done';
		if(-e $outputFile_OKflag)
		{
			unlink($outputFile_OKflag) or die "Cannot unlink $outputFile_OKflag";
		}	
		open(SUBMIT, '>', $submitFile) or die;
print SUBMIT qq(#!/bin/bash
#\$ -q phillippy.q
#\$ -l mem_free=25G
#\$ -N DIAMOND${jobN}
cd /data/projects/phillippy/projects/MetaMap/downloads/eggNOG/eggnog-mapper-1.0.3
python emapper.py -i $splitFile --output $outputFile -m diamond --usemem --cpu 8 && echo 1 > $outputFile_OKflag
);		
		close(SUBMIT);
		
		print SUBMITALL "qsub $submitFile \n";
	}	
	close(SUBMITALL);
	
	print "\n\nDone. Execute commands in $fn_submitall to submit.\n\n";
}
elsif($action eq 'collect')
{
	die "Input file split already? (flag file $flag_split present)" unless(-e $flag_split);
	
	my @existing_split_input_files = glob($prefix_split . '.i.*');
	
	foreach my $splitFile (@existing_split_input_files)
	{
		(my $outputFile = $splitFile) =~ s/\.i\./.o./;
		my $outputFile_OKflag = $outputFile . '.done';
		my $annotations_file = $outputFile . '.emapper.annotations';		
		warn "Output OK flag file $outputFile_OKflag not present" unless(-e $outputFile_OKflag);
		die "File $annotations_file not present" unless(-e $annotations_file);
	}
	
	open(OUTPUT, '>', $output) or die "Cannot open $output";
	print OUTPUT join("\t", qw/ProteinID GO_terms KEGG_KOs BiGG_reactions OGs COG_cat/), "\n";
	
	
	foreach my $splitFile (@existing_split_input_files)
	{
		(my $outputFile = $splitFile) =~ s/\.i\./.o./;
		my $annotations_file = $outputFile . '.emapper.annotations';
		open(CHUNK, '<', $annotations_file) or die "Cannot open $annotations_file";
		<CHUNK>;
		<CHUNK>;
		<CHUNK>;
		my $headerLine = <CHUNK>;
		chomp($headerLine);
		my @header_fields = split(/\t/, $headerLine);
		my %inHeader = map {$_ => 1} @header_fields;
		die unless($inHeader{'#query_name'});
		die unless($inHeader{'GO_terms'});
		die unless($inHeader{'KEGG_KOs'});
		die unless($inHeader{'BiGG_reactions'});
		die unless($inHeader{'OGs'});
		die unless($inHeader{'COG cat'});
		while(<CHUNK>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			next if(substr($line, 0, 1) eq '#');
			my @line_fields = split(/\t/, $line, -1);
			my %line = (mesh @header_fields, @line_fields);
			print OUTPUT join("\t", map {die "Line field $_ undefined\n$line $annotations_file" unless(defined $line{$_}); $line{$_}} ('#query_name', 'GO_terms', 'KEGG_KOs', 'BiGG_reactions', 'OGs', 'COG cat')), "\n";
			# die "Weird COG value $line{'COG cat'}" unless(
				# (not $line{'COG cat'}) or 
				# ($line{'COG cat'} =~ /^[\w]$/)
			# );
 		}
		close(CHUNK);
	}
	
	close(OUTPUT);

	print "\nGenerated file $output.\n\n";
}
else
{
	die;
}	
