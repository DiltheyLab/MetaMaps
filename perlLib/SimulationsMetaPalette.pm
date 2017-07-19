package SimulationsMetaPalette;

use strict;
use Data::Dumper;
use Cwd qw/abs_path/;

use taxTree;

sub doMetaPalette
{
	my $jobDir = shift;
	my $dbDir = shift;
	my $reads_fastq = shift;
	my $metaPalette_installation_dir = shift;
	my $jellyfish_2_bin = shift;
	my $fullTaxonomy = shift;
	
	die unless(defined $fullTaxonomy);
	
	my $metapalette_dir = $jobDir . '/metapalette/';

	if(-e $metapalette_dir)
	{
		system("rm -rf $metapalette_dir") and die "Couldn't delete $metapalette_dir";
	}
	unless(-e $metapalette_dir)
	{
		mkdir($metapalette_dir) or die "Cannot open $metapalette_dir";
	}
		
	$metapalette_dir = abs_path($metapalette_dir);
	die unless(-d $metapalette_dir); 
	my $jobDir_abs = abs_path($jobDir);
	
	$reads_fastq = abs_path($reads_fastq);		
	die "Simulated reads $reads_fastq not found" unless(-e $reads_fastq);
	
	my $metaPalette_src_dir = $metaPalette_installation_dir . '/src/Python';
	
	my $metaPalette_bacteria_DB_dir = $metaPalette_installation_dir . '/Bacteria';
	die "Expected bacteria DB dir $metaPalette_bacteria_DB_dir not found" unless(-e $metaPalette_bacteria_DB_dir);
	
	my $metaPalette_queryPerSeq = $metaPalette_installation_dir . '/src/QueryPerSeq/query_per_sequence';
	die "query_per_sequence not found" unless(-e $metaPalette_queryPerSeq);
	
	my $cwd_before = getcwd();
	
	chdir($metaPalette_src_dir) or die "Cannot chdir $metaPalette_src_dir";
	  
	# print "module load libs/hdf5/1.8.12\n";
	# print "module load libs/hdf5/1.8.12\n"; 

	
	my $cmd_first_quartile = qq(perl ${FindBin::Bin}/firstQuartileScore.pl $reads_fastq);
	my $return_first_quartile = `$cmd_first_quartile`;
	chomp($return_first_quartile);
	die Dumper("Weird output from $cmd_first_quartile", length($return_first_quartile), $return_first_quartile) unless(length($return_first_quartile) == 1);
	
	my $cmd_metaPalette = qq(eval 'module load libs/hdf5/1.8.12'; eval 'module load libs/scipy/0.18.1'; python Classify.py -d $metaPalette_bacteria_DB_dir -o $metapalette_dir -i $reads_fastq -Q ). quotemeta($return_first_quartile) . qq( -k sensitive -j $jellyfish_2_bin -q $metaPalette_queryPerSeq -t 16 -n);
	
	print $cmd_metaPalette, "\n";
	system($cmd_metaPalette) and die "Couldn't execute metapalette command $cmd_metaPalette";
	
	create_compatible_file_from_metapalette( 
		$jobDir_abs . '/results_metapalette.txt',
		$metapalette_dir . '/combined_reads.fastq.profile',
		$fullTaxonomy,
	);
	
	chdir($cwd_before) or die "Cannot chdir $cwd_before";
	
}

sub create_compatible_file_from_metapalette
{
	my $output_fn = shift;
	my $output_metapalette_fn = shift;
	my $fullTaxonomy = shift;
	die unless(defined $fullTaxonomy);
	
	my %S_byLevel;
		
	my $taxonomy_metaPalette = taxTree::readTaxonomy($fullTaxonomy);
	my $merged_href = taxTree::readMerged($fullTaxonomy);
	
	open(MP, '<', $output_metapalette_fn) or die "Cannot open $output_metapalette_fn";
	while(<MP>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		next if(substr($line, 0, 1) eq '#');
		next if(substr($line, 0, 1) eq '@');
		my @line_fields = split(/\t/, $line);
		die unless(scalar(@line_fields) == 5);
		my $taxonID = $line_fields[0];
		my $level = $line_fields[1];
		my $percentage = $line_fields[4];
		
		my $currentTaxonID = taxTree::findCurrentNodeID($taxonomy_metaPalette, $merged_href, $taxonID);
		
		die "Unknown taxonomy ID $taxonID" unless(exists $taxonomy_metaPalette->{$currentTaxonID});
		warn "Weird rank for taxonomy ID $currentTaxonID -- is $level expect $taxonomy_metaPalette->{$currentTaxonID}{rank}" unless($taxonomy_metaPalette->{$currentTaxonID}{rank} eq $level);
		
		my $realRank = $taxonomy_metaPalette->{$currentTaxonID}{rank};
		if($realRank)
		{
			$S_byLevel{$realRank}{$currentTaxonID} += ($percentage / 100);
			
			# my $name = taxTree::taxon_id_get_name($currentTaxonID, $taxonomy_metaPalette);
			# $S_byLevel{$realRank}{$name} += ($percentage / 100);
		}
	}
	close(MP);
	
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
	print OUTPUT join("\t", qw/AnalysisLevel taxonID Name Absolute PotFrequency/), "\n";	
	foreach my $l (keys %S_byLevel)
	{
		my $S_level = 0;
		foreach my $taxonID (keys %{$S_byLevel{$l}})
		{
			$S_level += $S_byLevel{$l}{$taxonID};
		}
		die unless(($S_level >= 0) and ($S_level <= 1));
		my $S_missing = 1 - $S_level;
		
		$S_byLevel{$l}{'Unclassified'} = $S_missing;
		
		foreach my $taxonID (keys %{$S_byLevel{$l}})
		{
			my $taxonID_for_print = $taxonID;
			my $name;
			if($taxonID eq 'Unclassified')
			{
				$name = 'Unclassified';
				$taxonID_for_print = 0;
			}
			else
			{
				$name = taxTree::taxon_id_get_name($taxonID, $taxonomy_metaPalette);
			}
			
			print OUTPUT join("\t", $l, $taxonID, $name, 0, $S_byLevel{$l}{$taxonID}), "\n";
		}
	}
	close(OUTPUT);
}

1;