package SimulationsKraken;

use strict;
use Data::Dumper;
use Cwd qw/abs_path getcwd/;
use List::Util qw/min max all/;

use taxTree;
use validation;

sub getKrakenBinPrefix
{
	return qq(/data/projects/phillippy/software/kraken-0.10.5-beta/bin/kraken);
}

sub getKraken2BinPrefix
{
	return qq(/data/projects/phillippy/software/kraken2-2.0.7-beta/bin/kraken2);
}

sub getBrackenDir
{
	return qq(/data/projects/phillippy/software/Bracken/);
}

sub getKrakenDBTemplate
{
	return '/data/projects/phillippy/projects/mashsim/src/krakenDBTemplate2/'; # make sure this is current!
}

sub getKraken2DBTemplate
{
	return '/data/projects/phillippy/projects/MetaMap/kraken2-dbtemplate/'; # make sure this is current!
}

sub getCentrifugeDir
{
	return qq(/data/projects/phillippy/software/centrifuge-1.0.4-beta);
}

sub getMasterTaxonomy
{
	return qq(/data/projects/phillippy/software/centrifuge-1.0.4-beta);
}


sub translateMetaMapToCentrifuge
{
	my $outputDir = shift;
	my $MetaMapDBDir = shift;
	my $centrifugeBinDir = shift;
	
	my $dbDir_abs = abs_path($MetaMapDBDir);
		

	
	my $MetaMapDBDir_abs = abs_path($MetaMapDBDir);
	
	unless(-e $outputDir)
	{
		mkdir($outputDir) or die "Cannot open $outputDir";
	}
	
	my $pre_chdir_cwd = getcwd();
	
	chdir($outputDir) or die; 
	
	my $fasta_for_mapping = abs_path($dbDir_abs . '/DB.fa');	
	
	my @files = map {'DB.' . $_ . '.cf'} 1 .. 3;
	if(all {-e $_} @files)
	{	
		my $max_time_since_centrifuge_index_change = max map {-M $_} @files;
		if($max_time_since_centrifuge_index_change < (-M $fasta_for_mapping))
		{
			warn "Centrifuge DB already present and up to date!";		
			chdir($pre_chdir_cwd) or die;		
			return;			
		}
		else
		{
			# die Dumper([$files[0], -M $files[0]], [$fasta_for_mapping, -M $fasta_for_mapping]);
		}
	}	
	
	if(-e $outputDir)
	{
		system("rm -rf $outputDir") and die "Cannot delete $outputDir";
		mkdir($outputDir) or die "Cannot open $outputDir";		
	}
	

	
	# if(-e 'DB')
	# {
		# system('rm -rf DB') and die "Cannot rm";
	# }
	# mkdir('DB') or die "Cannot mkdir DB";
	
	die "Required DB file $fasta_for_mapping not existing ($MetaMapDBDir)" unless(-e $fasta_for_mapping);	


	my $cmd_convert = qq(perl ${FindBin::Bin}/util/conversionTableForCentrifuge.pl --DB $MetaMapDBDir_abs );
	system($cmd_convert) and die "Could not execute command: $cmd_convert";
	
	die "Converted mashmap DB (mashmap -> centrifuge) missing!" unless(-e "${fasta_for_mapping}.centrifugeTranslation");
	my $taxonomy_names = $fasta_for_mapping . '.centrifugeTranslation.names.dmp';
	my $taxonomy_nodes = $fasta_for_mapping . '.centrifugeTranslation.nodes.dmp';
	die unless((-e $taxonomy_nodes) and (-e $taxonomy_names));
		
	# system("mv ${fasta_for_mapping}.centrifugeTranslation .") and die "Cannot move ${fasta_for_mapping}.centrifuge";
	 
	my $cmd_build_I = qq(${centrifugeBinDir}/centrifuge-build -p 32 --conversion-table ${fasta_for_mapping}.centrifugeTranslation --taxonomy-tree $taxonomy_nodes --name-table $taxonomy_names $fasta_for_mapping DB);
	system($cmd_build_I) and die "Could not execute command: $cmd_build_I";

	mkdir('taxonomy') or die "Cannot mkdir DB/taxonomy from " . getcwd();
	system("mv $taxonomy_names taxonomy/names.dmp") and die;
	system("mv $taxonomy_nodes taxonomy/nodes.dmp") and die;
	
	chdir($pre_chdir_cwd) or die;
}

sub translateMetaMapToKraken2
{
	my $kraken2_dir = shift;
	my $MetaMapDBDir = shift;
	my $kraken2DBTemplate = shift;
	my $kraken2_binPrefix = shift;
	
	my $dbDir_abs = abs_path($MetaMapDBDir);
	
	my $fasta_for_mapping = abs_path($dbDir_abs . '/DB.fa');	
	die "Required DB file $fasta_for_mapping not existing ($MetaMapDBDir)" unless(-e $fasta_for_mapping);	
	
	my $kraken2_main_DB_file = $kraken2_dir . '/DB/hash.k2d';
	if((-e $kraken2_main_DB_file) and ((-M $kraken2_main_DB_file) < (-M $fasta_for_mapping)))
	{
		warn "Kraken2 DB already present and up to date!";		
		return;			
	}

	if(-e $kraken2_dir)
	{
		system("rm -rf $kraken2_dir") and die "Cannot delete $kraken2_dir";
	}
	
	unless(-e $kraken2_dir)
	{
		mkdir($kraken2_dir) or die "Cannot open $kraken2_dir";
	}
	
	my $pre_chdir_cwd = getcwd();
	
	chdir($kraken2_dir) or die; 

	
	if(-e 'DB')
	{
		system('rm -rf DB') and die "Cannot rm";
	}
	
			
	
	my $cmd_copy_DB = qq(cp -r $kraken2DBTemplate DB);
	system($cmd_copy_DB) and die "Cannot cp $kraken2DBTemplate";
	die "DB missing" unless(-d 'DB');
	$kraken2DBTemplate = abs_path($kraken2DBTemplate);

	my $cmd_convert = qq(perl ${FindBin::Bin}/util/translateMashmapDBToKraken.pl --input $fasta_for_mapping --taxonomyDir ${dbDir_abs}/taxonomy --krakenTemplate_taxonomy ${kraken2DBTemplate}/taxonomy/);
	system($cmd_convert) and die "Could not execute command: $cmd_convert";
	die "Converted mashmap DB (mashmap -> kraken) missing!" unless(-e "${fasta_for_mapping}.kraken");
	
	system("mv ${fasta_for_mapping}.kraken .") and die "Cannot move ${fasta_for_mapping}.kraken";
	 
	my $cmd_build_II = qq(export PATH=\$PATH:/data/projects/phillippy/software/blast/ncbi-blast-2.6.0+/bin/; /usr/bin/time -v ${kraken2_binPrefix}-build --add-to-library DB.fa.kraken --threads 32 --db DB &> output_build_II.txt);
	system($cmd_build_II) and die "Could not execute command: $cmd_build_II";
	
	my $cmd_build_III = qq(/usr/bin/time -v ${kraken2_binPrefix}-build --build -threads 32 --db DB &> output_build_III.txt);
	system($cmd_build_III) and die "Could not execute command: $cmd_build_III";
	
	chdir($pre_chdir_cwd) or die;
}

sub translateMetaMapToKraken
{
	my $kraken_dir = shift;
	my $MetaMapDBDir = shift;
	my $krakenDBTemplate = shift;
	my $kraken_binPrefix = shift;
	my $Bracken_dir = shift;
	
	my $dbDir_abs = abs_path($MetaMapDBDir);

	my $fasta_for_mapping = abs_path($dbDir_abs . '/DB.fa');	
	die "Required DB file $fasta_for_mapping not existing ($MetaMapDBDir)" unless(-e $fasta_for_mapping);		
	
	my $kraken_present = 0;
	my $kraken_main_DB_file = $kraken_dir . '/DB/database.kdb';
	if((-e $kraken_main_DB_file) and ((-M $kraken_main_DB_file) < (-M $fasta_for_mapping)))
	{
		$kraken_present = 1;
		
	}

	my $bracken_present = 0;	
	my $bracken_main_DB_file = $kraken_dir . '/database75mers.kraken_cnts.bracken';
	if((-e $bracken_main_DB_file) and ((-M $bracken_main_DB_file) < (-M $kraken_main_DB_file)))
	{
		$bracken_present = 1;	
	}
	
	if($kraken_present and $bracken_present)
	{
		warn "Kraken/Bracken DB already present and up to date!";		
		return;		
	}
	
				
	if(-e $kraken_dir)
	{
		system("rm -rf $kraken_dir") and die "Cannot delete $kraken_dir";
	}
	
	unless(-e $kraken_dir)
	{
		mkdir($kraken_dir) or die "Cannot open $kraken_dir";
	}
	
	my $pre_chdir_cwd = getcwd();
	
	chdir($kraken_dir) or die; 

	
	if(-e 'DB')
	{
		system('rm -rf DB') and die "Cannot rm";
	}

	my $cmd_copy_DB = qq(cp -r $krakenDBTemplate DB);
	system($cmd_copy_DB) and die "Cannot cp $krakenDBTemplate";
	die "DB missing" unless(-d 'DB');

	$krakenDBTemplate = abs_path($krakenDBTemplate);
	
	my $cmd_convert = qq(perl ${FindBin::Bin}/util/translateMashmapDBToKraken.pl --input $fasta_for_mapping --taxonomyDir ${dbDir_abs}/taxonomy --krakenTemplate_taxonomy ${krakenDBTemplate}/taxonomy/);
	system($cmd_convert) and die "Could not execute command: $cmd_convert";
	die "Converted mashmap DB (mashmap -> kraken) missing!" unless(-e "${fasta_for_mapping}.kraken");
	
	system("mv ${fasta_for_mapping}.kraken .") and die "Cannot move ${fasta_for_mapping}.kraken";
	 
	my $cmd_build_II = qq(/usr/bin/time -v ${kraken_binPrefix}-build --add-to-library DB.fa.kraken --threads 32 --db DB &> output_build_II.txt);
	system($cmd_build_II) and die "Could not execute command: $cmd_build_II";
	
	my $cmd_build_III = qq(export PATH=/data/projects/phillippy/software/jellyfish-1.1.11/bin:\$PATH; /usr/bin/time -v ${kraken_binPrefix}-build --build --threads 32 --db DB &> output_build_III.txt);
	system($cmd_build_III) and die "Could not execute command: $cmd_build_III";

	
	my $cmd_Bracken_selfSimilarity = qq(bash -c '/usr/bin/time -v ${kraken_binPrefix} --db DB --fasta-input --threads=32 <( find -L DB/library \\( -name "*.fna"  -o -name "*.fa" -o -name "*.fasta" \\) -exec cat {} + ) > database_kraken');
	system($cmd_Bracken_selfSimilarity) and die "Could not execute command: $cmd_Bracken_selfSimilarity";

	my $cmd_Bracken_countkMers = qq(/usr/bin/time -v perl ${Bracken_dir}/count-kmer-abundances.pl --db=DB --read-length=2000 --threads=32 database_kraken > database75mers.kraken_cnts);
	system($cmd_Bracken_countkMers) and die "Could not execute command: $cmd_Bracken_countkMers";
	
	my $cmd_Bracken_kMerDist = qq(/usr/bin/time -v python ${Bracken_dir}/generate_kmer_distribution.py -i database75mers.kraken_cnts -o database75mers.kraken_cnts.bracken);
	system($cmd_Bracken_kMerDist) and die "Could not execute command: $cmd_Bracken_kMerDist";	
	
	chdir($pre_chdir_cwd) or die;
}

sub doKrakenOnExistingDB
{
	my $kraken_dir = shift;
	my $simulatedReads = shift;
	my $outputDir = shift;
	my $kraken_binPrefix = shift;
	my $Bracken_dir = shift;
	my $taxonID_original_2_contigs_href = shift;
	
	die unless(defined $taxonID_original_2_contigs_href);
	my $pre_chdir_cwd = getcwd();
	 
	chdir($kraken_dir) or die "Cannot chdir into $kraken_dir";  
	
	my $cmd_classify = qq(/usr/bin/time -v ${kraken_binPrefix} --preload --db DB $simulatedReads 1> $outputDir/reads_classified 2> $outputDir/kraken_resources);
	system($cmd_classify) and die "Could not execute command: $cmd_classify"; # todo
	
	my $cmd_report = qq(/usr/bin/time -v ${kraken_binPrefix}-report --db DB $outputDir/reads_classified 1> $outputDir/reads_classified_report 2> $outputDir/kraken_report_resources);
	system($cmd_report) and die "Could not execute command: $cmd_report"; # todo 
	
	foreach my $L (qw/S G F/)
	{ 
		my $cmd_Bracken_estAbundance = qq(/usr/bin/time -v python ${Bracken_dir}/est_abundance.py -i $outputDir/reads_classified_report -k database75mers.kraken_cnts.bracken -l $L -o $outputDir/reads_classified_report_bracken_${L} 2> $outputDir/bracken_resources_$L);
		system($cmd_Bracken_estAbundance) and die "Could not execute command: $cmd_Bracken_estAbundance"; # todo
	} 
	
	create_compatible_file_from_kraken(
		$outputDir . '/results_kraken.txt',
		'DB/taxonomy',
		$outputDir.'/reads_classified_report',
		$outputDir.'/reads_classified',	
		$taxonID_original_2_contigs_href
	);

	create_compatible_reads_file_from_kraken( 
		$outputDir . '/results_kraken.txt.reads2Taxon',
		'DB/taxonomy',
		$outputDir.'/reads_classified',
	);
		
	create_compatible_file_from_kraken_bracken(
		$outputDir . '/results_bracken.txt',
		'DB/taxonomy',
		$outputDir.'/reads_classified_report',
		$outputDir.'/reads_classified_report_bracken_S',
		$outputDir.'/reads_classified_report_bracken_G',
		$outputDir.'/reads_classified_report_bracken_F');
		
	chdir($pre_chdir_cwd) or die;			
}

sub doKraken2OnExistingDB
{
	my $kraken2_dir = shift;
	my $simulatedReads = shift;
	my $outputDir = shift;
	my $kraken2_binPrefix = shift;
	my $taxonID_original_2_contigs_href = shift;
	
	die unless(defined $taxonID_original_2_contigs_href);
	my $pre_chdir_cwd = getcwd();
	 
	chdir($kraken2_dir) or die "Cannot chdir into $kraken2_dir";  
	
	my $cmd_classify = qq(/usr/bin/time -v ${kraken2_binPrefix} --db DB --report $outputDir/reads_classified_report $simulatedReads 1> $outputDir/reads_classified 2> $outputDir/kraken_resources);
	system($cmd_classify) and die "Could not execute command: $cmd_classify"; # todo
	
	#my $cmd_report = qq(/usr/bin/time -v ${kraken2_binPrefix}-report --db DB $outputDir/reads_classified 1> $outputDir/reads_classified_report 2> $outputDir/kraken_report_resources);
	#system($cmd_report) and die "Could not execute command: $cmd_report"; # todo 

	create_compatible_file_from_kraken(
		$outputDir . '/results_kraken2.txt',
		'DB/taxonomy',
		$outputDir.'/reads_classified_report',
		$outputDir.'/reads_classified',	
		$taxonID_original_2_contigs_href
	);

	create_compatible_reads_file_from_kraken( 
		$outputDir . '/results_kraken2.txt.reads2Taxon',
		'DB/taxonomy',
		$outputDir.'/reads_classified',
	);

	chdir($pre_chdir_cwd) or die;			
}

sub doCentrifuge
{
	my $jobDir = shift;
	my $dbDir = shift;
	my $reads_fastq = shift;
	my $centrifugeBinDir = shift;
		
	my $simulatedReads = abs_path($reads_fastq);
	die unless(-e $simulatedReads);
	
	my $centrifuge_dir = $jobDir . '/centrifuge/';
	my $jobDir_abs = abs_path($jobDir);
	
	unless(-e $centrifuge_dir)
	{
		mkdir($centrifuge_dir) or "Cannot mkdir $centrifuge_dir";
	}
	
	my $centrifuge_DB_dir = $centrifuge_dir . '/DB/';
	unless(-e $centrifuge_DB_dir)
	{
		mkdir($centrifuge_DB_dir) or "Cannot mkdir $centrifuge_DB_dir";
	}	
	
	my %taxonID_original_2_contigs;
	my %contigLength;
	Util::read_taxonIDs_and_contigs($dbDir, \%taxonID_original_2_contigs, \%contigLength);
	
	SimulationsKraken::translateMetaMapToCentrifuge (
		$centrifuge_DB_dir,
		$dbDir, 
		$centrifugeBinDir,
	);

	SimulationsKraken::doCentrifugeOnExistingDB (
		$centrifuge_DB_dir,
		$simulatedReads,
		$centrifuge_dir,
		$centrifugeBinDir,
		\%taxonID_original_2_contigs
	);
}

sub doCentrifugeOnExistingDB
{
	my $centrifugeDBDir = shift;
	my $simulatedReads = shift;
	my $outputDir = shift;
	my $centrifugeBinDir = shift;
	my $taxonID_original_2_contigs_href = shift;
	
	die unless(defined $taxonID_original_2_contigs_href);
	my $pre_chdir_cwd = getcwd();
	 
	# chdir($kraken2_dir) or die "Cannot chdir into $kraken2_dir";  
	
	die "We have a directory problem - DB dir $centrifugeDBDir not present from " . getcwd() unless(-d $centrifugeDBDir);
	die "We have a directory problem - DB as specified in $centrifugeDBDir not present" unless(-e $centrifugeDBDir . '/DB.1.cf');
	my $cmd_classify = qq(/usr/bin/time -v ${centrifugeBinDir}/centrifuge  -x ${centrifugeDBDir}/DB -U $simulatedReads -S $outputDir/reads_classified --report-file $outputDir/reads_classified_report 2> $outputDir/centrifuge_resources);
	
	system($cmd_classify) and die "Could not execute command: $cmd_classify"; # todo

	
	my $cmd_report = qq(/usr/bin/time -v ${centrifugeBinDir}/centrifuge-kreport  -x ${centrifugeDBDir}/DB $outputDir/reads_classified > $outputDir/reads_classified_kreport 2> $outputDir/centrifuge_report_resources);
	system($cmd_report) and die "Could not execute command: $cmd_report"; # todo
	
	
	#my $cmd_report = qq(/usr/bin/time -v ${kraken2_binPrefix}-report --db DB $outputDir/reads_classified 1> $outputDir/reads_classified_report 2> $outputDir/kraken_report_resources);
	#system($cmd_report) and die "Could not execute command: $cmd_report"; # todo 
	
	create_compatible_composition_file_from_centrifuge(
		$outputDir . '/results_centrifuge.txt',
		"${centrifugeDBDir}/DB/taxonomy",
		$outputDir.'/reads_classified_kreport',
		$outputDir.'/reads_classified',	
		$taxonID_original_2_contigs_href
	);

	create_compatible_reads_file_from_centrifuge( 
		$outputDir . '/results_centrifuge.txt.reads2Taxon', 
		"${centrifugeDBDir}/DB/taxonomy",
		$outputDir.'/reads_classified',
	);

	chdir($pre_chdir_cwd) or die;			
}



sub doKraken
{
	my $jobDir = shift;
	my $dbDir = shift;
	my $reads_fastq = shift;
	
	my $krakenDBTemplate = shift;
	my $kraken_binPrefix = shift;
	my $Bracken_dir = shift;
	
	die unless(defined $Bracken_dir);
	

	my %taxonID_original_2_contigs;
	my %contigLength;
	Util::read_taxonIDs_and_contigs($dbDir, \%taxonID_original_2_contigs, \%contigLength);

	my $kraken_dir = $jobDir . '/kraken/';
	my $jobDir_abs = abs_path($jobDir);
	
	unless(-e $kraken_dir)
	{
		mkdir($kraken_dir) or "Cannot mkdir $kraken_dir";
	}
	
	translateMetaMapToKraken($kraken_dir, $dbDir, $krakenDBTemplate, $kraken_binPrefix, $Bracken_dir); # todo
	
	my $simulatedReads = abs_path($reads_fastq);
	die unless(-e $simulatedReads);
	
	my $outputPrefix = '';
	doKrakenOnExistingDB($kraken_dir, $simulatedReads, $jobDir_abs, $kraken_binPrefix, $Bracken_dir, \%taxonID_original_2_contigs);
	
}

sub doKraken2
{
	my $jobDir = shift;
	my $dbDir = shift;
	my $reads_fastq = shift;
	
	my $kraken2DBTemplate = shift;
	my $kraken2_binPrefix = shift;
		

	my %taxonID_original_2_contigs;
	my %contigLength;
	Util::read_taxonIDs_and_contigs($dbDir, \%taxonID_original_2_contigs, \%contigLength);

	my $kraken2_dir = $jobDir . '/kraken2/';
	my $jobDir_abs = abs_path($jobDir);
	
	unless(-e $kraken2_dir)
	{
		mkdir($kraken2_dir) or "Cannot mkdir $kraken2_dir";
	}
	translateMetaMapToKraken2($kraken2_dir, $dbDir, $kraken2DBTemplate, $kraken2_binPrefix); # todo
	
	my $simulatedReads = abs_path($reads_fastq);
	die unless(-e $simulatedReads);
	
	my $outputPrefix = '';
	doKraken2OnExistingDB($kraken2_dir, $simulatedReads, $jobDir_abs, $kraken2_binPrefix, \%taxonID_original_2_contigs);
	
}

sub create_compatible_file_from_kraken
{
	my $output_fn = shift;
	my $taxonomy_kraken_dir = shift;
	my $f_K = shift;
	my $f_reads = shift;
	my $create_compatible_file_from_kraken = shift;	
	die unless(defined $create_compatible_file_from_kraken);
	
	my $taxonomy_kraken = taxTree::readTaxonomy($taxonomy_kraken_dir);
	
	my $output_fn_2 = $output_fn . '.ignoreUnclassified';
	my %S_byLevel;
	my $n_unclassified;
	my $n_root;
	open(KRAKEN, '<', $f_K) or die "Cannot open $f_K";
	while(<KRAKEN>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		if($f[5] eq 'unclassified')
		{
			die if(defined $n_unclassified);
			$n_unclassified = $f[1];
			next;
		}
		elsif($f[5] eq 'root')
		{
			die if(defined $n_root);
			$n_root = $f[1];
			next;
		}
	}
	
	$n_unclassified = 0 unless(defined $n_unclassified);
	my $n_total_reads = $n_unclassified + $n_root;
	
	my $n_unclassified_check = 0;
	
	my %_getLightning_cache;
	my $getLightning = sub {
		my $taxonID = shift;
		if(exists $_getLightning_cache{$taxonID})
		{
			return $_getLightning_cache{$taxonID};
		}
		else
		{
			my $lightning = validation::getAllRanksForTaxon_withUnclassified($taxonomy_kraken, $taxonID, $create_compatible_file_from_kraken);
			$_getLightning_cache{$taxonID} = $lightning;
			return $lightning;
		}
	};
	
	if($n_total_reads == 0)
	{
		die "Weird - $n_total_reads reads? Output $output_fn, kraken output $f_K";
	}
	
	my @evaluateAccuracyAtLevels = validation::getEvaluationLevels();
	
	my %reads_at_levels;
	open(KRAKEN, '<', $f_reads) or die "Cannot open $f_reads";
	while(<KRAKEN>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		my $classified = $f[0];
		my $readID = $f[1];
		my $taxonID = $f[2];
		die "Weird classification symbol in $f_reads: '$classified' (line $. of $f_reads in ".getcwd(). ")" unless(($classified eq 'C') or ($classified eq 'U'));
		if($classified eq 'C')
		{
			my $lightning = $getLightning->($taxonID);
			$reads_at_levels{'definedAndHypotheticalGenomes'}{$taxonID}++;
			RANK: foreach my $rank (@evaluateAccuracyAtLevels)
			{
				die unless(defined $lightning->{$rank});
				$reads_at_levels{$rank}{$lightning->{$rank}}++;
			}			
		}
		else
		{
			$n_unclassified_check++;
		}
	}
	close(KRAKEN);
	die "Inconsistency w.r.t. unclassified reads -- $n_unclassified_check vs $n_unclassified" unless($n_unclassified_check == $n_unclassified);
	
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
	open(OUTPUT2, '>', $output_fn_2) or die "Cannot open $output_fn_2";
	print OUTPUT join("\t", qw/AnalysisLevel ID Name Absolute PotFrequency/), "\n";		
	print OUTPUT2 join("\t", qw/AnalysisLevel ID Name Absolute PotFrequency/), "\n";		
	
	foreach my $level ('definedAndHypotheticalGenomes', @evaluateAccuracyAtLevels)
	{
		$reads_at_levels{$level}{"Unclassified"} = 0 if(not exists $reads_at_levels{$level}{"Unclassified"});
		# $reads_at_levels{$level}{"Undefined"} = 0 if(not exists $reads_at_levels{$level}{"Undefined"});
		
		$reads_at_levels{$level}{"Unclassified"} += $n_unclassified;
		my $reads_all_taxa = 0;
		my $reads_all_taxa_ignoreUnclassified = 0;
		foreach my $taxonID (keys %{$reads_at_levels{$level}})
		{
			my $taxonID_for_print = $taxonID;
			my $name; 
			if($taxonID eq 'Unclassified')
			{
				$name = 'Unclassified';
				$taxonID_for_print = 0;
			}
			elsif($taxonID eq 'Undefined')
			{
				die;
				$name = 'NotLabelledAtLevel'; 
				$taxonID_for_print = -1;				
			}
			elsif($taxonID eq 'NotLabelledAtLevel')
			{
				die;
				$name = 'NotLabelledAtLevel'; 
				$taxonID_for_print = -1;				
			}			
			else
			{
				$name = taxTree::taxon_id_get_name($taxonID, $taxonomy_kraken);
			}
			
			my $nReads  = $reads_at_levels{$level}{$taxonID};
			
			print OUTPUT join("\t", $level, $taxonID_for_print, $name, $nReads, $nReads / $n_total_reads), "\n";
			print OUTPUT2 join("\t", $level, $taxonID_for_print, $name, ($taxonID eq 'Unclassified') ? ($nReads - $n_unclassified): $nReads, $nReads / $n_root), "\n";
			
			$reads_all_taxa += $nReads;
			$reads_all_taxa_ignoreUnclassified += (($taxonID eq 'Unclassified') ? ($nReads - $n_unclassified): $nReads);
		}
		
		die unless($reads_all_taxa == $n_total_reads);
		die unless($reads_all_taxa_ignoreUnclassified == $n_root);
	}
	close(OUTPUT);
	close(OUTPUT2);
}

sub create_compatible_composition_file_from_centrifuge
{
	my $target_output_fn = shift;
	my $taxonomy_kraken_dir = shift;
	my $f_distribution_krakenFormat = shift;
	my $f_reads = shift;
	my $create_compatible_file_from_kraken = shift;	
	die unless(defined $create_compatible_file_from_kraken);
	
	my $taxonomy_kraken = taxTree::readTaxonomy($taxonomy_kraken_dir);
	
	my $target_output_fn_2 = $target_output_fn . '.ignoreUnclassified';
	my %S_byLevel;
	my $n_unclassified;
	my $n_root;
	
	open(KRAKEN, '<', $f_distribution_krakenFormat) or die "Cannot open $f_distribution_krakenFormat";
	while(<KRAKEN>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		if($f[5] eq 'unclassified')
		{
			die if(defined $n_unclassified);
			$n_unclassified = $f[1];
			next;
		}
		elsif($f[5] eq 'root')
		{
			die if(defined $n_root);
			$n_root = $f[1];
			next;
		}
	}
	
	$n_unclassified = 0 unless(defined $n_unclassified);
	my $n_total_reads = $n_unclassified + $n_root;
	
	my $n_unclassified_check = 0;
	
	my %_getLightning_cache;
	my $getLightning = sub {
		my $taxonID = shift;
		if(exists $_getLightning_cache{$taxonID})
		{
			return $_getLightning_cache{$taxonID};
		}
		else
		{
			unless(($taxonID eq '0') or (exists $taxonomy_kraken->{$taxonID}))
			{
				warn "Not sure what to do with taxon ID $taxonID (file $f_reads)";
			}
			my $lightning = validation::getAllRanksForTaxon_withUnclassified($taxonomy_kraken, $taxonID, $create_compatible_file_from_kraken);
			$_getLightning_cache{$taxonID} = $lightning;
			return $lightning;
		}
	}; 
	
	if($n_total_reads == 0)
	{
		die "Weird - $n_total_reads reads? Output $target_output_fn, kraken output $f_distribution_krakenFormat";
	}
	
	my @evaluateAccuracyAtLevels = validation::getEvaluationLevels();
	
	my %reads_at_levels;
	open(CENTRIFUGE, '<', $f_reads) or die "Cannot open $f_reads";
	my $headerLine = <CENTRIFUGE>;
	my @header_fields = split(/\t/, $headerLine);
	die unless($header_fields[0] eq 'readID');
	die unless($header_fields[1] eq 'seqID');
	die unless($header_fields[2] eq 'taxID');
	while(<CENTRIFUGE>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		my $readID = $f[0];
		my $seqID = $f[1];
		my $taxID = $f[2];
		if(($seqID eq 'unclassified') or ($taxID eq '0'))
		{
			$n_unclassified_check++;		
		}
		else
		{
			die "Line $. of file $f_reads - taxID is $taxID, but not unclassified?\n$line" if($taxID eq '0');
			my $lightning = $getLightning->($taxID);
			$reads_at_levels{'definedAndHypotheticalGenomes'}{$taxID}++;
			RANK: foreach my $rank (@evaluateAccuracyAtLevels)
			{
				die unless(defined $lightning->{$rank});
				$reads_at_levels{$rank}{$lightning->{$rank}}++;
			}
		}
	}
	close(CENTRIFUGE);

	die "Inconsistency w.r.t. unclassified reads -- $n_unclassified_check vs $n_unclassified" unless($n_unclassified_check == $n_unclassified);
	
	open(OUTPUT, '>', $target_output_fn) or die "Cannot open $target_output_fn";
	open(OUTPUT2, '>', $target_output_fn_2) or die "Cannot open $target_output_fn_2";
	print OUTPUT join("\t", qw/AnalysisLevel ID Name Absolute PotFrequency/), "\n";		
	print OUTPUT2 join("\t", qw/AnalysisLevel ID Name Absolute PotFrequency/), "\n";		
	
	foreach my $level ('definedAndHypotheticalGenomes', @evaluateAccuracyAtLevels)
	{
		$reads_at_levels{$level}{"Unclassified"} = 0 if(not exists $reads_at_levels{$level}{"Unclassified"});
		# $reads_at_levels{$level}{"Undefined"} = 0 if(not exists $reads_at_levels{$level}{"Undefined"});
		
		$reads_at_levels{$level}{"Unclassified"} += $n_unclassified;
		my $reads_all_taxa = 0;
		my $reads_all_taxa_ignoreUnclassified = 0;
		foreach my $taxonID (keys %{$reads_at_levels{$level}})
		{
			my $taxonID_for_print = $taxonID;
			my $name; 
			if($taxonID eq 'Unclassified')
			{
				$name = 'Unclassified';
				$taxonID_for_print = 0;
			}
			elsif($taxonID eq 'Undefined')
			{
				die;
				$name = 'NotLabelledAtLevel'; 
				$taxonID_for_print = -1;				
			}
			elsif($taxonID eq 'NotLabelledAtLevel')
			{
				die;
				$name = 'NotLabelledAtLevel'; 
				$taxonID_for_print = -1;				
			}			
			else
			{
				$name = taxTree::taxon_id_get_name($taxonID, $taxonomy_kraken);
			}
			
			my $nReads  = $reads_at_levels{$level}{$taxonID};
			
			print OUTPUT join("\t", $level, $taxonID_for_print, $name, $nReads, $nReads / $n_total_reads), "\n";
			print OUTPUT2 join("\t", $level, $taxonID_for_print, $name, ($taxonID eq 'Unclassified') ? ($nReads - $n_unclassified): $nReads, $nReads / $n_root), "\n";
			
			$reads_all_taxa += $nReads;
			$reads_all_taxa_ignoreUnclassified += (($taxonID eq 'Unclassified') ? ($nReads - $n_unclassified): $nReads);
		}
		
		die unless($reads_all_taxa == $n_total_reads);
		die unless($reads_all_taxa_ignoreUnclassified == $n_root);
	}
	close(OUTPUT);
	close(OUTPUT2);
}

sub create_compatible_reads_file_from_kraken
{
	my $output_fn = shift;
	my $taxonomy_kraken_dir = shift;
	my $f_reads = shift;
	
	my $taxonomy_kraken = taxTree::readTaxonomy($taxonomy_kraken_dir);
	
	my $output_fn_unclassified = $output_fn . '.unclassified';
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";	
	open(OUTPUT_UNCL, '>', $output_fn_unclassified) or die "Cannot open $output_fn_unclassified";	
	open(KRAKEN, '<', $f_reads) or die "Cannot open $f_reads";
	while(<KRAKEN>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		my $classified = $f[0];
		my $readID = $f[1];
		my $taxonID = $f[2];
		die unless(($classified eq 'C') or ($classified eq 'U'));
		if($classified eq 'C')
		{
			print OUTPUT $readID, "\t", $taxonID, "\n";
		}
		else
		{
			print OUTPUT $readID, "\t", 0, "\n";		 
			print OUTPUT_UNCL $readID, "\t", 'Unclassified', "\n";		
		}
	}
	close(KRAKEN);
	close(OUTPUT);
	close(OUTPUT_UNCL);
	
}

sub create_compatible_reads_file_from_centrifuge
{
	my $output_fn = shift;
	my $taxonomy_kraken_dir = shift;
	my $f_reads = shift;
	
	my $taxonomy_kraken = taxTree::readTaxonomy($taxonomy_kraken_dir);
	
	my $output_fn_unclassified = $output_fn . '.unclassified';
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";	
	open(OUTPUT_UNCL, '>', $output_fn_unclassified) or die "Cannot open $output_fn_unclassified";	
	open(CENTRIFUGE, '<', $f_reads) or die "Cannot open $f_reads";
	my $headerLine = <CENTRIFUGE>;
	my @header_fields = split(/\t/, $headerLine);
	die unless($header_fields[0] eq 'readID');
	die unless($header_fields[1] eq 'seqID');
	die unless($header_fields[2] eq 'taxID');
	while(<CENTRIFUGE>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		my $readID = $f[0];
		my $seqID = $f[1];
		my $taxID = $f[2];
		if($seqID eq 'unclassified')
		{
			print OUTPUT $readID, "\t", 0, "\n";		 
			print OUTPUT_UNCL $readID, "\t", 'Unclassified', "\n";			
		}
		else
		{
			print OUTPUT $readID, "\t", $taxID, "\n";
		}
	}
	close(CENTRIFUGE);
	close(OUTPUT);
	close(OUTPUT_UNCL);
	
}

sub create_compatible_file_from_kraken_bracken
{
	my $output_fn = shift;
	my $taxonomy_kraken_dir = shift;
	
	my $f_K = shift;
	my $f_S = shift;
	my $f_G = shift;
	my $f_F = shift;	
	
	my $taxonomy_kraken = taxTree::readTaxonomy($taxonomy_kraken_dir);

	my $output_fn_2 = $output_fn . '.ignoreUnclassified';

	my $n_unclassified;
	my $n_root;
	open(KRAKEN, '<', $f_K) or die "Cannot open $f_K";
	while(<KRAKEN>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @f = split(/\t/, $line);
		if($f[5] eq 'unclassified')
		{
			die if(defined $n_unclassified);
			$n_unclassified = $f[1];
			next;
		}
		elsif($f[5] eq 'root')
		{
			die if(defined $n_root);
			$n_root = $f[1];
			next;
		}
	}
	
	my $n_total_reads = $n_unclassified + $n_root;
	
	print "Reads unclassified by Kraken $n_unclassified / classified $n_root\n";

	my $read_S = sub {
		my $fn = shift;
		my $rank = shift;
		my $ignoreUnclassfied = shift;

		my $n_reads_classified_level = 0;
		my %S;
		open(S, '<', $fn) or die "Cannot open $fn";
		my $header_line = <S>;
		chomp($header_line);
		my @header_fields = split(/\t/, $header_line);
		die unless($header_fields[1] eq 'taxonomy_id');
		die unless($header_fields[5] eq 'new_est_reads');
		die unless($header_fields[6] eq 'fraction_total_reads');
		while(<S>)
		{
			my $line = $_;
			chomp($line);
			next unless($line);
			my @f = split(/\t/, $line);
			my $taxonID = $f[1];
			my $nReads = $f[5];
			my $fraction = $f[6];
			die "Unknown taxonomy ID $taxonID" unless(exists $taxonomy_kraken->{$taxonID});
			die "Weird rank for taxonomy ID $taxonID" unless($taxonomy_kraken->{$taxonID}{rank} eq $rank);
			
			if($ignoreUnclassfied)
			{
				$S{$taxonID}[0] += $nReads;
				$S{$taxonID}[1] += ($nReads / $n_root);
			}
			else
			{
				$S{$taxonID}[0] += $nReads;
				$S{$taxonID}[1] += ($nReads / $n_total_reads);			
			}
			$n_reads_classified_level += $nReads;
		}
		close(S);
		
		my $n_reads_unclassified_level = ($ignoreUnclassfied) ? ($n_root - $n_reads_classified_level) : ($n_total_reads - $n_reads_classified_level);
		# die unless($n_reads_classified_level >= $n_unclassified);

		$S{'Unclassified'}[0] = $n_reads_unclassified_level;
		$S{'Unclassified'}[1] = ($ignoreUnclassfied) ? ($n_reads_unclassified_level / $n_root) : ($n_reads_unclassified_level / $n_total_reads);
		
		return \%S;
	};
	
	my $results_species = $read_S->($f_S, 'species');
	my $results_genus = $read_S->($f_G, 'genus');
	my $results_family = $read_S->($f_F, 'family');
	
	my $results_species_ignoreUnclassified = $read_S->($f_S, 'species', 1);
	my $results_genus_ignoreUnclassified = $read_S->($f_G, 'genus', 1);
	my $results_family_ignoreUnclassified = $read_S->($f_F, 'family', 1);
	
	open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
	open(OUTPUT2, '>', $output_fn_2) or die "Cannot open $output_fn_2";
	
	print OUTPUT join("\t", qw/AnalysisLevel ID Name Absolute PotFrequency/), "\n";
	print OUTPUT2 join("\t", qw/AnalysisLevel ID Name Absolute PotFrequency/), "\n";
	
	my $print_S = sub {
		my $S_href = shift;
		my $level = shift;
		my $printToOutput2 = shift;
		
		foreach my $taxonID (keys %$S_href)
		{
			my $taxonID_for_print = $taxonID;
			my $name;
			if($taxonID eq 'Unclassified')
			{
				$name = 'Unclassified';
				$taxonID_for_print = 0;
			}
			elsif($taxonID eq 'Undefined')
			{
				die;
				$name = 'NotLabelledAtLevel'; 
				$taxonID_for_print = -1;				
			}
			elsif($taxonID eq 'NotLabelledAtLevel')
			{
				die;
				$name = 'NotLabelledAtLevel'; 
				$taxonID_for_print = -1;				
			}			
			else
			{
				$name = taxTree::taxon_id_get_name($taxonID, $taxonomy_kraken);
			}
				
			if($printToOutput2)
			{
				print OUTPUT2 join("\t", $level, $taxonID_for_print, $name, $S_href->{$taxonID}[0], $S_href->{$taxonID}[1]), "\n";
			}
			else
			{
				print OUTPUT join("\t", $level, $taxonID_for_print, $name, $S_href->{$taxonID}[0], $S_href->{$taxonID}[1]), "\n";
			}
		}
	};
	
	$print_S->($results_species, 'species');
	$print_S->($results_genus, 'genus');
	$print_S->($results_family, 'family');

	$print_S->($results_species_ignoreUnclassified, 'species', 1);
	$print_S->($results_genus_ignoreUnclassified, 'genus', 1);
	$print_S->($results_family_ignoreUnclassified, 'family', 1);
	
	close(OUTPUT);
	close(OUTPUT2);
	
	print "\nCreated file $output_fn \n\n";
}


1;

