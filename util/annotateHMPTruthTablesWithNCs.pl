use strict;
use Data::Dumper;
use List::MoreUtils qw/mesh/;
use FindBin;
use lib "$FindBin::Bin/../perlLib";
use Util;
use taxTree;
$| = 1;

my $DB = 'databases/miniSeq+H';
my @files = ('/data/projects/phillippy/projects/MetaMap/tmp/truthHMP7_bwa_nanopore.distribution_genomes', '/data/projects/phillippy/projects/MetaMap/tmp/truthHMP7_bwa_pacbio.distribution_genomes');
my $mappingBaseFile = '/data/projects/phillippy/projects/mash_map/Jobs/blasr/hmp/targetAll/mock.all.genome.fa';

my $masterTaxonomy_dir = '/data/projects/phillippy/projects/MetaMap/downloads/taxonomy';
my $MetaMap_taxonomy = taxTree::readTaxonomy($masterTaxonomy_dir);
my $MetaMap_taxonomy_merged = taxTree::readMerged($masterTaxonomy_dir);

# also see /data/projects/phillippy/projects/MetaMap/util/truthForHMP.pl, where this file being generated
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


my $taxonomyDir = $DB . '/taxonomy';
my $DBtaxonomy = taxTree::readTaxonomy($taxonomyDir);


my %gi_2_contigID;
my %taxonIDs_in_mapping_ref;
my %taxonIDs_old_in_mapping_ref;
open(MAPPINGBASE, '<', $mappingBaseFile) or die "Cannot open $mappingBaseFile";
while(<MAPPINGBASE>)
{
	if(substr($_, 0, 1) eq '>')
	{
		chomp;
		substr($_, 0, 1) = '';
		my $contigID = $_;
		die "Invalid contig ID - no GI! $contigID " unless($contigID =~ /gi\|(\d+)\|/);		
		my $gi = $1;
		push(@{$gi_2_contigID{$1}}, $contigID);
		die Dumper("No GI-2-Taxon for GI $gi", [keys %gi_2_taxon]) unless(exists $gi_2_taxon{$gi});		
		my $taxonID_original = $gi_2_taxon{$gi};
		$taxonIDs_old_in_mapping_ref{$taxonID_original}++;
		
		my $taxonID_current = taxTree::findCurrentNodeID($MetaMap_taxonomy, $MetaMap_taxonomy_merged, $taxonID_original);
		$taxonIDs_in_mapping_ref{$taxonID_current}++;
	}
}
close(MAPPINGBASE);

my %taxon_2_gi;
foreach my $gi (keys %gi_2_taxon)
{
	my $taxonID = $gi_2_taxon{$gi};
	my $taxonID_current = taxTree::findCurrentNodeID($MetaMap_taxonomy, $MetaMap_taxonomy_merged, $taxonID);	
	push(@{$taxon_2_gi{$taxonID_current}}, $gi);
}

my %taxonID_2_contigs;
my %contigLength;
Util::read_taxonIDs_and_contigs($DB, \%taxonID_2_contigs, \%contigLength);

foreach my $file (@files)
{
	my %printed_taxonID;
	my $f_out = $file . '.withGI';
	open(F, '<', $file) or die "Cannot open $file";
	open(FOUT, '>', $f_out) or die "Cannot open $f_out for writing";
	my $headerLine = <F>;
	chomp($headerLine);
	my @headerFields = split(/\t/, $headerLine);
	print FOUT join("\t", @headerFields, 'GIs'), "\n";
	while(<F>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @lineFields = split(/\t/, $line);
		die unless($#lineFields == $#headerFields);
		my %line = (mesh @headerFields, @lineFields);
		my $taxonID = $line{taxonID};
		
		die unless($MetaMap_taxonomy->{$taxonID});
		die Dumper("Unknown taxon ID $taxonID", $taxonIDs_old_in_mapping_ref{$taxonID}) unless($taxonIDs_in_mapping_ref{$taxonID});
		die "Unknown taxon ID $taxonID (II)" unless($taxon_2_gi{$taxonID});
				
		my @gis = @{$taxon_2_gi{$taxonID}};
		
		print FOUT join("\t", @lineFields, join(',', @gis)), "\n";
		
		# my @contigs = keys %{$taxonID_2_contigs{$taxonID}};
		# my %NCs = map {$_ => 1} map {die unless($_ =~ /\|(NC_.+)/); $1} @contigs;
		# my @NCs = (keys %NCs);
		# print FOUT join("\t", @lineFields, join(',', @NCs)), "\n";
		
		$printed_taxonID{$taxonID} = 1;
	}
	close(F);
	
	foreach my $mappingTaxonID (keys %taxonIDs_in_mapping_ref)
	{	
		next if($printed_taxonID{$mappingTaxonID});
		my $name = taxTree::taxon_id_get_name($mappingTaxonID, $MetaMap_taxonomy);
		my @gis = @{$taxon_2_gi{$mappingTaxonID}};		
		my @line_fields = ($mappingTaxonID, $name, 0, 0, 0, 0, join(',', @gis));
		print FOUT join("\t", @line_fields), "\n";

	}
	close(FOUT);
	print "Produced $f_out\n";
}	