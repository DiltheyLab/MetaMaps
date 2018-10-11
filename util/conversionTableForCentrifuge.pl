#!/usr/bin/perl

# test command: perl translateMashmapDBToKraken.pl --input /data/projects/phillippy/projects/mashsim/src/simulations/0/DB.fa
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/../perlLib";

use taxTree;
use Util;
use SimulationsKraken;

my $suitableMasterTaxonomy = SimulationsKraken::getKrakenDBTemplate() . '/taxonomy';
die "$suitableMasterTaxonomy does not exist" unless(-d $suitableMasterTaxonomy);

# my $krakenTemplate_taxonomy_names = '/data/projects/phillippy/projects/mashsim/src/krakenDBTemplate/taxonomy/names.dmp';

my $DB;
GetOptions (
	'DB:s' => \$DB,
);

unless($DB)
{
	die "Please specify --DB";
}

my $DB_fasta = $DB . '/DB.fa';
die "File $DB_fasta does not exist" unless(-e $DB_fasta);

my $output_fn = $DB_fasta . '.centrifugeTranslation';
my $output_fn_names = $DB_fasta . '.centrifugeTranslation.names.dmp';
my $output_fn_nodes = $DB_fasta . '.centrifugeTranslation.nodes.dmp';

# read taxonomy

my $taxonomy = taxTree::readTaxonomy($DB . '/taxonomy');

my %centrifuge_known_ID;
open(KRAKENIDs, '<', $suitableMasterTaxonomy . '/names.dmp') or die "Cannot open names in $suitableMasterTaxonomy";
while(<KRAKENIDs>)
{	
	chomp;
	next unless($_);
	die unless($_ =~ /^(\d+)/);
	$centrifuge_known_ID{$1}++;
}
close(KRAKENIDs);
my %merged;
open(MERGED, '<', $suitableMasterTaxonomy . '/merged.dmp') or die "Cannot open merged in $suitableMasterTaxonomy";
while(<MERGED>)
{	
	chomp;
	next unless($_);
	die unless($_ =~ /^(\d+)\s+\|+\s+(\d+)/);
	my $from = $1;
	my $to = $2;
	die if(defined $merged{$from});
	$merged{$from} = $to;
}
close(MERGED);

# read input / substitute IDs

open(INPUT, '<', $DB_fasta) or die "Cannot open $DB_fasta";
open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
while(<INPUT>)
{
	if(substr($_, 0, 1) eq '>')
	{
		chomp;
		die "Invalid space in line $. of $DB_fasta" if($_ =~ /\s/);
		my $taxonID = Util::extractTaxonID($_, $DB_fasta, $.);
		
		die "Taxon ID not defined" unless(exists $taxonomy->{$taxonID});
		
		my $newID = $taxonID;
		if(substr($newID, 0, 1) eq 'x')
		{		
			my @ancestors = taxTree::get_ancestors($taxonomy, $taxonID);
			$newID = $ancestors[0];
			die unless($newID =~ /^\d+$/);
			print "Substitute $taxonID -> $newID \n";			
		}
		
		while(exists $merged{$newID})
		{
			$newID = $merged{$newID};
		}	
		
		die "Taxon ID $newID is unknown to centrifuge - $taxonID - $newID" unless($centrifuge_known_ID{$newID});
		
		die unless($_ =~ /^>(.+?\|.+?)\|/);
		my $id_for_centrifuge = $1;
		
		print OUTPUT join("\t", $id_for_centrifuge, $newID), "\n";
	}
}
close(OUTPUT);
close(INPUT);

system("cp ${suitableMasterTaxonomy}/names.dmp $output_fn_names") and die "Cannot copy names to $output_fn_names";
system("cp ${suitableMasterTaxonomy}/nodes.dmp $output_fn_nodes") and die "Cannot copy nodes to $output_fn_nodes";

print "\nProduced file $output_fn\n\n";
