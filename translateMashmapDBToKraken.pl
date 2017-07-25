#!/usr/bin/perl

# test command: perl translateMashmapDBToKraken.pl --input /data/projects/phillippy/projects/mashsim/src/simulations/0/DB.fa
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin/perlLib";

use taxTree;
use Util;

my $taxonomyDir = '/data/projects/phillippy/projects/mashsim/NCBI/refseq/taxonomy/';
my $krakenTemplate_taxonomy = '/data/projects/phillippy/projects/mashsim/src/krakenDBTemplate/taxonomy/';
my $krakenTemplate_taxonomy_names = '/data/projects/phillippy/projects/mashsim/src/krakenDBTemplate/taxonomy/names.dmp';

$taxonomyDir = '';
$krakenTemplate_taxonomy = '';

my $input_fn;
my $really;
GetOptions (
	'input:s' => \$input_fn,
	'taxonomyDir:s' => \$taxonomyDir,
	'krakenTemplate_taxonomy:s' => \$krakenTemplate_taxonomy,
);
my $output_fn = $input_fn . '.kraken';

# read taxonomy

my $taxonomy = taxTree::readTaxonomy($taxonomyDir);

# read IDs in kraken template DB

my %kraken_known_ID;
open(KRAKENIDs, '<', $krakenTemplate_taxonomy . '/names.dmp') or die "Cannot open names in $krakenTemplate_taxonomy";
while(<KRAKENIDs>)
{	
	chomp;
	next unless($_);
	die unless($_ =~ /^(\d+)/);
	$kraken_known_ID{$1}++;
}
close(KRAKENIDs);
my %merged;
open(MERGED, '<', $krakenTemplate_taxonomy . '/merged.dmp') or die "Cannot open merged in $krakenTemplate_taxonomy";
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

open(INPUT, '<', $input_fn) or die "Cannot open $input_fn";
open(OUTPUT, '>', $output_fn) or die "Cannot open $output_fn";
while(<INPUT>)
{
	if(substr($_, 0, 1) eq '>')
	{
		my $taxonID = Util::extractTaxonID($_, $input_fn, $.);
		
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
		
		die "Taxon ID $newID is unknown to kraken - $taxonID - $newID" unless($kraken_known_ID{$newID});
		
		$_ =~ s/kraken:taxid\|([x\d]+?)\|/kraken:taxid\|$newID\|/;	
	}
	
	print OUTPUT $_;
}
close(OUTPUT);
close(INPUT);

print "\nProduced file $output_fn\n\n";

