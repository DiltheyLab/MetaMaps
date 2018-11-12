#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my $inputFA;
my $outputFA;
my $taxonID;
GetOptions (
	'inputFA:s' => \$inputFA, 
	'outputFA:s' => \$outputFA, 
	'taxonID:s' => \$taxonID, 
);

unless($inputFA and $outputFA and $taxonID)
{
	print_help();
}

unless(-e $inputFA)
{
	die "File $inputFA (from parameter --inputFA) does not exist";
}

my $contigI = 0;
open(FIN, '<', $inputFA) or die "Cannot open $inputFA";
open(FOUT, '>', $outputFA) or die "Cannot open $outputFA for writing";
while(<FIN>)
{
	if(substr($_, 0, 1) eq '>')
	{
		$contigI++; 
		$_ =~ s/\s.+//;
		substr($_, 0, 1) = '>kraken:taxid|' . $taxonID . '|C' . $contigI . '|';
	}
	print FOUT $_;
}
close(FIN);
close(FOUT);

sub print_help
{
	print qq(
addTaxonIDToFasta.pl

  Add taxon ID in MetaMap format to sequence identifiers in FASTA file.
  
Usage:

  perl addTaxonIDToFasta.pl --inputFA inputFasta.fa --outputFA outputFasta.fa --taxonID 9604
  
Parameters:

  inputFA
      
	  Input fasta file.
  
  outputFA
  
      Output fasta file.
      
  taxonID
      
      Taxon ID to be added, e.g. 9604 for human.
	  
	  
);
exit;
}

