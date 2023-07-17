use strict;
use warnings;
use Data::Dumper;

unless(scalar(@ARGV) == 1)
{
	print_help();
}

my $DB = $ARGV[0];

my %contigTranslation;

my $taxons_in = $DB . '/taxonInfo.txt';
my $taxons_out = $DB . '/taxonInfo.txt.2';
open(TAXONS, '<', $taxons_in) or die;
open(TAXONOUT, '>', $taxons_out) or die;
while(<TAXONS>)
{
	my $line = $_;
	chomp($line);
	next unless($line);
	my @fields = split(/ /, $line);
	my $taxonID = shift(@fields);
	my $remainingLine = join(' ', @fields);
	my @contigs_pre = split(/(\=\d+;)/, $remainingLine);
	my @contigs;
	for(my $i = 0; $i <= $#contigs_pre; $i++)
	{
		if(($i % 2) == 0)
		{
			push(@contigs, $contigs_pre[$i]);
		}	
		else
		{
			die unless($contigs_pre[$i] =~ /;$/);
			$contigs_pre[$i] =~ s/;$//;			
			$contigs[$#contigs] .= $contigs_pre[$i];
		}
	}
	my @newContigFields;
	foreach my $contigID (@contigs)
	{
		die unless($contigID =~ /^(.+)=(\d+)$/);
		my $realContigID = $1;
		my $contigLength = $2;
		die unless($realContigID =~ /^(\S+)\s/);
		my $firstThingBeforeSpace = $1;	
		die if(exists $contigTranslation{$realContigID});
		$contigTranslation{$realContigID} = $firstThingBeforeSpace;
		push(@newContigFields, $firstThingBeforeSpace . '=' . $contigLength);
	}
	
	print TAXONOUT $taxonID, " ", join(';', @newContigFields), "\n";
}
close(TAXONS);
close(TAXONOUT);

my $DB_in = $DB . '/DB.fa';
my $DB_out = $DB . '/DB.fa.2';
open(DB, '<', $DB_in) or die;
open(DBOUT, '>', $DB_out) or die;
while(<DB>)
{
	my $line = $_;
	next unless($line);
	if(substr($line, 0, 1) eq '>')
	{	
		chomp($line);
		my $contigID = substr($line, 1);	
		die "Undefined contig ID $contigID" unless(exists $contigTranslation{$contigID});
		print DBOUT '>', $contigTranslation{$contigID}, "\n";
	}
	else
	{
		print DBOUT $line;
	}
}
close(DB);
close(DBOUT);

print "\n\nProduced files\n\t$taxons_out\n\t$DB_out\n\nBackup the original files, replace the original files with the .2 files, and validate the database.\n\n";

sub print_help
{
	print qq(
shortenContigIDs.pl

  Strip down contig IDs in a database.
  
Usage:

  perl shortenContigIDs.pl dbNAME
  
Example:

  perl shortenContigIDs.pl databases/refseq
  
  
);
exit;
}
