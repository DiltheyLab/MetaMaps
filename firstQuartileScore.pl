use strict;
use warnings;
use Data::Dumper;

if(scalar(@ARGV) == 0)
{
	die "Please provide FASTQ as first argument";
}

my $FASTQ = $ARGV[0];

my %Q;
my $totalQ = 0;
open(F, '<', $FASTQ) or die "Cannot open $FASTQ";
while(<F>)	
{
	my $l = $_;
	die unless(substr($l, 0, 1) eq '@');
	<F>;
	my $ll = <F>;
	die unless(substr($ll, 0, 1) eq '+');
	my $q = <F>;
	chomp($q);
	my @Q = split(//, $q);
	for(@Q)
	{
		$Q{$_}++;
		$totalQ++;
	}	
	
	# last if($. > 100);
}
close(F);

my $quartileQ;
my $quantile = 0;
foreach my $q (sort {$a cmp $b} keys %Q)
{
	$quantile += ($Q{$q}/$totalQ);
	print STDERR join("\t", ord($q),  $q, $Q{$q}, $quantile), "\n";
	if((not defined $quartileQ) and ($quantile >= 0.25))
	{
		$quartileQ = $q;
	}
}
unless(defined $quartileQ)
{
	die "Could not determine quartile quality";
}

print $quartileQ, "\n";