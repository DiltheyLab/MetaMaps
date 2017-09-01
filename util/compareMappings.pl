#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

unless(scalar(@ARGV) == 2)
{
	die "USAGE: perl compareMappings.pl FILE1 FILE2";
}

my $m_f1 = get_mappings_href($ARGV[0]);
my $m_f2 = get_mappings_href($ARGV[1]);

my @s = setstats($m_f1, $m_f2);

print "Statistics:\n\n\tS1: $ARGV[0]\n\tS2: $ARGV[1]\n\n";
print "s1-exclusive: $s[0]", "\n";
print "UNION       : $s[1]", "\n";
print "s2-exclusive: $s[2]", "\n";

print "\n";

sub setstats
{
	my $h1 = shift;
	my $h2 = shift;
	
	my $a = 0;
	my $b = 0;
	my $c = 0;
	
	my %_u = (map {$_ => 1} ((keys %$h1), (keys %$h2)));
	my @u = keys %_u;
	
	for(@u)
	{
		if((exists $h1->{$_}) and (exists $h2->{$_}))
		{
			$b++;
		}
		elsif((exists $h1->{$_}) and (not exists $h2->{$_}))
		{
			$a++;
		}
		elsif((not exists $h1->{$_}) and (exists $h2->{$_}))
		{
			$c++;
		}
		else
		{
			die Dumper($_, (exists $h1->{$_}), (exists $h2->{$_}));
			die;
		}
	}
	return ($a, $b, $c);
}	

sub get_mappings_href
{	
	my $fn = shift;
	my $forReturn_href = {};
	open(F, '<', $fn) or die "Cannot open $fn";
	while(<F>)
	{
		my $l = $_;
		chomp($l);
		my @f = split(/ /, $l);
		my $k = join('_', @f[0, 2, 3, 4, 5, 7, 8]);
		$forReturn_href->{$k}++;
	}	
	close(F);
	return $forReturn_href;
}