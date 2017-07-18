package Util;

use strict;

sub getGenomeLength
{
	my $taxonID = shift;
	my $taxon_2_contig = shift;
	my $contig_2_length = shift;
	
	die unless(defined $contig_2_length);
	
	my $gL = 0;
	die "Cannot determine genome length for taxon ID $taxonID" unless(defined $taxon_2_contig->{$taxonID});
	
	my @contigIDs;
	if(ref($taxon_2_contig->{$taxonID}) eq 'ARRAY')
	{
		@contigIDs = @{$taxon_2_contig->{$taxonID}}
	}
	elsif(ref($taxon_2_contig->{$taxonID}) eq 'HASH')
	{
		@contigIDs = keys %{$taxon_2_contig->{$taxonID}}
	}
	else
	{
		die;
	}
	foreach my $contigID (@contigIDs)
	{
		die unless(defined $contig_2_length->{$contigID});
		$gL += $contig_2_length->{$contigID};
	}
	
	return $gL;
}	

sub mean
{
	my $s = 0;
	die unless(scalar(@_));
	foreach my $v (@_)
	{
		$s += $v;
	}
	return ($s / scalar(@_));
}

sub sd
{
	die unless(scalar(@_));
	my $m = mean(@_);
	my $sd_sum = 0;
	foreach my $e (@_)
	{
		$sd_sum += ($m - $e)**2;
	}
	my $sd = sqrt($sd_sum);
	return $sd;
}

sub readFASTA
{
	my $file = shift;	
	my $cut_sequence_ID_after_whitespace = shift;
	
	my %R;
	
	open(F, '<', $file) or die "Cannot open $file";
	my $currentSequence;
	while(<F>)
	{		
		my $line = $_;
		chomp($line);
		$line =~ s/[\n\r]//g;
		if(substr($line, 0, 1) eq '>')
		{
			if($cut_sequence_ID_after_whitespace)
			{
				$line =~ s/\s+.+//;
			}
			$currentSequence = substr($line, 1);
			$R{$currentSequence} = '';
		}
		else
		{
			die "Weird input in $file" unless (defined $currentSequence);
			$R{$currentSequence} .= uc($line);
		}
	}	
	close(F);
		
	return \%R;
}


1;
