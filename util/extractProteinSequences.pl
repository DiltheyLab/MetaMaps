use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Set::IntervalTree;
use List::MoreUtils qw/mesh all/;
use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request::Common;

$| = 1;

my $DB = '/data/projects/phillippy/projects/MetaMap/databases/refseq_with_annotations';
my $proteinDB = '/data/projects/phillippy/projects/MetaMap/downloads/NCBI_proteins';
GetOptions (
	'DB:s' => \$DB, 
	'proteinDB:s' => \$proteinDB, 
);

die unless($DB and $proteinDB);

my $annotations_file = $DB . '/DB_annotations.txt';
die "DB $DB does not contain gene-level annotations" unless(-e $annotations_file);

my $mappings_file = $proteinDB . '/idmapping_selected.tab';
die "Protein DB $proteinDB does not contain ID mappings file $mappings_file" unless(-e $mappings_file);

my @gz_files = glob($proteinDB . '/*.fsa_aa.gz');
unless(@gz_files)
{
	die "Couldn't find any protein sequence files *.fsa_aa.gz in protein DB directory $proteinDB";
}	

my $fn_out_ids = $annotations_file . '.proteinIDs';
my $require_total_proteins_N = 0;
# unless(-e $fn_out_ids)
{
	my %proteinIDs;
	open(ANNOTATIONS, '<', $annotations_file) or die "Cannot open $annotations_file";
	my $annotations_header_line = <ANNOTATIONS>;
	chomp($annotations_header_line);
	my @header_fields = split(/\t/, $annotations_header_line);
	die unless($header_fields[$#header_fields - 1] eq 'CDSProduct');
	while(<ANNOTATIONS>)
	{
		my $line = $_;
		chomp($line);
		next unless($line);
		my @line_fields = split(/\t/, $line, -1);
		die unless($#line_fields == $#header_fields);
		my %line = (mesh @header_fields, @line_fields);
		die unless(defined $line{CDSProduct});
		my $proteinID = $line{CDSProduct};
		next unless($proteinID);
		$proteinIDs{$proteinID}++;
		# print "Want: '$proteinID'\n";
		# last if(scalar(keys %proteinIDs) > 100);
	}
	close(ANNOTATIONS);

	print "\nFound ", scalar(keys %proteinIDs), " proteins of interest.\n";

	open(PROTEINIDS, '>', $fn_out_ids) or die "Cannot open $fn_out_ids";
	print PROTEINIDS join("\n", keys %proteinIDs), "\n";
	close(PROTEINIDS);

	print "\nProduced $fn_out_ids\n\n";
	
	$require_total_proteins_N = scalar(keys %proteinIDs);
}

my %have_protein_sequences;
my $protein_sequences_fa = $annotations_file . '.proteinSequences.fa';
if(-e $protein_sequences_fa)
{
	open(EXISTING, '<', $protein_sequences_fa) or die "Cannot open $protein_sequences_fa";
	while(<EXISTING>)
	{
		if(substr($_, 0, 1) eq '>')
		{
			my $protein_id = substr($_, 1);
			chomp($protein_id);
			$protein_id =~ s/\s.*//;
			die "Duplicate protein ID" if($have_protein_sequences{$protein_id});
			$have_protein_sequences{$protein_id}++;
		}
	}
	close(EXISTING);
}
print "\nAlready have ", scalar(keys %have_protein_sequences), " protein sequences.\n";

print "Reading IDs from file $fn_out_ids \n";
my @running_wanted_IDs;
my %want_protein_sequences;
# acquire_protein_sequences(['abc123'], {'abc123' => 1}, {}, 'bla');
open(PROTEINIDS, '<', $fn_out_ids) or die "Cannot open $fn_out_ids";
while(<PROTEINIDS>)
{
	chomp;
	next unless($_);
	my $proteinID = $_;
	next if($have_protein_sequences{$proteinID});
	unless($want_protein_sequences{$proteinID})
	{
		$want_protein_sequences{$proteinID}++;	
		push(@running_wanted_IDs, $proteinID);
		if(scalar(@running_wanted_IDs) >= 5000)
		{
			acquire_protein_sequences(\@running_wanted_IDs, \%want_protein_sequences, \%have_protein_sequences, $protein_sequences_fa);
			
			print "Have ", scalar(keys %have_protein_sequences), " / ", $require_total_proteins_N, " protein sequences.\n";
		}
	}
}	
close(PROTEINIDS);
acquire_protein_sequences(\@running_wanted_IDs, \%want_protein_sequences, \%have_protein_sequences, $protein_sequences_fa);

sub acquire_protein_sequences
{
	my $proteins_aref = shift;
	my $want_protein_sequences_href = shift;
	my $have_protein_sequences_href = shift;
	my $output_fn = shift;
	
	print "Now query ", scalar(@$proteins_aref), " protein sequences.\n";
	
	die unless(all {$want_protein_sequences_href->{$_}} @$proteins_aref);
	
	my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';	
	my $ua = new LWP::UserAgent;
	$ua->agent("elink/1.0 " . $ua->agent);

	my $url = $base . '/efetch.fcgi';
	
	my $request = HTTP::Request::Common::POST( $url, [ db => 'protein', rettype => 'fasta', id => join(',', @$proteins_aref), email => 'alexander.dilthey@med.uni-duesseldorf.de'] );
	$request->content_type('application/x-www-form-urlencoded');

	my $response = $ua->request($request); 
	my $output = $response->content;

	my @lines = split(/\n/, $output);
	foreach my $line (@lines)
	{
		if(substr($line, 0, 1) eq '>')
		{
			my $protein_id = substr($line, 1);
			chomp($protein_id);
			$protein_id =~ s/\s.*//;
			
			die "Protein ID $protein_id was not requested" unless($want_protein_sequences_href->{$protein_id});
			die "Duplicate protein ID" if($have_protein_sequences_href->{$protein_id});
			delete $want_protein_sequences_href->{$protein_id};
			$have_protein_sequences_href->{$protein_id}++;
		}
	}

	open(OUTPUT, '>>', $output_fn) or die "Cannot open $output_fn";
	print OUTPUT $output;
	close(OUTPUT);
	
	# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=WP_012689726.1&usehistory=y
	
	my $got_proteins = 0;
	foreach my $proteinID (@$proteins_aref)
	{
		if(exists $have_protein_sequences_href->{$proteinID})
		{
			$got_proteins++;
		}
	}
	
	unless(scalar(@$proteins_aref) == $got_proteins)
	{
		warn "Of " . scalar(@$proteins_aref) . " requested protein sequences, got " . $got_proteins . "\n";
	}
	@$proteins_aref = ();
}

__END__

my %mapping_refseq_to_pir;
my %mapping_refseq_to_embl_cds;
my %mapping_pir_to_refseq;
open(TRANSLATIONS, '<', $mappings_file) or die "Cannot open $mappings_file";
while(<TRANSLATIONS>)
{
	my $line = $_;
	chomp($line);
	my @line_fields = split(/\t/, $line, -1);
	my $refseq_id = $line_fields[3];
	my $pir_id = $line_fields[11];
	$pir_id =~ s/\s//g;
	my @pir_ids = split(/;/, $pir_id);
	if($refseq_id and $pir_id and (exists $proteinIDs{$refseq_id}))
	{
		$mapping_refseq_to_pir{$refseq_id} = \@pir_ids;	
		foreach my $onePirID (@pir_ids)
		{
			#die Dumper("Inconsistency", "Existing", $mapping_pir_to_refseq{$onePirID}, "New", $refseq_id) if((exists $mapping_pir_to_refseq{$onePirID}) and ($mapping_pir_to_refseq{$onePirID} ne $refseq_id));
			$mapping_pir_to_refseq{$onePirID}{$refseq_id}++;
		}
	}
	my $embl_cds_id = $line_fields[17];
	$embl_cds_id =~ s/\s//g;
	my @embl_cds_ids = split(/;/, $embl_cds_id);
	if($refseq_id and $embl_cds_id and (exists $proteinIDs{$refseq_id}))
	{
		$mapping_refseq_to_embl_cds{$refseq_id} = \@embl_cds_ids;	
		foreach my $oneCDSID (@embl_cds_ids)
		{
			#die Dumper("Inconsistency", "Existing", $mapping_pir_to_refseq{$onePirID}, "New", $refseq_id) if((exists $mapping_pir_to_refseq{$onePirID}) and ($mapping_pir_to_refseq{$onePirID} ne $refseq_id));
			$mapping_pir_to_refseq{$oneCDSID}{$refseq_id}++;
		}
	}	
} 
close(TRANSLATIONS);
			 
my $output_file = $annotations_file . '.proteinSequences';
open(OUTPUT, '>', $output_file) or die "Cannot open $output_file";
open(OUTPUT2, '>', '_some_random_proteinSequences') or die "Cannot open _some_random_proteinSequences";
open(OUTPUTALLIDS, '>', "_all_protein_ids") or die "Cannot open _all_protein_ids";
my %proteinIDs_found;
my %proteinIDs_ambiguous;
print "Now processing files in $proteinDB...\n";
my $fileI = 0;
foreach my $file (@gz_files)
{
	$fileI++;
	print "\r\tProcessing file $fileI of ", scalar(@gz_files), "     ";
	
	my $currentRunningSequence;
	my $currentRunningId = '';
	my $currentRunningFullId = '';
	my $processOneSequence = sub {
		if($currentRunningId)
		{
			print OUTPUTALLIDS $currentRunningFullId, "\n";
			#print "Have: '$currentRunningId'", "\n";	
			my @possible_refseq_keys = (exists $mapping_pir_to_refseq{$currentRunningId}) ? (keys %{$mapping_pir_to_refseq{$currentRunningId}}) : ();
			my @lookup_keys = ($currentRunningId, @possible_refseq_keys);
			my @possible_keys = grep {exists $proteinIDs{$_}} @lookup_keys;
			foreach my $foundKey (@possible_keys)
			{
				die "Illegal ID with whitespace: $foundKey" if($foundKey =~ /\s/);
				print OUTPUT '>', $foundKey, ' ', $currentRunningFullId, "\n";
				print OUTPUT $currentRunningSequence, "\n";
				$proteinIDs_found{$foundKey}++;			
			}

			if(rand(1) < 0.001)
			{
				print OUTPUT2 '>', $currentRunningFullId, "\n";
				print OUTPUT2 $currentRunningSequence, "\n";				
			}
		}	
	};
	open(PIPE, "zless $file |") or die;
	while(<PIPE>)
	{
		next unless($_);
		if(substr($_, 0, 1) eq '>')
		{
			$processOneSequence->();
			$currentRunningSequence = '';
			$currentRunningFullId = substr($_, 1);
			chomp($currentRunningFullId);
			($currentRunningId = $currentRunningFullId)=~ s/\s.*//;
		}
		else
		{
			chomp;
			$currentRunningSequence .= $_;
		}
	}
	close(PIPE);
	$processOneSequence->();
}
close(OUTPUT);
close(OUTPUT2);
close(OUTPUTALLIDS);


print "\n";

print "Found ", scalar(grep {$proteinIDs_found{$_}} keys %proteinIDs), " of ", scalar(keys %proteinIDs), " desired proteins.\n\n";
my @missing_protein_ids = grep {not $proteinIDs_found{$_}} keys %proteinIDs;
print "Produced $output_file\n\n";
if(@missing_protein_ids)
{
	my $maxIndex = ($#missing_protein_ids >= 9) ? 9 : $#missing_protein_ids;
	print "Selection of missing protein IDs:\n";
	print join("\n", @missing_protein_ids[0 .. $maxIndex]), "\n";
}
my @ambiguous_translations = grep {(exists $proteinIDs_found{$_}) and ($proteinIDs_found{$_} > 1)} keys %proteinIDs;

print "Proteins with more than one sequence: ", scalar(@ambiguous_translations), "\n";
if(@ambiguous_translations)
{
	my $maxIndex = ($#ambiguous_translations >= 9) ? 9 : $#ambiguous_translations;
	print "... selection of these:\n";
	print join("\n", @ambiguous_translations[0 .. $maxIndex]), "\n";
}
