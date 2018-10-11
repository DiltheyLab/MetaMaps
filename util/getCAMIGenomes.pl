use strict;
use Data::Dumper;
use LWP::Simple;
use LWP::UserAgent;
use HTTP::Request::Common;

my $genomesOfOrigin = '/data/projects/phillippy/projects/MetaMap/tmp/truthCAMI.genomesOfOrigin';
my @genomeIDs;
open(GENOMES, '<', $genomesOfOrigin) or die;
while(<GENOMES>)
{
	chomp;
	next unless($_);
	push(@genomeIDs, $_);
}
close(GENOMES);


print "Now query for ", scalar(@genomeIDs), " protein sequences.\n";


my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';	
my $ua = new LWP::UserAgent;
$ua->agent("elink/1.0 " . $ua->agent);

my $url = $base . '/efetch.fcgi';

my $request = HTTP::Request::Common::POST( $url, [ db => 'nucleotide', rettype => 'fasta', id => join(',', @genomeIDs), email => 'alexander.dilthey@med.uni-duesseldorf.de'] );
$request->content_type('application/x-www-form-urlencoded');

my $response = $ua->request($request); 
my $output = $response->content;


open(OUTPUT, '>', $genomesOfOrigin . '.fasta') or die "Cannot open $genomesOfOrigin" . '.fasta';
print OUTPUT $output;
close(OUTPUT);

print "\nProduced $genomesOfOrigin.fasta \n";
