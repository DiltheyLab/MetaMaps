use strict;
use Storable qw/dclone store retrieve/;
use Data::Dumper;

# my $simulationStore_dir = 'databases/miniSeq+H/simulations_i100_specifiedFrequencies_limitedMemory/0';
# my $simulationStore_dir = 'databases/miniSeq_100/simulations_logNormal/0';

my $simulation_href = retrieve $simulationStore_dir . '/simulationStore';
$simulation_href->{dbDirs_metamap} = [$simulationStore_dir . '/DB_fullDB'];
$simulation_href->{readsFastq} = $simulationStore_dir . '/reads.fastq';
$simulation_href->{outputDirectory} = $simulationStore_dir ;

foreach my $dir (@{$simulation_href->{dbDirs_metamap}})
{
	die unless(-e $dir);
}
die unless(-e $simulation_href->{readsFastq});
die unless(-e $simulation_href->{outputDirectory});

store $simulation_href, $simulationStore_dir . '/simulationStore';

print "\n\nOK - fixed store in $simulationStore_dir\n\n";
		  

