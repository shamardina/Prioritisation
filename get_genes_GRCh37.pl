use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Transcript;


my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -port => 3337
    );

my $file = $ARGV[0] or die "Provide gene-transcriptID tab-separated file\n";
open(my $data, '<', $file) or die "Could not open '$file' $!\n";
my %genes;
while (<$data>) {
    chomp();
    my @fields = split("\t");
    # # Transcript version is now removed when processing gene list, hence commented out:
    # $fields[1] =~ s/\.\d+$//;
    $genes{$fields[0]} = $fields[1];
}

my $transcript_adaptor = $registry->get_adaptor('Human', 'Core', 'Transcript');

for my $gene (keys(%genes)) {
    my $transcript = $transcript_adaptor->fetch_by_stable_id($genes{$gene});
    foreach my $exon (@{$transcript->get_all_Exons()}) {
	print(join("\t", $exon->seq_region_name(), $exon->start() - 1, $exon->end(), $transcript->stable_id(), $gene), "\n");
    }
}

