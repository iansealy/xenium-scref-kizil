#!/usr/bin/env perl

# PODNAME: dump_ensembl_genes.pl
# ABSTRACT: Get all Ensembl genes

use warnings;
use strict;
use autodie;
use Carp;
use Try::Tiny;

use Getopt::Long;
use Pod::Usage;
use Readonly;
use Bio::EnsEMBL::Registry;

=head1 DESCRIPTION


=head1 EXAMPLES

    perl \
        -Ibranch-ensembl-74/ensembl/modules \
        dump_ensembl_genes.pl \
        | sort > data/danio_rerio_e74.txt

=cut

# Constants
Readonly our $ENSEMBL_SPECIES => 'danio_rerio';

# Default options
my $ensembl_dbhost = 'ensembldb.ensembl.org';
my $ensembl_dbport;
my $ensembl_dbuser = 'anonymous';
my $ensembl_dbpass;
my $slice_regexp;
my ( $debug, $help, $man );

# Get and check command line options
get_and_check_options();

# Connnect to Ensembl database
Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => $ensembl_dbhost,
    -port => $ensembl_dbport,
    -user => $ensembl_dbuser,
    -pass => $ensembl_dbpass,
);

# Get genebuild version
my $genebuild_version = 'e' . Bio::EnsEMBL::ApiVersion::software_version();
warn 'Genebuild version: ', $genebuild_version, "\n" if $debug;

# Get Ensembl adaptors
my $sa =
  Bio::EnsEMBL::Registry->get_adaptor( $ENSEMBL_SPECIES, 'core', 'Slice' );

# Ensure database connection isn't lost; Ensembl 64+ can do this more elegantly
## no critic (ProhibitMagicNumbers)
if ( Bio::EnsEMBL::ApiVersion::software_version() < 64 ) {
## use critic
    Bio::EnsEMBL::Registry->set_disconnect_when_inactive();
}
else {
    Bio::EnsEMBL::Registry->set_reconnect_when_lost();
}

# Get all slices
my $slices = $sa->fetch_all('toplevel');
warn scalar @{$slices}, " slices\n" if $debug;

foreach my $slice ( @{$slices} ) {
    next
      if defined $slice_regexp
      && $slice->seq_region_name !~ m/$slice_regexp/xms;
    warn 'Slice: ', $slice->name, "\n" if $debug;

    # Get all genes
    my $genes = $slice->get_all_Genes( undef, 'core' );
    warn scalar @{$genes}, " genes\n" if $debug;
    foreach my $gene ( @{$genes} ) {
        printf "%s\n", join "\t", $gene->stable_id,
          ( $gene->external_name || $gene->stable_id );
    }
}

# Get and check command line options
sub get_and_check_options {

    # Get options
    GetOptions(
        'ensembl_dbhost=s' => \$ensembl_dbhost,
        'ensembl_dbport=i' => \$ensembl_dbport,
        'ensembl_dbuser=s' => \$ensembl_dbuser,
        'ensembl_dbpass=s' => \$ensembl_dbpass,
        'slice_regexp=s'   => \$slice_regexp,
        'debug'            => \$debug,
        'help'             => \$help,
        'man'              => \$man,
    ) or pod2usage(2);

    # Documentation
    if ($help) {
        pod2usage(1);
    }
    elsif ($man) {
        pod2usage( -verbose => 2 );
    }

    return;
}

=head1 USAGE

    dump_ensembl_genes.pl
        [--ensembl_dbhost host]
        [--ensembl_dbport port]
        [--ensembl_dbuser username]
        [--ensembl_dbpass password]
        [--slice_regexp regexp]
        [--debug]
        [--help]
        [--man]

=head1 OPTIONS

=over 8

=item B<--ensembl_dbhost HOST>

Ensembl MySQL database host.

=item B<--ensembl_dbport PORT>

Ensembl MySQL database port.

=item B<--ensembl_dbuser USERNAME>

Ensembl MySQL database username.

=item B<--ensembl_dbpass PASSWORD>

Ensembl MySQL database password.

=item B<--slice_regexp REGEXP>

Regular expression for limiting slices.

=item B<--debug>

Print debugging information.

=item B<--help>

Print a brief help message and exit.

=item B<--man>

Print this script's manual page and exit.

=back

=cut
