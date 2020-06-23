#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_n $opt_q $opt_l $opt_h $opt_L $opt_H $opt_c );
use strict;

&init();
my ($pos, $ref, $alt, $id, $f, $gq, $depth, $A0, $A1, 
    $nREF, $nALT, $pIBD0, $pIBD1, $pIBD2, $l, $bin);
my $bin_start = 1;
my $aggr_ibd0 = 1;
my $aggr_ibd1 = 1;
my $aggr_ibd2 = 1;
my ( $bin_end );

# Initialize bins
my $num_sites = 0;

open( COM, $opt_c ) or die( "$opt_c: $!\n" );
while( chomp( $l =<COM> ) ) {
    if ( $l =~ /^#/ ) { # skip header / comment lines
	next;
    }
    # Parse an input line
    ($pos, $ref, $alt, $id, $f, $gq, $depth, $A0, $A1, 
     $nREF, $nALT, $pIBD0, $pIBD1, $pIBD2) 
	= split( "\t", $l );

    # Check if it passes user-defined filters
    if ( &passes_filters( $gq, $depth, ($nREF + $nALT) ) ) {

	# If we've seen a bin's worth of sites
	if ( $num_sites == $opt_n ) {
	    # output the data in the bin
	    &output_bin( $bin_start, $bin_end, 
			 $aggr_ibd0, $aggr_ibd1, $aggr_ibd2,
			 $num_sites );
	    # ...then reset everything for the next bin
	    $num_sites = 0;
	    $bin_start = $pos;
	    $aggr_ibd0 = 1;
	    $aggr_ibd1 = 1;
	    $aggr_ibd2 = 1;
	}
	# Update aggregate likelihoods with data from this site
	$aggr_ibd0 *= $pIBD0;
	$aggr_ibd1 *= $pIBD1;
	$aggr_ibd2 *= $pIBD2;
	$bin_end = $pos;
	$num_sites++;
    }
}
close( COM );

if ( $num_sites > 0 ) {
    &output_bin( $bin_start, $bin_end, 
		 $aggr_ibd0, $aggr_ibd1, $aggr_ibd2,
		 $num_sites );
}

sub output_bin {
    my $bin_start = shift;
    my $bin_end   = shift;
    my $aggr_ibd0 = shift;
    my $aggr_ibd1 = shift;
    my $aggr_ibd2 = shift;
    my $num_sites = shift;
    printf( "%d %d %8.7e %8.7e %8.7e %d\n", 
	    $bin_start,
	    $bin_end,
	    $aggr_ibd0, 
	    $aggr_ibd1, 
	    $aggr_ibd2,
	    $num_sites );
}

sub passes_filters {
    my $gq    = shift;
    my $depth = shift;
    my $nOBS  = shift;

    if ( $gq < $opt_q ) {
	return 0;
    }
    if ( ($depth < $opt_l) || ($depth > $opt_h) ) {
	return 0;
    }
    if ( ($nOBS < $opt_L) || ($nOBS > $opt_H) ) {
	return 0;
    }
    return 1;
}
    
sub init {
    my $q_DEF = 40;
    my $l_DEF = 5;
    my $h_DEF = 12;
    my $L_DEF = 1;
    my $H_DEF = 4;
    my $n_DEF = 100;
    getopt( 'n:q:l:h:L:H:c:' );
    unless( -f $opt_c ) {
	print( "summarize-comparison.pl Summarizes output of compare-vcf2bam.pl\n" );
	print( "Writes the aggregate probabilities across bins of the input file.\n" );
	print( "These bins can be plotted to see regional trends.\n" );
	print( "Input data must be sorted by position.\n" );
	print( "Output format is: POS_START POS_END AGGR_IBD0 AGGR_IBD1 AGGR_IBD2 NUM_SITES\n" );
	print( "-n <number of sites per bin; default = $n_DEF>\n" );
	print( "-q <Genotype Quality cutoff; default = $q_DEF>\n" );
	print( "-l <Low coverage cutoff for reference genotype; default = $l_DEF>\n" );
	print( "-h <High coverage cutoff for reference genotype; default = $h_DEF>\n" );
	print( "-L <Low coverage cutoff for comparison data; default = $L_DEF>\n" );
	print( "-H <High coverage cutoff for comparison data; default = $H_DEF>\n" );
	print( "-c <input comparison table file>\n" );
	exit( 0 );
    }
    unless( defined($opt_n) ) {
	$opt_n = $n_DEF;
    }
    unless( defined($opt_q) ) {
	$opt_q = $q_DEF;
    }
    unless( defined($opt_l) ) {
	$opt_l = $l_DEF;
    }
    unless( defined($opt_h) ) {
	$opt_h = $h_DEF;
	}
    unless( defined($opt_L) ) {
	    $opt_L = $L_DEF;
    }
    unless( defined($opt_H) ) {
	$opt_H = $H_DEF;
    }
}

