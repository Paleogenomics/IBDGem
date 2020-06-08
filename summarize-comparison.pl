#!/usr/local/bin/perl

use Getopt::Std;
use vars qw( $opt_b $opt_q $opt_l $opt_h $opt_L $opt_H $opt_c );
use strict;

&init();
my ($pos, $ref, $alt, $id, $f, $gq, $depth, $A0, $A1, $nREF, $nALT, $pIBD0, $pIBD1, $pIBD2,
    $l, $bin);
my (@ibd0_bins, @ibd1_bins, @ibd2_bins, @nsites_bins);

# Initialize bins
my $num_bins = int(250000000/$opt_b);
for( $bin = 0; $bin <= $num_bins; $bin++ ) {
    $ibd0_bins[$bin] = 1;
    $ibd1_bins[$bin] = 1;
    $ibd2_bins[$bin] = 1;
    $nsites_bins[$bin] = 0;
}

open( COM, $opt_c ) or die( "$opt_c: $!\n" );
while( chomp( $l =<COM> ) ) {
    if ( $l =~ /^#/ ) {
	next;
    }
    ($pos, $ref, $alt, $id, $f, $gq, $depth, $A0, $A1, $nREF, $nALT, $pIBD0, $pIBD1, $pIBD2) 
	= split( "\t", $l );
    if ( &passes_filters( $gq, $depth, ($nREF + $nALT) ) ) {
	$bin = int( $pos/$opt_b );
	$ibd0_bins[$bin] *= $pIBD0;
	$ibd1_bins[$bin] *= $pIBD1;
	$ibd2_bins[$bin] *= $pIBD2;
	$nsites_bins[$bin]++;
    }
}
close( COM );

for( $bin = 0; $bin <= $#ibd0_bins; $bin++ ) {
    printf( "%8.7e %8.7e %8.7e %d\n", 
	    $ibd0_bins[$bin],
	    $ibd1_bins[$bin],
	    $ibd2_bins[$bin],
	    $nsites_bins[$bin] );
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
    my $b_DEF = 1000000;
    my $q_DEF = 40;
    my $l_DEF = 5;
    my $h_DEF = 12;
    my $L_DEF = 1;
    my $H_DEF = 4;
    getopt( 'b:q:l:h:L:H:c:' );
    unless( -f $opt_c ) {
	print( "summarize-comparison.pl Summarizes output of compare-vcf2bam.pl\n" );
	print( "Writes the aggregate probabilities across bins of the input file.\n" );
	print( "These bins can be plotted to see regional trends.\n" );
	print( "-b <bin size; default = $b_DEF>\n" );
	print( "-q <Genotype Quality cutoff; default = $q_DEF>\n" );
	print( "-l <Low coverage cutoff for reference genotype; default = $l_DEF>\n" );
	print( "-h <High coverage cutoff for reference genotype; default = $h_DEF>\n" );
	print( "-L <Low coverage cutoff for comparison data; default = $L_DEF>\n" );
	print( "-H <High coverage cutoff for comparison data; default = $H_DEF>\n" );
	print( "-c <input comparison table file>\n" );
	exit( 0 );
    }
    unless( defined($opt_b) ) {
	$opt_b = $b_DEF;
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

