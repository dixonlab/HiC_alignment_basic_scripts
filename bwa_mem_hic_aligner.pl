#!/usr/bin/perl

use strict;

MAIN : {

    my ($fastq_file, $ref, $threads) = @ARGV;
    if ((not defined $fastq_file) ||
	(not defined $ref) ||
	(not defined $threads)) {
	die ("Usage: ./bwa_mem_hic_aligner.pl <fastq file> <reference genome> <threads>\n");
    }

    my $root = $fastq_file;
    $root =~ s/\.fastq.+$//g;
    
    open(FILE,"bwa mem -t $threads $ref $fastq_file |");
    
    my $line = <FILE>;
    my $next_line = <FILE>;
    chomp $line;
    chomp $next_line;

    while ($line =~ m/^\@/) {
	print $line . "\n";
	$line = $next_line;
	$next_line = <FILE>;
	chomp $next_line;
    }
    
    my ($id1) = split(/\t/,$line);
    my ($id2) = split(/\t/,$next_line);
	
    while () {

	my @array;
	push(@array,$line);

	my $quit_flag = 0;
	
	while ($id1 eq $id2) {
	    push(@array,$next_line);
	    $line = $next_line;
	    $next_line = <FILE>;
	    if ((not defined $line) ||
		(not defined $next_line)) {
		$quit_flag = 1;
		last;
	    }
	    chomp $next_line;
	    ($id1) =split(/\t/,$line);
	    ($id2) = split(/\t/,$next_line);
	}

	if ($quit_flag == 1) {
	    last;
	}
       
	if (scalar(@array) == 1) {
	    print $array[0] . "\n";
	} else {
	    my $hash;
	    for (my $i = 0; $i < scalar(@array); $i++) {
		my @read_array = split(/\t/,$array[$i]);
		my $dist = get_dist_from_cigar($read_array[1],$read_array[5]);
		if (not defined $hash->{$dist}) {
		    $hash->{$dist} = $i;
		} else {
		    die ("I am assuming this doesn't happen\n");
		}
	    }
	    
	    my @key_array = sort {$a <=> $b} keys %$hash;
	    my @read_array = split(/\t/,$array[$hash->{$key_array[0]}]);
	    my @bin = split(//,reverse(dec2bin($read_array[1])));
	    $bin[11] = 0;
	    $read_array[1] = bin2dec(join("",reverse(@bin)));
	    print(join("\t",@read_array) . "\n");

	}

	$line = $next_line;
	$next_line = <FILE>;
	if ((not defined $line) ||
	    (not defined $next_line)) {
	    last;
	}
	chomp $next_line;
	($id1) =split(/\t/,$line);
	($id2) = split(/\t/,$next_line);
    }
    close(FILE);

}

sub get_dist_from_cigar {
    
    my $flag = $_[0];
    my $cigar = $_[1];

    my @bin = split(//,reverse(dec2bin($flag)));

    my @cig_array = ($cigar =~ m/\d+\D+/g);

    my $sum = 0;

    if ($bin[4] == 0) {
	my $i = 0;
	until (($cig_array[$i] =~ m/M/) ||
	       ($i == (scalar(@cig_array) - 1))) {
	    if ($cig_array[$i] =~ m/[HS]/) {
		my $cig = $cig_array[$i];
		chop $cig;
		$sum += $cig;
	    }
	    $i++;
	}
    } else {
	my $i = scalar(@cig_array) - 1;
	until (($cig_array[$i] =~ m/M/) ||
	       ($i == 0)) {
	    if ($cig_array[$i] =~ m/[HS]/) {
		my $cig = $cig_array[$i];
		chop $cig;
		$sum +=$cig;
            }
	    $i += -1;
	}
    }

    return($sum);

}

sub dec2bin {

    my $str = unpack("B32", pack("N", shift));
    return $str;

}

sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
