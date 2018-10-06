#!/usr/bin/perl

use strict;
use Getopt::Long;

MAIN : {

    my $read1_bam;
    my $read2_bam;
    my $qual_limit = 0;
    my $single_flag;
    my $bed_file;
    my $two_bed;
    my $frag_file;
    my $defined_chr;

    GetOptions('file1=s' => \$read1_bam,
	       'file2=s' => \$read2_bam,
	       'qual=i' => \$qual_limit,
	       'single' => \$single_flag,
	       'bed_file=s' => \$bed_file,
	       'paired_bed=s' => \$two_bed,
	       'frag_file=s' => \$frag_file,
	       'chr_file=s' => \$defined_chr);
    
    my $error1 = "\nUsage: ./two_read_bam_combiner.pl --file1 read1_bam_file --file2 read2_bam_file\n\n";
    my $error2 = "Options:\n";
    my $error3 = "--qual N [Reads pairs with mapq less than N will be filtered out, default is 0]\n";
    my $error4 = "--single [Will also print single ended reads with mapq greater than quality limit]\n";
    my $error5 = "--bed_file filter_bed_file [Using locations provided in the bed file, will keep only pairs in locations in the bed file]\n";
    my $error6 = "--paired_bed two_way_bed_file [format: chr1 start1 end1 chr2 start2 end2 - this will filter out any reads mapping between these loci]\n";
    my $error7 = "--frag_file fragment_file [Will add a label to each read for the fragment and its pair]\n";
    my $error8 = "--chr_file defined_chr_file [Will only considering chr in this file as mapped]\n";

    if ((not defined $read1_bam) ||
	(not defined $read2_bam)) {
	die ($error1,$error2,$error3,$error4,$error5,$error6,$error7,$error8,"\n");
    }

    if ($single_flag) {
	print STDERR "Will print read pairs and single reads with mapping quality greater than $qual_limit\n";
    } else {
	print STDERR "Will print read pairs with mapping quality greater than $qual_limit\n";
    }

    my $hash;
    my $lookup;
    if (defined $bed_file) {
	print STDERR "Will filter out reads based on $bed_file\n";
	open(FILE,$bed_file);
	while (my $line = <FILE>) {
	    chomp $line;
	    my ($chr, $start, $end, $n) = split(/\t/,$line);
	    push(@{$hash->{$chr}},$start);
	    $lookup->{$chr}->{$start} = $end;
	}
	close(FILE);
    } else {
	print STDERR "No bed file specified for filtering\n";
    }
   
    my $def_chr_hash;
    if (defined $defined_chr) {
	print STDERR "Will only consider reads mapping to chromosomes in $defined_chr as mapped\n";
	open(FILE,$defined_chr);
	while (my $line = <FILE>) {
	    chomp $line;
	    my ($chr, $size) = split(/\t/,$line);
	    $def_chr_hash->{$chr} = 1;
	}
	close(FILE);
    }

    my $for_frag_hash;
    my $rev_frag_hash;
    my $for_frag_lookup;
    my $rev_frag_lookup;

    if (defined $frag_file) {
	open(FILE,$frag_file);
	
	while (my $line = <FILE>) {
	    chomp $line;
	    my ($chr, $start, $end) = split(/\t/,$line);
	    push(@{$for_frag_hash->{$chr}},$start);
	    push(@{$rev_frag_hash->{$chr}},$end);
	    $for_frag_lookup->{$chr}->{$start} = $end;
	    $rev_frag_lookup->{$chr}->{$end} = $start;
	}

	close(FILE);
    }

    open(FILE1,"samtools view -h $read1_bam |");
    open(FILE2,"samtools view $read2_bam |");

    my $line1 = <FILE1>;
    my $line2 = <FILE2>;

    until ($line1 !~ m/^\@/) {
	print $line1;
	$line1 = <FILE1>;
    }

    my $counter = 0;
    my $new_counter = 0;

    while (defined $line1) {
	$counter++;
	if ($counter == ($new_counter + 100000)) {
	    print STDERR $counter . "\n";
	    $new_counter = $counter;
	}
	chomp $line1;
	chomp $line2;
	
	my ($id1, $flag1, $chr_from1, $loc_from1, $mapq1, $cigar1, $d1_1, $d2_1, $d3_1, $read1, $read_qual1, @rest1) = split(/\t/,$line1);
	my ($id2, $flag2, $chr_from2, $loc_from2, $mapq2, $cigar2, $d1_2, $d2_2, $d3_2, $read2, $read_qual2, @rest2) = split(/\t/,$line2);

	my ($new_id1, $new_id2) = id_test($id1,$id2);
       
	my $bin1 = reverse(dec2bin($flag1));
	my $bin2 = reverse(dec2bin($flag2));

	my @binary1 = split(//,$bin1);
	my @binary2 = split(//,$bin2);

	my $trouble = 0;
	if (($binary1[2] == 1) && ($mapq1 >= 10)) {
	    $trouble = 1;
	}
	if (($binary2[2]== 1) && ($mapq2 >= 10)) {
            $trouble = 1;
        }

	my $qual_test1 = 0;
	my $qual_test2 = 0;
	if ($qual_limit == 0) {
	    if ($mapq1 == 0) {
		$qual_test1 = 1;
	    }
	    if ($mapq2 == 0) {
		$qual_test2 = 1;
	    }
	} else {
	    if ($mapq1 < $qual_limit) {
		$qual_test1 = 1;
	    }
	    if ($mapq2 < $qual_limit) {
		$qual_test2 = 1;
            }
	}

	my $def_chr1 = 0;
	my $def_chr2 = 0;
	if ($chr_from1 eq "*") {
	    $def_chr1 = 1;
	}
	if ($chr_from2 eq "*") {
	    $def_chr2 = 1;
	}
	if (defined $defined_chr) {
	    if (not defined $def_chr_hash->{$chr_from1}) {
		$def_chr1 = 1;
	    }
	    if (not defined $def_chr_hash->{$chr_from2}) {
		$def_chr2 = 1;
            }	    
	}

	my $proper_pair1;
	my $proper_pair2;
	my $dist1;
	my $dist2;

	if (($binary1[2] == 0) && ($binary2[2] == 0)) {
	    $proper_pair1 = 1;
	    $proper_pair2 = 1;
	    if ($chr_from1 eq $chr_from2) {
		my $dist = abs($loc_from1 - $loc_from2);
		if ($loc_from1 >= $loc_from2) {
		    $dist1 = -1*$dist;
		    $dist2 = $dist;
		} else {
		    $dist1 = $dist;
		    $dist2 = -1*$dist;
		}
	    } else {
		$dist1 = 0;
		$dist2 = 0;
	    }
	} else {
	    $proper_pair1 = 0;
            $proper_pair2 = 0;
	    $dist1 = 0;
	    $dist2 = 0;
	}

	my $out_chr1;
	my $out_chr2;
	if ($chr_from1 eq $chr_from2) {
	    if (($chr_from1 ne "*") &&
		($chr_from2 ne "*")) {
		$out_chr1 = "=";
		$out_chr2 = "=";
	    } else {
		$out_chr1 = "*";
		$out_chr2 = "*";
	    }
	} else {
	    $out_chr1 = $chr_from2;
	    $out_chr2 = $chr_from1;
	}
	
       
	my $new_bin1 = join("","000000000000000000000",$binary1[10],$binary1[9],$binary1[8],"0","1",$binary2[4],$binary1[4],$binary2[2],$binary1[2],$proper_pair1,"1");
	my $new_bin2 = join("","000000000000000000000",$binary2[10],$binary2[9],$binary2[8],"1","0",$binary1[4],$binary2[4],$binary1[2],$binary2[2],$proper_pair2,"1");

	my $new_flag1 = bin2dec($new_bin1);
	my $new_flag2 = bin2dec($new_bin2);

	unless ($trouble > 0) {

	    if (($qual_test1 == 0) &&
		($qual_test2 == 0) && 
		($def_chr1 == 0) &&
		($def_chr2 == 0)) {

		my $test_sum = 0;

		if (defined $bed_file) {

		    my $test1 = bin_search($chr_from1,$loc_from1,$hash,$lookup);
		    my $test2 = bin_search($chr_from2,$loc_from2,$hash,$lookup);

		    $test_sum += $test1;
		    $test_sum += $test2;

		}

		if ($test_sum == 0) {

		    if (defined $frag_file) {
			my $frag1 = frag_search($chr_from1,$loc_from1,$binary1[4],$cigar1,$for_frag_hash,$rev_frag_hash,$for_frag_lookup,$rev_frag_lookup);
			my $frag2 = frag_search($chr_from2,$loc_from2,$binary2[4],$cigar2,$for_frag_hash,$rev_frag_hash,$for_frag_lookup,$rev_frag_lookup);
			print(join("\t",$new_id1,$new_flag1,$chr_from1,$loc_from1,$mapq1,$cigar1,$out_chr1,$loc_from2,$dist1,$read1,$read_qual1,@rest1,"XI:Z:" . $frag1,"XA:Z:" . $frag2) . "\n");
			print(join("\t",$new_id2,$new_flag2,$chr_from2,$loc_from2,$mapq2,$cigar2,$out_chr2,$loc_from1,$dist2,$read2,$read_qual2,@rest2,"XI:Z:" . $frag2,"XA:Z:" . $frag1) . "\n");
		    } else {
		    
			print(join("\t",$new_id1,$new_flag1,$chr_from1,$loc_from1,$mapq1,$cigar1,$out_chr1,$loc_from2,$dist1,$read1,$read_qual1,@rest1) . "\n");
			print(join("\t",$new_id2,$new_flag2,$chr_from2,$loc_from2,$mapq2,$cigar2,$out_chr2,$loc_from1,$dist2,$read2,$read_qual2,@rest2) . "\n");

		    }

		}

	    } elsif (($qual_test1 == 0) && ($qual_test2 > 0) && ($single_flag) && ($def_chr1 == 0)) {

		my $test_sum = 0;

		if (defined $bed_file) {

		    my $test1 = bin_search($chr_from1,$loc_from1,$hash,$lookup);

		    $test_sum += $test1;
		    
		}

		if ($test_sum == 0) {

		    if (defined $frag_file) {
                        my $frag1 = frag_search($chr_from1,$loc_from1,$binary1[4],$cigar1,$for_frag_hash,$rev_frag_hash,$for_frag_lookup,$rev_frag_lookup);	
			print $line1 . "\tXI:Z:" . $frag1 . "\n";  
		    } else {

			print $line1 . "\n";

		    }

		}

	    } elsif (($qual_test1 > 0) && ($qual_test2 == 0) && ($single_flag) && ($def_chr2 == 0)) {
		
		my $test_sum = 0;

		if (defined $bed_file) {

		    my $test2 = bin_search($chr_from2,$loc_from2,$hash,$lookup);

		    $test_sum += $test2;
		    
		} 

		if ($test_sum == 0) {

		    if (defined $frag_file) {
			my $frag2 = frag_search($chr_from2,$loc_from2,$binary2[4],$cigar2,$for_frag_hash,$rev_frag_hash,$for_frag_lookup,$rev_frag_lookup);
			print $line2 . "\tXI:Z:" . $frag2 . "\n"; 
		    } else {

			print $line2 . "\n";
		    }
		}

	    }

	}

	$line1 = <FILE1>;
	$line2 = <FILE2>;

    }

}

sub frag_search {
    
    my $chr_from = $_[0];
    my $loc_from = $_[1];
    my $bin = $_[2];
    my $cigar = $_[3];
    my $for_frag_hash = $_[4];
    my $rev_frag_hash = $_[5];
    my $for_frag_lookup = $_[6];
    my $rev_frag_lookup = $_[7];

    if ($bin == 0) {
	my $i = 0;
	my $j = scalar(@{$for_frag_hash->{$chr_from}}) - 1;
	my $mid = int(($i + $j)/2);
	while (($j - $i) > 1) {
	    if ($loc_from > $for_frag_hash->{$chr_from}->[$mid]) {
		$i = $mid;
		$mid = int(($i + $j)/2);
	    } else {
		$j = $mid;
		$mid = int(($i + $j)/2);
	    }
	}
	my $start;
	if ($loc_from < $for_frag_hash->{$chr_from}->[$j]) {
	    $start = $for_frag_hash->{$chr_from}->[$i];
	} else {
	    $start = $for_frag_hash->{$chr_from}->[$j];	    
	}
	my $end = $for_frag_lookup->{$chr_from}->{$start};
	my $out_frag = $chr_from . "," . $start . "," . $end;
	return($out_frag);
    } else {
	my $i = 0;
        my $j = scalar(@{$rev_frag_hash->{$chr_from}}) - 1;
        my $mid = int(($i + $j)/2);
        while (($j - $i) > 1) {
            if (($loc_from + get_read_length_cigar($cigar)) > $rev_frag_hash->{$chr_from}->[$mid]) {
                $i = $mid;
                $mid = int(($i + $j)/2);
            } else {
                $j = $mid;
                $mid = int(($i + $j)/2);
            }
        }
	my $end;
        if (($loc_from + get_read_length_cigar($cigar))< $rev_frag_hash->{$chr_from}->[$i]) {
            $end = $rev_frag_hash->{$chr_from}->[$i];
        } else {
            $end = $rev_frag_hash->{$chr_from}->[$j];
        }
        my $start = $rev_frag_lookup->{$chr_from}->{$end};
        my $out_frag = $chr_from . "," . $start . "," . $end;
        return($out_frag);
    }

}

sub get_read_length_cigar {

    my $cigar = $_[0];
    my @cig_array = ($cigar =~ m/\d+D+/g);
    my $sum = 0;
    foreach my $cig (@cig_array) {
	if ($cig =~ m/[MDN]/) {
	    chop $cig;
	    $sum += $cig;
	}
    }
    return($sum);

}

sub id_test {

    my $id1 = $_[0];
    my $id2 = $_[1];

    $id1 =~ s/\/[12]$//g;
    $id2 =~ s/\/[12]$//g;

    if ($id1 ne $id2) {
	die ("The IDs don't match up\n");
    } else {
	return($id1,$id2);
    }

}

sub bin_search {

    my $chr = $_[0];
    my $loc = $_[1];
    my $hash = $_[2];
    my $lookup = $_[3];

    if ((not defined $hash->{$chr}) ||
	(not defined $lookup->{$chr})) {
	return(1);
    } else {
	my $i = 0;
	my $j = scalar(@{$hash->{$chr}}) - 1;
	my $mid = int(($i + $j)/2);
	while (($j - $i) > 1) {
	    if ($loc < $hash->{$chr}->[$mid]) {
		$j = $mid;
		$mid = int(($i + $j)/2);
	    } else {
		$i = $mid;
		$mid = int(($i + $j)/2);
	    }
	}
	if (($loc >= $hash->{$chr}->[$i]) &&
	    ($loc <= $lookup->{$chr}->{$hash->{$chr}->[$i]})) {
	    return(1);
	} else {
	    return(0);
	}

    }

}

sub dec2bin {

    my $str = unpack("B32", pack("N", shift));
    return $str;

}

sub bin2dec {
    return unpack("N", pack("B32", substr("0" x 32 . shift, -32)));
}
    
