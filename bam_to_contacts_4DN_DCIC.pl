#!/usr/bin/perl

use strict;

MAIN : {

    my ($bam_file, $prefix) = @ARGV;
    if ((not defined $bam_file) ||
	(not defined $prefix)) {
	die ("Usage: ./bam_to_contacts_4DN_DCIC.pl <bam file (must be name sorted)> <prefix [for temp files]>\n");
    }
       
    print "## pairs format v1.0\n";
    print "#columns: readID chr1 position1 chr2 position2 strand1 strand2\n";
    
    my $hash;
    my $counter = 0;

    my @chr_array;
    my $chr_keys;

    my $t1 = time();

    my $file_hash;

    open(FILE,"samtools view -H $bam_file |");
    while (my $liner = <FILE>) {
	chomp $liner;

	if ($liner =~ m/^@/) {

	    
	    if ($liner =~ m/^\@SQ/) {
		my ($id, $chr_name, $size) = split(/\t/,$liner);
		$chr_name =~ s/SN://g;
                if (($chr_name !~ m/random/) &&
                    ($chr_name !~ m/chrUn/) &&
                    ($chr_name !~ m/decoy/) &&
                    ($chr_name !~ m/EBV/) &&
		    ($chr_name !~ m/lambda/)) {
                    push(@chr_array,$chr_name);
                }
	    }
	}
    }
    close(FILE);    

    if (not defined $chr_keys->{"chr1"}) {
	if (scalar(@chr_array) < 1) {
	    die ("There are no chromosomes in the header file\n");
	}
	for (my $i = 0; $i < scalar(@chr_array); $i++) {
	    $chr_keys->{$chr_array[$i]} = $i;
	    for (my $j = $i; $j < scalar(@chr_array); $j++) {
#		print STDERR $i . "\t" . $j . "\n";
		my $name = "$i\_$j";
		open($file_hash->{$name},">$prefix.$name.txt") || die ("Can't open file\n");
	    }
	}

    }

#    exit;

    open(FILE,"samtools view $bam_file |");
    my $line1 = <FILE>;
    chomp $line1;

    while () {

	unless ($counter % 1000000) {
	    print STDERR $counter . "\n";
	    my $val = scalar(keys %$hash);
	}
	$counter++;

	my ($id1, $flag1, $chr_from1, $loc_from1, $mapq1, $cigar1, $chr_to1, $loc_to1) = split(/\t/,$line1);
#	print STDERR $chr_from1 . "\t" . $chr_to1 . "\n";
	if (($chr_from1 =~ m/random/) ||
	    ($chr_from1 =~ m/chrUn/) ||
	    ($chr_from1 =~ m/decoy/) ||
	    ($chr_from1 =~ m/EBV/) ||
	    ($chr_from1 =~ m/lambda/) ||
	    ($chr_to1 =~ m/random/) ||
	    ($chr_to1 =~ m/chrUn/) ||
	    ($chr_to1 =~ m/decoy/) ||
	    ($chr_to1 =~ m/EBV/) ||
	    ($chr_to1 =~ m/lambda/)) {
	    $line1 = <FILE>;
	    chomp $line1;
	    if (not defined $line1) {
		last;
	    }
	    next;
	}
	
	if ($chr_to1 eq "*") {
	    $line1 = <FILE>;
            chomp $line1;
	    if (not defined $line1) {
		last;
	    }
	    next;
	}
	
	my $chr1 = $chr_from1;
	my $pos1;
	my $strand1;
	my $b1 = reverse(dec2bin($flag1));
	my $st1 = substr($b1,4,1);

	if ($st1 == 0) {
	    $pos1 = $loc_from1;
	    $strand1 = "+";
	} else {
	    my @array = ($cigar1 =~ m/\d+\D/g);
	    my $length = 0;
	    foreach my $val (@array) { 
		if ($val =~ m/[MD]/) {
		    chop $val; 
		    $length += $val;
		}
	    }
		    
	    $pos1 = $loc_from1 + $length;
	    $strand1 = "-";
	}

	my $line2 = <FILE>;

	if (not defined $line2) {
	    last;
	}

	chomp $line2;
	my ($id2, $flag2, $chr_from2, $loc_from2, $mapq2, $cigar2, $chr_to2, $loc_to2) = split(/\t/,$line2);

	if ($id1 ne $id2) {
	    $line1 = $line2;
	    next;
#	    print STDERR $id1 . "\t" . $id2 . "\n";
#	    die ("This IDs don't match!\n");
	}

	my $chr2 = $chr_from2;
	my $pos2;
	my $strand2;
	my $b2 = reverse(dec2bin($flag2));
        my $st2 = substr($b2,4,1);

	if ($st2 == 0) {
	    $pos2 = $loc_from2;
	    $strand2 = "+";
	} else {
	    my @array = ($cigar2 =~ m/\d+\D/g);
	    my $length = 0;
	    foreach my $val (@array) {
		if ($val =~ m/[MD]/) {
		    chop $val;
		    $length += $val;
		}
	    }
		    
	    $pos2 = $loc_from2 + $length;
	    $strand2 = "-";
	}

	my $n1 = $chr_keys->{$chr1};
	my $n2 = $chr_keys->{$chr2};
		
	if ($chr1 eq $chr2) {
	    my $out_name = "$n1\_$n2";
	    if (abs($pos1 - $pos2) >= 1000) {
		if ($pos1 < $pos2) {		    
		    print { $file_hash->{$out_name} } join("\t",$id1,$chr1,$pos1,$chr2,$pos2,$strand1,$strand2) . "\n";
		} else {

		    print { $file_hash->{$out_name} } join("\t",$id1,$chr2,$pos2,$chr1,$pos1,$strand2,$strand1) . "\n";
		}
	    }
	} else {
		    
	    if ($chr_keys->{$chr1} < $chr_keys->{$chr2}) {
		
		my $out_name = "$n1\_$n2";
	
		print { $file_hash->{$out_name} } join("\t",$id1,$chr1,$pos1,$chr2,$pos2,$strand1,$strand2) . "\n";

	    } else {

		my $out_name = "$n2\_$n1";

#		print STDERR $out_name . "\n";
#		print STDERR $chr1 . "\t" . $chr2 . "\n";

		if (defined $file_hash->{$out_name}) {
			
		    print { $file_hash->{$out_name} } join("\t",$id1,$chr2,$pos2,$chr1,$pos1,$strand2,$strand1) . "\n";

		}
			
	    }
		
	}
	$line1 = <FILE>;
	chomp $line1;
	if (not defined $line1) {
	    last;
	}
    }
    close(FILE);

    for (my $i = 0; $i < scalar(@chr_array); $i++) {
	for (my $j = $i; $j < scalar(@chr_array); $j++) {
	    my $name = "$i\_$j";
	    close($file_hash->{$name});

	    open(FILE,"<$prefix.$name.txt") || die ("Can't open file\n");
	    while (my $line = <FILE>) {
		chomp $line;
		print $line . "\n";
	    }
	    close(FILE);
	    
	    system("rm $prefix.$name.txt");

	}
    }

    my $t2 = time();
    my $interval = $t2 - $t1;
    print STDERR "Took $interval seconds\n";

    
    
}

sub dec2bin {

    my $str = unpack("B32", pack("N", shift));
    return $str;

}
