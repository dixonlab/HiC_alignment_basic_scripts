#!/usr/bin/perl

use strict;


MAIN : {

    if (scalar(@ARGV) < 2) {
	die ("Usage: ./merge_contact_files.pl [at least two contact files]\n");
    }

    my $file_hash;

    for (my $i = 0; $i < scalar(@ARGV); $i++) {
	
	print STDERR "Going through " . $ARGV[$i] . " now\n";

	if ($ARGV[$i] =~ m/gz$/) {
	    open(FILE,"zcat $ARGV[$i] |");
	} else {
	    open(FILE,$ARGV[$i]);
	}
	
	while (my $line = <FILE>) {
	    chomp $line;
	    if ($line =~ m/\#/) {
		if ($i == 0) {
		    print $line . "\n";
		}
		next;
	    } 
	    my ($id, $chr1, $pos1, $chr2, $pos2, $str1, $str2) = split(/\t/,$line);

	    my $name = $chr1 . "_" . $chr2;

	    if (-e "$name.txt") {
		
	    } else {
		open($file_hash->{$name},">$name.txt") || die ("Can't open file\n");
	    }

	    print { $file_hash->{$name} } join("\t",$id,$chr1,$pos1,$chr2,$pos2,$str1,$str2) . "\n";

	}
    
	close(FILE);

    }

    foreach my $name (keys %$file_hash) {
	
	close($file_hash->{$name});

	open(FILE,"<$name.txt") || die ("Can't open file\n");
	while (my $line = <FILE>) {
	    chomp $line;
	    print $line . "\n";
	}
	close(FILE);
	    
	system("rm $name.txt");


    }

}
