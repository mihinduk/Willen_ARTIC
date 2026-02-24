#!/usr/bin/perl -w
# N_count_and_clean.pl
# 2023_09_19
use strict;

# This script will take A file of sequences, count the Ns and -s in the sequence, print a file ONLY containing sequences with <= 9000 Ns and -s

# $ARGV[0] = CLEANED fasta file
# $ARGV[1] = number of Ns
# $ARGV[2] = genome size

my ($N_count, $N, $desired, $seq_count, $final_count) = (0,0,0,0,0);
my ($name, $N_current);
my ($file, $rest) = split (/\./, scalar $ARGV[0], 2);
my $N_max = scalar $ARGV[1];
my $N_too_many = $N_max + 1;
my $outfile_desired = "$file". "_$N_max"."_N_or_less.fasta";
my @Ns = ();
my $genome_length = scalar $ARGV[2]; print "\$genome_length = $genome_length\n";
my $acceptable_length = $genome_length - $N_max; print "\$acceptable_length = $acceptable_length\n";

# Make the seq by seq file
open (IN, $ARGV[0]) or die "Couldn't open your cleaned fasta file $!";
open (OUT_desired, ">$outfile_desired") or die "Couldn't create your file of desired seqs $!";
while (my $line = <IN>) {
	next if ($line =~/^\n/);
	chomp ($line);
	if ($line=~/^\>/) { # This skips lines that start with >.
	$name = $line;
	$N_count = 0;
	$seq_count = 0;
	$final_count = 0;
	next;
	}
	elsif ($line !~ m/^>/) {
	++$N_count while $line =~m/N/ig;
	++$seq_count while $line =~m/[acgt]/ig; #print "\$seq_count = $seq_count\n";
	 if ($N_max < $N_count ){print "Too many Ns: $name\n"; next;}
	 $final_count = $seq_count + $N_max;
	 if ($final_count < $acceptable_length){ print "Too short: $name\t$final_count\n"; next;}
	 else {print OUT_desired "$name\n$line\n"; $desired ++;}
	 next;
	}
}
close (IN);
close (OUT_desired);


print "\nArs longa, vita brevis - Art is long, life is short - Hippocrates \nFinished counting your N's \n Sequences with $N_max or fewer Ns: $desired\nYour sequences are in: $outfile_desired\n";

