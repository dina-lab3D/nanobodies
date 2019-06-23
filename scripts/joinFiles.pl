#!/usr/bin/perl -w

use strict;

if ($#ARGV != 6) {
  print "Usage: joinFiles.pl <file1> <file2> <file3> <header_col> <seq_col> <cdr3_col> <affinity_col>\n";
  exit;
}

my $sequence_file1 = $ARGV[0];
my $sequence_file2 = $ARGV[1];
my $sequence_file3 = $ARGV[2];

my $header_col = $ARGV[3];
my $seq_col = $ARGV[4];
my $cdr3_col = $ARGV[5];
my $affinity_col = $ARGV[6];

my %map_id1;
my %map_id2;
my %map_id3;

open(DATA, $sequence_file1);
while(my $line = <DATA>) {
    chomp $line;
    my @tmp=split(',',$line);
    my $header = $tmp[$header_col]; #id
    if ($header ne "id") { $map_id1{$header} = $line; }
}
close DATA;

open(DATA, $sequence_file2);
while(my $line = <DATA>) {
    chomp $line;
    my @tmp=split(',',$line);
    my $header = $tmp[$header_col]; #id
    if ($header ne "id") { $map_id2{$header} = $line; }
}
close DATA;

open(DATA, $sequence_file3);
while(my $line = <DATA>) {
    chomp $line;
    my @tmp=split(',',$line);
    my $header = $tmp[$header_col]; #id
    if ($header ne "id") { $map_id3{$header} = $line; }
}
close DATA;

my $counter = 0;

foreach my $key (keys %map_id1) {
    # check if exists in 2 or 3
    my $affinity2 = -1;
    my $affinity3 = -1;
    my @tmp = split(',', $map_id1{$key});
    my $affinity1 = $tmp[$affinity_col];
    my $seq = $tmp[$seq_col];
    my $cdr3 = $tmp[$cdr3_col];

    if(exists $map_id2{$key}) {
      my @tmp = split(',', $map_id2{$key});
      $affinity2 = $tmp[$affinity_col];
    }
    if(exists $map_id3{$key}) {
      my @tmp = split(',', $map_id3{$key});
      $affinity3 = $tmp[$affinity_col];
    }

    print "$counter,$key,$seq,$cdr3,$affinity1,$affinity2,$affinity3 \n";
    $counter++;
}

foreach my $key (keys %map_id2) {
    # check if exists in 1 or 3
    my $affinity3 = -1;
    my @tmp = split(',', $map_id2{$key});
    my $affinity2 = $tmp[$affinity_col];
    my $seq = $tmp[$seq_col];
    my $cdr3 = $tmp[$cdr3_col];

    if(exists $map_id1{$key}) { next; }
    if(exists $map_id3{$key}) {
      my @tmp = split(',', $map_id3{$key});
      $affinity3 = $tmp[$affinity_col];
    }

    print "$counter,$key,$seq,$cdr3,-1,$affinity2,$affinity3 \n";
    $counter++;
}

foreach my $key (keys %map_id3) {
    # check if exists in 1 or 2
    if(exists $map_id1{$key}) { next; }
    if(exists $map_id2{$key}) { next; }

    my @tmp = split(',', $map_id3{$key});
    my $affinity3 = $tmp[$affinity_col];
    my $seq = $tmp[$seq_col];
    my $cdr3 = $tmp[$cdr3_col];

    print "$counter,$key,$seq,$cdr3,-1,-1,$affinity3 \n";
    $counter++;
}
