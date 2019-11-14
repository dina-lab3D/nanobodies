#!/usr/bin/perl -w

use strict;

if ($#ARGV != 5) {
  print "Usage: cluster_analysis.pl <map_id_to_clusters.csv> <joined.csv> <num_col> <header_col> <seq_col> <cdr3_col> \n";
  exit;
}

my $cluster_file = $ARGV[0];

my $num_col = $ARGV[2];
my $header_col = $ARGV[3];
my $seq_col = $ARGV[4];
my $cdr3_col = $ARGV[5];

my %map_id_cluster;

open(DATA, $cluster_file);
while(my $line = <DATA>) {
    chomp $line;
    my @tmp=split(',',$line);
    my $id = int $tmp[0]; #id
    #$id =~ s/^\s+|\s+$//g;
    if ($id ne "id") { $map_id_cluster{$id} = $tmp[1]; }
}
close DATA;

my @neg_counter_array = (0.0,0.0,0.0,0.0,0.0);
my @pos_counter_array = (0,0,0,0,0);
my @rings_counter_array = (0,0,0,0,0);
my @cluster_size_counter_array = (0,0,0,0,0);
my @cdr3_length_counter_array = (0.0,0.0,0.0,0.0,0.0);
my @cdr3_cys_counter_array = (0.0,0.0,0.0,0.0,0.0);

open(DATA, $ARGV[1]);
while(my $line = <DATA>) {
    chomp $line;
    my @tmp=split(',',$line);
    my $id = $tmp[$num_col]; #id
    $id =~ s/^\s+|\s+$//g;

    if ($id ne "id") {
        if(exists $map_id_cluster{$id}) {
            my $cluster = int $map_id_cluster{$id};

            #print "$cluster,$line\n";

            $cluster_size_counter_array[$cluster-1]++;

            my $cdr3 = $tmp[$cdr3_col];
            #my $cys = int $tmp[5];
            #print "$cluster $cdr3 $cys\n";
            my $count_negative = () = $cdr3 =~ /D|E/g;
            my $count_positive = () = $cdr3 =~ /K|R/g;
            my $count_rings = () = $cdr3 =~ /W|Y/g;
            $neg_counter_array[$cluster-1] += $count_negative;
            $pos_counter_array[$cluster-1] += $count_positive;
            $rings_counter_array[$cluster-1] += $count_rings;
            $cdr3_length_counter_array[$cluster-1] += length $cdr3;
            #$cdr3_cys_counter_array[$cluster-1] += $cys;
        }
    }
}
print "Cluster# size -   +  net_charge total_charge WY cdr3_length \n";
for(my $i=0; $i<=4; $i++) {
    my $cl_num = $i+1;
    print "$cl_num $cluster_size_counter_array[$i]";
    my $neg_ratio = $neg_counter_array[$i]/$cluster_size_counter_array[$i];
    my $pos_ratio = $pos_counter_array[$i]/$cluster_size_counter_array[$i];
    my $net_charge = $pos_ratio-$neg_ratio;
    my $total_charge = $pos_ratio+$neg_ratio;
    my $rings_ratio = $rings_counter_array[$i]/$cluster_size_counter_array[$i];
    my $cdr3_length = $cdr3_length_counter_array[$i]/$cluster_size_counter_array[$i];
    my $cys_count = $cdr3_cys_counter_array[$i]/$cluster_size_counter_array[$i];
    #print " $neg_ratio $pos_ratio $net_charge $total_charge $rings_ratio $cdr3_length $cys_count\n";
    printf(" %.2f %.2f %.2f %.2f %.2f %.2f\n", $neg_ratio, $pos_ratio, $net_charge, $total_charge,$rings_ratio,$cdr3_length);
}
