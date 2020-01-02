#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";


if ($#ARGV != 5) {
  print "Usage: iterate_dock.pl <nanobodies_csv_file1> <num_col> <header_col> <seq_col> <cdr3_col> <next_file>\n";
  exit;
}

my $sequence_file = $ARGV[0];
my $num_col = $ARGV[1];
my $header_col = $ARGV[2];
my $seq_col = $ARGV[3];
my $cdr3_col = $ARGV[4];


my $seq_counter = 0;
my $folder_counter = 1;
my $job_counter = 0;

my %map_id_folder;

open(DATA, $ARGV[0]);
while(<DATA>) {

  # run folder management - up to 1,000 in a folder
  my $rem = $seq_counter % 1000;
  if($rem == 0) {
      if($seq_counter > 0) { # not first one
          chdir "..";
      }
      my $folder = "run_" . Util::get_counter_string_limit($folder_counter, 1000);
      #print "$folder \n";
      mkdir $folder;
      chdir $folder;
      $folder_counter++;
  }
  if($seq_counter >= 0) {
    chomp;
    my @tmp=split(',',$_);
    my $cdr3 = $tmp[$cdr3_col];$cdr3 =~ s/^\s+//;
    my $seqnum = int $tmp[$num_col]; #id
    my $header = $tmp[$header_col]; $header =~ s/^\s+//;
    my $sequence = $tmp[$seq_col]; $sequence =~ s/^\s+//;

    my $dirname = "seq_" . $seqnum;
    $map_id_folder{$header} = $dirname;
    #print "$dirname $seqnum $header $cdr3 $sequence\n";
  }
  $seq_counter++;
}
close DATA;
chdir "..";


open(DATA2, $ARGV[5]);
while(<DATA2>) {
    chomp;
    my @tmp=split(',',$_);
    my $id1 = $tmp[0]; $id1 =~ s/^\s+//;
    my $id2 = $tmp[1]; $id2 =~ s/^\s+//;
    #print "$id1 $id2 : ";
    my $folder = "";
    if(exists $map_id_folder{$id1}) {
        $folder = $map_id_folder{$id1};
    } else {
        if(exists $map_id_folder{$id2}) {
            $folder = $map_id_folder{$id2};
        }
    }
    if (length($folder) > 0) {
        #print "Docking done $id1 $id2 $folder\n";
    } else {
        #print "Missing, $seq_counter, $id1, $id2, $tmp[2], $tmp[3], $tmp[4], $tmp[5]\n";
        print "$seq_counter, $_\n";
        $seq_counter++;
    }
}
close DATA2;
