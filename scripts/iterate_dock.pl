#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";


if ($#ARGV != 4) {
  print "Usage: iterate_dock.pl <nanobodies_csv_file> <num_col> <header_col> <seq_col> <cdr3_col>\n";
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

open(DATA, $ARGV[0]);
while(<DATA>) {

  # run folder management - up to 1,000 in a folder
  my $rem = $seq_counter % 1000;
  if($rem == 0) {
      if($seq_counter > 0) { # not first one
          chdir "..";
      }
      my $folder = "run_" . Util::get_counter_string_limit($folder_counter, 1000);
      print "$folder \n";
      mkdir $folder;
      chdir $folder;
      $folder_counter++;
  }
  if($seq_counter > 0) {
    chomp;
    my @tmp=split(',',$_);
    my $cdr3 = $tmp[$cdr3_col];
    my $seqnum = int $tmp[$num_col]; #id
    my $header = $tmp[$header_col]; #id
    my $sequence = $tmp[$seq_col];
    my $dirname = "seq_" . $seqnum;
    print "$dirname $seqnum $header $cdr3 $sequence\n";

    if(-e "$dirname/nb_loop_0.pdb") {
      chdir $dirname;
      if(not (-e "soap_score4.res" and (-s "soap_score4.res" > 200))) {
          print "Missing docking $dirname\n";
          my $cmd = "$home/dock.pl $dirname";
          print "$cmd\n";
          `$cmd`;
           $job_counter++;
          my $hours3 = 3600*5;
          if($job_counter > 500) { $job_counter = 0; sleep($hours3); }

      } else {
          print "Docking done $dirname\n";
      }
      chdir "..";
    }
  }
  $seq_counter++;
}
