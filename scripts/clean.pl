#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";


if ($#ARGV != 4) {
  print "Usage: clean.pl <nanobodies_csv_file> <num_col> <header_col> <seq_col> <cdr3_col>\n";
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
  # folder management - up to 1,000 in a folder
  my $rem = $seq_counter % 1000;
  if($rem == 0) {
      if($seq_counter > 0) { # not first one
          chdir "..";
      }
      my $folder = "run_" . Util::get_counter_string_limit($folder_counter, 1000);
      print "$folder \n";
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

    if(-e "$dirname/NANO.BL00990001.pdb") {
        chdir $dirname;
        my $pdbs = `grep LOOP model.log | sort -nk4 | head -n5 | awk '{print \$2}' | tr '\n' ' '`;
        print $pdbs;
        my @list = split(' ', $pdbs);


        for (my $j=0; $j < 5; $j++) {
            my $loopfile = "nb_loop_" . $j . ".pdb";
            `cp $list[$j] $loopfile`;
        }

        `rm -f NANO.B*.pdb`;
        `rm -f default gnuplotfile NANO.sch NANO.lrsr NANO.ini slurm*`;
        chdir "..";
    }
  }
  $seq_counter++;
}
