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
  if($seq_counter >= 0) {
    chomp;
    my @tmp=split(',',$_);
    my $cdr3 = $tmp[$cdr3_col];$cdr3 =~ s/^\s+//;
    my $seqnum = int $tmp[$num_col]; #id
    my $header = $tmp[$header_col]; $header =~ s/^\s+//;
    my $sequence = $tmp[$seq_col]; $sequence =~ s/^\s+//;
    my $cdr1 = $tmp[$cdr3_col+1];
    my $cdr2 = $tmp[$cdr3_col+2];

    my $dirname = "seq_" . $seqnum;
    print "$dirname $seqnum $header $cdr3 $sequence $cdr1 $cdr2\n";


    if(-e "$dirname/nb_loop_0.pdb") {
      chdir $dirname;
      if(not (-e "soap_score4.res" and (-s "soap_score4.res" > 200))) {
          writeCDR3File($cdr3, $sequence);
          writeFrameFile($cdr1, $cdr2, $cdr3, $sequence);
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

sub writeCDR3File {
    my $cdr3 = shift;
    my $sequence = shift;
    my $cdr3_start = index($sequence, $cdr3);
    my $cdr3_end = $cdr3_start + length($cdr3);
    $cdr3_start +=1;
    print "CDR3 loop: $cdr3_start $cdr3_end\n";
    open OUT, ">mycdrs3";
    for(my $i=$cdr3_start; $i<$cdr3_end; $i++) {
        print OUT "$i H\n";
    }
    close OUT;
}

sub writeFrameFile {
    my $cdr1 = shift;
    my $cdr2 = shift;
    my $cdr3 = shift;
    my $sequence = shift;

    my $cdr1_start = index($sequence, $cdr1);
    my $cdr1_end = $cdr1_start + length($cdr1);
    $cdr1_start +=1;

    my $cdr2_start = index($sequence, $cdr2);
    my $cdr2_end = $cdr2_start + length($cdr2);
    $cdr2_start +=1;

    my $cdr3_start = index($sequence, $cdr3);
    my $cdr3_end = $cdr3_start + length($cdr3);
    $cdr3_start +=1;

    open OUT, ">myframe";
    for(my $i=1; $i<$cdr1_start-3; $i++) { print OUT "$i H\n"; }
    for(my $i=$cdr1_end+3; $i<$cdr2_start-3; $i++) { print OUT "$i H\n"; }
    for(my $i=$cdr2_end+3; $i<$cdr3_start-3; $i++) { print OUT "$i H\n"; }
    for(my $i=$cdr3_end+3; $i<length($sequence); $i++) { print OUT "$i H\n"; }
    close OUT;
}
