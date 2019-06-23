#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";

my $modeller_home = "/cs/staff/dina/modeller9.18/bin";

if ($#ARGV != 4) {
  print "Usage: iterate_model_nanobodies.pl <nanobodies_csv_file> <num_col> <header_col> <seq_col> <cdr3_col>\n";
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
    my $cdr3 = $tmp[$cdr3_col]; $cdr3 =~ s/^\s+//;
    my $seqnum = int $tmp[$num_col]; #id
    my $header = $tmp[$header_col]; $header =~ s/^\s+//;
    my $sequence = $tmp[$seq_col]; $sequence =~ s/^\s+//;
    my $dirname = "seq_" . $seqnum;
    print "$dirname $seqnum $header $cdr3 $sequence\n";
    mkdir $dirname; # or die "Can't make directory $dirname\n";
    chdir $dirname;

    # recovery
    if(not (-e "NANO.BL00990001.pdb" || -e "nb_loop_0.pdb")) {
        print "Missing models $dirname\n";

      open OUT, ">seq";
      print OUT ">$header\n$sequence\n";
      close OUT;

      open OUT2, ">nano.ali";
      print OUT2 ">P1;NANO\n";
      print OUT2 "sequence:NANO:::::::0.00: 0.00\n";
      print OUT2 "$sequence*\n";
      close OUT2;

      my $cdr3_start = index($sequence, $cdr3);
      my $cdr3_end = $cdr3_start + length($cdr3);
      $cdr3_start +=1;
      print "CDR3 loop: $cdr3_start $cdr3_end\n";

      my $cmd = "$modeller_home/modpy.sh python $home/modelNanobody.py seq $cdr3_start $cdr3_end >& model.log";
      print "$cmd\n";
      my $currdir = cwd;
      open OUT3, ">script.sh";
      print OUT3 "#!/bin/tcsh\n";
      print OUT3 "cd $currdir\n";
      print OUT3 "$cmd\n";
      $cmd = "rm -f NANO.DL0* NANO.D0* NANO.V99* NANO.rsr NANO.IL00000001.pdb";
      print OUT3 "$cmd\n";


      `sbatch --mem=6000m -c1 --time=2:0:0 script.sh`;
       $job_counter++;
      my $hours5 = 3600*5;
      if($job_counter > 2000) { $job_counter = 0; sleep($hours5); }

    } else {
        print "Models done $dirname\n";
        #my $cmd = "grep NANO.B99 model.log | grep -v Open | tail -n1";
        #my $res = `$cmd`;
        #print $res;
    }
    chdir "..";
  }
  $seq_counter++;
}
