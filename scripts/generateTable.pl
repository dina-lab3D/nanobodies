#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";


if ($#ARGV != 4) {
  print "Usage: generateTable.pl <nanobodies_csv_file> <num_col> <header_col> <seq_col> <cdr3_col>\n";
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
      chdir $folder or die "No such folder $folder\n";
      $folder_counter++;
  }
  if($seq_counter >= 0) {
    chomp;
    my @tmp=split(',',$_);
    my $cdr3 = $tmp[$cdr3_col];$cdr3 =~ s/^\s+//;
    my $seqnum = int $tmp[$num_col]; #id
    my $header = $tmp[$header_col]; $header =~ s/^\s+//;
    my $sequence = $tmp[$seq_col];$sequence =~ s/^\s+//;
    my $dirname = "seq_" . $seqnum;
#print "$dirname $dirname/soap_score4.res\n";
    if(-e "$dirname/soap_score4.res" and (-s "$dirname/soap_score4.res" > 200) and -e "$dirname/epi.csv") {
        chdir $dirname;
        print "$dirname $seqnum $header $cdr3 $sequence\n";
        if(not -e "top10.txt") {
            `grep "|" soap_score?.res | grep -v SOAP | awk '{print \$1" "\$2" "\$4}' | sort -nk3 | head -10 > top10.txt`;
        }
        my $soap = `awk '{sum+=\$3} END {print sum/NR}' top10.txt`; chomp $soap;
        my $cdr3len = length $cdr3;
        my $cmd =  "echo $_, $dirname, $cdr3len, $soap, > tmp";
        `$cmd`;
        $cmd = "cat epi.csv >> tmp";
        `$cmd`;
        `tr '\n' ' ' < tmp > stat.csv`;
        `sed -i -e '\$a\\' stat.csv`;
        #`echo \n >> stat.csv`;
        #print "$dirname $header $cdr3 $cdr3len $cmd1 $cmd10 $cmd100 $cmd1000  $label\n";

        `rm -f soap_score*.pdb log tmp 1n5u.pdb`;
        `cat stat.csv >> ../../epi_stat.csv`;
        chdir "..";
    } else {
        print "Docking not done $dirname\n";
    }

  }
  $seq_counter++;
}
