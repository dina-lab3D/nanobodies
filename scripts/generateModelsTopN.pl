#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";

if ($#ARGV >= 2) {
  print "Usage: generateModelsFromTop10txt.pl <N=10> <prefix>\n";
  exit;
}

my $N = 10;
if($#ARGV == 0) { $N = $ARGV[0]; }
my $prefix = "";
if($#ARGV == 1) { $prefix = $ARGV[1]; }


my $soapfiles = $prefix . "soap_score0.res " . $prefix . "soap_score1.res " . $prefix . "soap_score2.res " . $prefix . "soap_score3.res " . $prefix . "soap_score4.res";
my $cmd = `grep "|" $soapfiles | grep -v SOAP | awk '{print \$1" "\$2" "\$4}' | sort -nk3 | head -$N`;

foreach (split(/\n/,$cmd)) {
  #print "new $_\n";
  my @tmp = split (' ', $_);
  chop $tmp[0];
  my $soap_file = $tmp[0];
  my $res_num = $tmp[1];
  my $cmdt = "/cs/staff/dina/scripts/transOutput.pl $soap_file $res_num $res_num l";
  `$cmdt`;
  print "$cmdt\n";
}
my $res_files = `ls *soap_score*.pdb`;
print "Res: $res_files\n";
