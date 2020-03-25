#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);
#use File::Basename qw();
my ($folder_name, $path, $suffix) = File::Basename::fileparse(Cwd::cwd());
print "$folder_name\n";

my $home = "$FindBin::Bin";
require "$home/Util.pm";

if ($#ARGV >= 2) {
  print "Usage: generateModelsFromTop10.pl <N=10> <prefix>\n";
  exit;
}

my $N = 10;
if($#ARGV >= 0) { $N = $ARGV[0]; }
my $prefix = "";
if($#ARGV == 1) { $prefix = $ARGV[1]; }


my $soapfiles = $prefix . "soap_score0.res " . $prefix . "soap_score1.res " . $prefix . "soap_score2.res " . $prefix . "soap_score3.res " . $prefix . "soap_score4.res";
my $cmd = `grep "|" $soapfiles | grep -v SOAP | awk '{print \$1" "\$2" "\$4}' | sort -nk3 | head -$N`;
#print "grep \"|\" $soapfiles | grep -v SOAP | awk '{print \$1\" \"\$2\" \"\$4}' | sort -nk3 | head -$N";

my $counter = 1;
my $res_files = "";
my $soap = 0;
my $firstPDB = $folder_name . "_" . $prefix . "_1l.pdb";
foreach (split(/\n/,$cmd)) {
  #print "new $_\n";
  my @tmp = split (' ', $_);
  chop $tmp[0];
  my $soap_file = $tmp[0];
  my $res_num = $tmp[1];
  my $score = $tmp[2];
  my $cmdt = "/cs/staff/dina/scripts/transOutput.pl $soap_file $res_num $res_num l";
  `$cmdt`;
  print "$cmdt\n";
  my $pdbfile = $soap_file . "." . $res_num . ".pdb";
  if($counter > 1) { # transform to first to account for symmetry
    `pdb_trans 1.373 0.3337 2.863 23.77 39.6 -37.76 < $pdbfile > tmp.pdb`;
    my $res1 = `rmsd $firstPDB $pdbfile | tail -n1`; chomp $res1; $res1 += 0.0;
    my $res2 = `rmsd $firstPDB tmp.pdb | tail -n1`;  chomp $res2; $res2 += 0.0;
    print "$res1 $res2\n";
    if($res2 < $res1) {
      `mv tmp.pdb $pdbfile`;
    }
  } else {
    `cp $pdbfile $firstPDB`;
  }
  my $newFileName = $folder_name . "_" . $prefix . "_" . $counter . ".pdb";
  `cp 1dug.pdb $newFileName`;
  `cat $pdbfile >> $newFileName`;
  $counter++;
  $res_files .= $newFileName;
  $res_files .= " ";
  $soap += $score;
}
#my $res_files = `ls *soap_score*.pdb`;
print "Res: $res_files\n";
my $precision = `~/scripts/rmsdAll.pl $res_files | awk '{sum+=\$3} END {print sum/NR}'`;
chomp $precision;
print "Precision = $precision\n";
$soap /= $counter;
print "SOAP = $soap\n";
$cmd = "$home/../epiDock/interface AB H 6.0 $res_files";
`$cmd`;
open OUT, ">stat.csv";
print OUT "$folder_name , $precision, $soap\n";
close OUT;
