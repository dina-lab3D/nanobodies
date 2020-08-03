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


my $soapfiles = $prefix . "soap_score0.res " . $prefix . "soap_score1.res " . $prefix . "soap_score2.res " . $prefix . "soap_score3.res " . $prefix . "soap_score4.res " . $prefix . "soap_score5.res " . $prefix . "soap_score6.res " . $prefix . "soap_score7.res " . $prefix . "soap_score8.res " . $prefix . "soap_score9.res";
my $cmd = `grep "|" $soapfiles | grep -v SOAP | awk '{print \$1" "\$2" "\$4}' | sort -nk3 | head -$N`;
#print "grep \"|\" $soapfiles | grep -v SOAP | awk '{print \$1\" \"\$2\" \"\$4}' | sort -nk3 | head -$N";

my $counter = 1;
my $res_files = "";
my $soap = 0;
foreach (split(/\n/,$cmd)) {
  #print "new $_\n";
  my @tmp = split (' ', $_);
  chop $tmp[0];
  my $soap_file = $tmp[0];
  my $res_num = $tmp[1];
  my $score = $tmp[2];
  my $cmdt = "/cs/staff/dina/scripts/transOutput.pl $soap_file $res_num $res_num ";
  `$cmdt`;
  print "$cmdt\n";
  my $pdbfile = $soap_file . "." . $res_num . ".pdb";
  my $newFileName = $folder_name . "_" . $prefix . "_" . $counter . ".pdb";
  `mv $pdbfile $newFileName`;
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
print "$home/calcInterfaceArea.pl $res_files | tail -n1\n";
my $area = `$home/calcInterfaceArea.pl $res_files | tail -n1`;
chomp $area;

#my $xlnum = `wc -l dist_constraints | awk '{ print \$1 }'`; chomp $xlnum;
#my $solnum = `wc best_trans0 best_trans1 best_trans2 best_trans3 best_trans4 | tail -n1 | awk '{ print \$1 }'`;
#chomp $solnum;
#my $best_xl_thr = `grep \"|\" xldocking0.res xldocking1.res xldocking2.res xldocking3.res xldocking4.res | grep -v rmsd | cut -d '|' -f13 | sort | tail -n1`;
#chomp $best_xl_thr;
#my $cdr3length = `wc -l cdrs3 | awk '{ print \$1 }'`; chomp $cdr3length;

open OUT, ">stat.csv";
#print OUT "$folder_name , $cdr3length, $solnum, $xlnum, $best_xl_thr, $precision, $soap, $area";
print OUT "$folder_name , $precision, $soap, $area";
close OUT;
`cat epi.csv >> stat.csv`;
