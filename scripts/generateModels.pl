#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";


if ($#ARGV != 1 and $#ARGV != 2) {
  print "Usage: generateModels.pl <dirname> <antigen.pdb (full path)> [prefix]\n";
  exit;
}

my $dirname = $ARGV[0];
my $antigenPDB = $ARGV[1];
my $prefix = "";
if($#ARGV == 2) { $prefix = $ARGV[2]; }
my $N = 10;
chdir $dirname;# || die "NO SUCH Directory: $ARGV[0]";
opendir(DIR, ".");

while (my $dirname = readdir(DIR)) {
  print "$dirname\n";
  if(-d $dirname and $dirname =~ /seq/) {
    print "$dirname\n";
    my $soapfile4 = $prefix . "soap_score4.res";
    if((-e "$dirname/$soapfile4") and (-s "$dirname/$soapfile4" > 200)) {
        chdir $dirname;
        `cp $antigenPDB .`;
        my $soapfiles = $prefix . "soap_score0.res " . $prefix . "soap_score1.res " . $prefix . "soap_score2.res " . $prefix . "soap_score3.res " . $prefix . "soap_score4.res";
        #my $cmd = `grep "|" $soapfiles | grep -v SOAP | awk '{print \$1" "\$2" "\$4}' | sort -nk3 | head -10`;
        my $cmd = `grep "|" $soapfiles | grep -v SOAP | awk '{print \$1" "\$2" "\$4}' | sort -nk3 | head -$N`;
        print $dirname;
        `grep "|" $soapfiles | grep -v SOAP | awk '{print \$1" "\$2" "\$4}' | sort -nk3 | head -$N > top$N.txt`;
        my $soap = `awk '{sum+=\$3} END {print sum/NR}' top$N.txt`; chomp $soap;

        foreach (split(/\n/,$cmd)) {
            #print "new $_\n";
            my @tmp = split (' ', $_);
            chop $tmp[0];
            my $soap_file = $tmp[0];
            my $res_num = $tmp[1];
            my $cmdt = "/cs/staff/dina/scripts/transOutput.pl $soap_file $res_num";
            `$cmdt`;
            print "$cmdt\n";
        }
        my $res_files = `ls *soap_score*.pdb`;
        print "Res: $res_files\n";
        if(length $res_files > 0) {
            my $cmd = "$home/../epiDock/interface AB H 6.0 soap_score*.pdb";
            print "$cmd\n";
            `$cmd`;
        }

        my $cdr3len = `wc cdrs3 | cut -d ' ' -f2`; chomp $cdr3len;
        $cmd =  "echo $dirname, $cdr3len, $soap, > tmp";
        `$cmd`;
        $cmd = "cat epi.csv >> tmp";
        `$cmd`;
        `tr '\n' ' ' < tmp > stat.csv`;
        `sed -i -e '\$a\\' stat.csv`;
        `echo \n >> stat.csv`;
        #print "$dirname $header $cdr3 $cdr3len $cmd1 $cmd10 $cmd100 $cmd1000  $label\n";

        `rm -f soap_score*.pdb log `;
        `cat stat.csv >> ../../epi_stat.csv`;
        chdir "..";


    } else {
        print "Docking not done $dirname\n";
    }

  }

}
