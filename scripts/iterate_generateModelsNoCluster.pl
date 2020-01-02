#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";
my $centroid_home = "/cs/staff/dina/utils/srcs/centroid";

if ($#ARGV < 0) {
  print "Usage: iterate_generateModels.pl <dir1> <dir2> ...\n";
  exit;
}

my $currdir = cwd;

for(my $i=0; $i<$#ARGV+1; $i++) {
    chdir $ARGV[$i];
    print "$ARGV[$i]\n";
    `rm -f soap_score*.pdb xl_soap*.pdb best_xl*.pdb`;
    if(-e "dist_constraints") {
        print "$ARGV[$i] xl\n";
        `$home/generateModelsTopN.pl 10`;
        `$home/generateModelsTopN.pl 10 best_xl`;
        my $center_file = $ARGV[$i] . "_soap_centers.pdb";
        `$centroid_home/centroid soap_score*.pdb > $center_file`;
        my $center_file2 = $ARGV[$i] . "xl_centers.pdb";
        `$centroid_home/centroid best_xl*.pdb > $center_file2`;
        `$home/generateModelsTopN.pl 30`;
        my $center_file3 = $ARGV[$i] . "_soap_centers30.pdb";
        `$centroid_home/centroid soap_score*.pdb > $center_file3`;
    }
    chdir "..";
}
