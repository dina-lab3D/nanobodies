#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
require "$home/Util.pm";


if ($#ARGV < 0) {
  print "Usage: iterate_generateModels.pl <antigen.pdb (full path)> <dir1> <dir2> ...\n";
  exit;
}

my $currdir = cwd;

for(my $i=1; $i<$#ARGV+1; $i++) {
    my $script_name = "sc".$ARGV[$i].".sh";
    open OUT, ">$script_name";
    print OUT "#!/bin/tcsh\n";
    print OUT "cd $currdir\n";
    print OUT "$home/generateModels.pl $ARGV[$i] $ARGV[0]\n";
    close OUT;
    my $cmd = "sbatch --time=8:0:0 $script_name";
    print "$cmd\n";
    `$cmd`;
}
