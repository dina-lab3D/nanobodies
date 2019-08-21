#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use Cwd qw(cwd);

my $home = "$FindBin::Bin";
my $pd_home = "/cs/staff/dina/projects2/PatchDock/";
my $imp_home = "/cs/labs/dina/dina/libs/imp_fast/";

if ($#ARGV < 0) {
  print "Usage: dockXL.pl <antigen_pdb> <folder1> <folder2>...\n";
  exit;
}

my $antigen = $ARGV[0];
my $antigenPDB = basename($antigen);
for(my $i=1; $i<$#ARGV+1; $i++) {
  my $dirname = $ARGV[$i];
  `cp $antigen $dirname/`; # copy the antigen
  chdir $dirname;
  print $dirname;
  my $pdbs = `grep LOOP model.log | sort -nk4 | head -n5 | awk '{print \$2}' | tr '\n' ' '`;
  print $pdbs;
  my @list = split(' ', $pdbs);

  my $currdir = cwd;
  open OUT, ">dscript.sh";
  print OUT "#!/bin/tcsh\n";
  print OUT "cd $currdir\n";

  for (my $j=0; $j < 5; $j++) {
    my $loopfile = "nb_loop_" . $j . ".pdb";

    # set nanobody chain id to H
    my $cmd = "/cs/staff/dina/scripts/chainger.pl $loopfile ' ' 'H'";
    print OUT "$cmd\n";

    # parameter file
    $cmd = "$pd_home/buildParamsXlinksAA.pl $loopfile $antigenPDB 4.0 AA";
    print OUT "$cmd\n";
    $cmd = "mv params.txt params$j.txt";
    print OUT "$cmd\n";
    $cmd = "cp mycdrs3 cdrs3";
    print OUT "$cmd\n";
    $cmd = "cp myframe frame";
    print OUT "$cmd\n";

    # patchdock
    if(not (-e "xldocking$j.res" and (-s "xldocking$j.res" > 1000))) {
      $cmd = "$pd_home/patch_dock.Linux params$j.txt xldocking$j.res";
      print OUT "$cmd\n";
      $cmd = "$pd_home/PatchDockOut2Trans.pl xldocking$j.res > trans$j";
      print OUT "$cmd\n";
    }

    # soap score
    if(not (-e "xlsoap_score$j.res" and (-s "xlsoap_score$j.res" > 200))) {
      $cmd = "$imp_home/setup_environment.sh $imp_home/bin/soap_score $antigenPDB $loopfile trans$j -o xlsoap_score$j.res";
      print OUT "$cmd\n";
    }
  }

  close OUT;

  `sbatch --time=48:0:0 dscript.sh`;
  chdir "..";
}
