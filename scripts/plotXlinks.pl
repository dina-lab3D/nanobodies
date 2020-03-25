#!/usr/bin/perl -w

use strict;
use FindBin;

my $home = "$FindBin::Bin";

if ($#ARGV != 1) {
  print "plotFit.pl <fit_file1> <dist_thr> \n";
  exit;
}

my $fitFile = $ARGV[0];
my $outFile = trimExtension($fitFile) . ".eps";
my $dist_thr = $ARGV[1];

open OUT, ">gnuplot_fit.txt";
print OUT "set terminal postscript eps size 3.0,2.0 color enhanced  linewidth 2.5 font 'Helvetica,18';  set output '$outFile';\n";
print OUT "set encoding iso_8859_1;set xlabel 'distance (\305)';set ylabel 'Frequency' offset 0.5\n";

print OUT "set style line 11 lc rgb '#808080' lt 1; set border 3 back ls 11;set xtics nomirror out scale 0.75;set ytics nomirror scale 0.75;set format y '%.0f%%';set format x '%2.f'\n";
print OUT "set style line 11 lc rgb '#808080' lt 1;set border 3 back ls 11\n";

print OUT "set yrange [0:*]\n";
print OUT "set xrange [0:*]\n";

print OUT "set linetype 1 lc rgb 'red'\n";
print OUT "set linetype 2 lc rgb 'blue'\n";

print OUT "set boxwidth 0.7;set style fill solid\n";

my $cmd = "$home/histogram.pl 3.0 5 < $ARGV[0] > hist_tmp";
print "$cmd\n";
`$cmd`;

print OUT "mycolor(x) = x > $dist_thr ? 1 : 2\n";
print OUT "plot 'hist_tmp' u 1:(\$2*100):(mycolor(\$1)) with boxes linecolor variable notitle \n";

`gnuplot gnuplot_fit.txt`;

sub trimExtension {
  my $str = shift;
  $str =~ s/\.[^.]+$//;
  return $str;
}
