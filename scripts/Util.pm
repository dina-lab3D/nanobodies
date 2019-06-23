#!/usr//bin/perl -w

package Util;

use strict;


sub read_fasta_sequences {
  my $input_file = shift;
  my @sequences = ();
  open(IN, $input_file) or die "Cannot open $input_file\n";
  my $curr_sequence = '';
  my $curr_header;
  while(my $line = <IN>) {
    if($line =~ /^>/) {
      $curr_header = $line;
      #print "Header = $line\n";
      if(length($curr_sequence) > 0) {
        push(@sequences, ($curr_header, $curr_sequence));
        $curr_sequence = '';
      }
    } else {
      chomp($line);
      $curr_sequence .= $line;
    }
  }
  if(length($curr_sequence) > 0) {
    push(@sequences, ($curr_header, $curr_sequence));
  }
  close(IN);
  return @sequences;
}


sub get_counter_string {
    my $seq_counter = shift;
    my $seq_counter_string;
    if($seq_counter < 10000) { $seq_counter_string .= "0" ; }
    if($seq_counter < 1000) { $seq_counter_string .= "0" ; }
    if($seq_counter < 100) { $seq_counter_string .= "0"; }
    if($seq_counter < 10) { $seq_counter_string .= "0"; }
    $seq_counter_string .= $seq_counter;
    return $seq_counter_string;
}

sub get_counter_string_limit {
  my $seq_counter = shift;
  my $limit = shift;
  my $seq_counter_string;
  while ($limit / 10 > 1) {
    $seq_counter_string .= "0" ;
    $limit = $limit / 10;
  }

  $seq_counter_string .= $seq_counter;
  return $seq_counter_string;
}

1;
