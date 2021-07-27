#!/usr/bin/env perl
use strict;
use warnings;
use feature 'say';

my ($name, $path) = @ARGV;
die "Usage: $0 name inputPath" unless $name and $path;
open(my $fh, "<", $path) or die "$!: $path";

while(<$fh>){
  if ($. == 1){
    die "Header not recognised: $_" unless $_ eq "taxon\tcoverage\tcpm\ttaxon_num_reads\ttaxon_num_markers\n";
    say join ("\t", "taxon", "$name");
  } else{
    chomp;
    my ($taxon, $coverage, $cpm, $taxonNumReads, $taxonNumMarkers) = split "\t";
    if($taxonNumReads >= 4.0 && $taxonNumMarkers >= 2){
      say join("\t", $taxon, $cpm);
    }
  }

}
