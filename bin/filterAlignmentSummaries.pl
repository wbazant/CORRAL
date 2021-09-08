#!/usr/bin/env perl
use strict;
use warnings;
use feature 'say';

my ($name, $path) = @ARGV;
die "Usage: $0 name inputPath" unless $name and $path;
open(my $fh, "<", $path) or die "$!: $path";

while(<$fh>){
  if ($. == 1){
    chomp;
    my @header = split "\t", $_;
    die "Header not recognised: $_" unless $header[0] eq "taxon" && $header[2] eq "cpm" && $header[4] eq "taxon_num_markers";
    say join ("\t", "taxon", "$name");
  } else{
    chomp;
    my ($taxon, $_1, $cpm, $_3, $taxonNumMarkers) = split "\t";
    if($taxonNumMarkers >= 2){
      say join("\t", $taxon, $cpm);
    }
  }

}
