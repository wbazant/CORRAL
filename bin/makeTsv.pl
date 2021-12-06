#!/usr/bin/env perl
use strict;
use warnings;
use feature 'say';
use List::MoreUtils qw/first_index/;
use File::Basename;


my ($dir, $pattern, $chosenColumn) = @ARGV;
die "Usage: $0 dir pattern <column|cpm>" unless -d $dir && $pattern;

$chosenColumn //= 'cpm';

opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
my @files = grep {$_=~m{$pattern}} readdir($dh);
closedir $dh;

die "No files in input directory" unless @files;

my %sampleLabels;

my %result;
for my $file (@files){
  open(my $fh, "<", "$dir/$file") or die "$!: $file";
  my $name = $file;
  if( -l "$dir/$file" ){
    my $f = readlink "$dir/$file";
    $name = basename $f;
  }
  my ($sample) = $name =~ m{(.*)$pattern};
  $sampleLabels{$sample}++;
  die $file unless $sample;
  my $chosenColumnIndex;
  while(<$fh>){
    chomp;
    my @line = split "\t";

    if ($. == 1){
      $chosenColumnIndex = first_index {$_ eq $chosenColumn} @line;
      die "$file Header not recognised: $_" unless $chosenColumnIndex > -1;
    } else{
      my $taxon = $line[0];
      my $chosenValue = $line[$chosenColumnIndex];
      $result{$taxon}{$sample} = $chosenValue;
    }
  }
}

my @sampleLabels = sort keys %sampleLabels;
my $n1 = scalar @files;
my $n2 = scalar @sampleLabels;
die "Number of input files ($n1) and different value column labels ($n2) doesn't match"
  unless $n1 == $n2;


say join "\t", "", @sampleLabels;
for my $key (sort keys %result){
  my %h = %{$result{$key}};
  say join "\t", $key, map {$_ // ""} @h{@sampleLabels};
}


