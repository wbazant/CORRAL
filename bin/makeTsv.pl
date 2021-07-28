#!/usr/bin/env perl
use strict;
use warnings;
use feature 'say';


my ($dir, $pattern) = @ARGV;
die "Usage: $0 dir" unless -d $dir;

opendir(my $dh, $dir) || die "Can't opendir $dir: $!";
my @files = grep {$_=~m{$pattern}} readdir($dh);
closedir $dh;

die "No files in input directory" unless @files;

my %keyColumnLabels;
my %valueColumnLabels;

my %result;

for my $file (@files){
  open(my $fh, "<", $file) or die "$!: $file";
  my $h = <$fh>;
  chomp $h;
  my ($keyColumnLabel, $valueColumnLabel) = split "\t", $h;
  $keyColumnLabels{$keyColumnLabel}++;
  $valueColumnLabels{$valueColumnLabel}++;

  while(<$fh>){
    chomp;
    my ($key, $value) = split "\t";
    $result{$key}{$valueColumnLabel} = $value;
  }
}

my @valueColumnLabels = sort keys %valueColumnLabels;
my $n1 = scalar @files;
my $n2 = scalar @valueColumnLabels;
die "Number of input files ($n1) and different value column labels ($n2) doesn't match"
  unless $n1 == $n2;

my ($k, @ks) = keys %keyColumnLabels;
die "Ambiguous key labels: ".join(", ", $k, @ks) unless $k and not @ks;


say join "\t", $k, @valueColumnLabels;
for my $key (sort keys %result){
  my %h = %{$result{$key}};
  say join "\t", $key, map {$_ // ""} @h{@valueColumnLabels};
}


