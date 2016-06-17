#!/usr/bin/perl -w
# This file is part of CAPTRE
########################################
#
# File Name:
#   genAllelicGeneReadCountMatrix.pl
# 
# Description:
#   generate allelic gene read count matrix from 
#   individual samples' read count files
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
# 
# Date:
#   Mon Aug 26 22:08:47 CEST 2013
#
########################################

use strict;
my $usage = "$0 <infile_prefix> <infile_suffix> <sample_map_file> <geneID_file> <outfile>\n";
my $infile_prefix = shift || die $usage;
my $infile_suffix = shift || die $usage;
my $sample = shift || die $usage;
my $genefile = shift || die $usage;
my $outfile = shift || die $usage;
open(GEN, $genefile) || die "cannot open $genefile for reading\n";
my @gene = <GEN>; 
my %counts;
my @samples;
open(IN, $sample) || die "Can't open $sample for reading!\n";
while(<IN>) { 
  chomp;
  my @a = split ":", $_;
  $a[1] = $a[0] if(scalar(@a) == 1);
  push @samples, $a[1];
  for(my $i=0;$i<@gene;$i++) { 
    chomp $gene[$i];
    $counts{$gene[$i]}{$a[1]} = 0;
  }
  open(CNT, $infile_prefix.$a[0].$infile_suffix) || die "cannot open $infile_prefix.$a[0].$infile_suffix for reading\n";
  while(<CNT>) { 
    chomp; 
    my @c = split;
    #$counts{$c[0]}{$a[1]} = $c[1];
    $counts{$c[0]}{$a[1]} += $c[1];
  }
}


open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

for(my $j=0;$j<@samples;$j++) {
  print OUT "$samples[$j]\t";
}
print OUT "\n";
for(my $i=0;$i<@gene;$i++) { 
  chomp $gene[$i];
  print OUT "$gene[$i]\t";
  for(my $j=0;$j<@samples;$j++) {
    print OUT "$counts{$gene[$i]}{$samples[$j]}\t";
  }
  print OUT "\n";
}
close IN;
close OUT;
