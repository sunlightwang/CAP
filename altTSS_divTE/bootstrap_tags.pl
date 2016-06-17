#!/usr/bin/perl -w
## This file is part of CAPTRE
########################################
#
# File Name:
#   bootstrap_tags
# 
# Description:
#   
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
# 
# Date:
#   Wed Jun  3 22:08:01 CEST 2015
#
########################################

use Math::Random::OO::Bootstrap;
use strict;
my $usage = "$0 <infile> <outfile>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

#$prng->seed(42);

my @tags;
while(<IN>) { 
  chomp;
  my @a = split; 
  push @tags, ($a[0]) x $a[1];
}

#print STDERR "Reading done!\n";

my $bt_tmp = &bootstrap(\@tags);
my @bt_tags = @{$bt_tmp};

#print STDERR "Bootstraping done!\n";

my %bt_hash; 
for(my $i=0; $i<@bt_tags; $i++) { 
  $bt_hash{$bt_tags[$i]} ++;
}

#print STDERR "Collapsing done!\n";

foreach my $k (sort {my @aa=split "_",$a; my @bb=split "_",$b; $aa[0] cmp $bb[0] || $aa[2] cmp $bb[2] || $aa[1] <=> $bb[1]} keys %bt_hash) { 
  print OUT $k,"\t",$bt_hash{$k},"\n"; 
}

close IN;
close OUT;

sub bootstrap {
  my @s = @{$_[0]};
  my $b = Math::Random::OO::Bootstrap->new(@s);
  my @ret; 
  for(my $i=0; $i<@s; $i++) {
    push @ret, $b->next(); # draws randomly from the sample
  }
  return(\@ret);
}
