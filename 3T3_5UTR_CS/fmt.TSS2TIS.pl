#!/usr/bin/perl -w
## This file is part of CAPTRE
########################################
#
# File Name:
#   fmt.TSS2TIS.pl
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
#   Mon Jul 13 10:42:32 CEST 2015
#
########################################

use strict;
use List::Util qw(min max sum);

my $min_cnt = 8; # min summing cnt per splice site
my $next_min_cnt = 20; # next to min summing cnt per splice site
my $min_major_frac = 0.90; # min splice-event fraction per splice site for the major event
my $min_major_frac_below_min_cnt = 0.80; # min splice-event fraction per splice site for the major event if cnt < $next_min_cnt
my $max_discard_frac = 0.05; # a splice event up to this fraction will be discarded

my $usage = "$0 <tss.summit> <tis.in> <splice.in> <outfile.bed>\n";
my $tss_file = shift || die $usage;
my $tis_file = shift || die $usage;
my $splice_file = shift || die $usage;
my $outfile = shift || die $usage;
my (%tss, %tis, %ss_L, %ss_R); 

open(IN, $tss_file) || die "Can't open $tss_file for reading!\n";
while(<IN>) { 
  chomp; 
  my @a = split; 
  next unless $a[3] =~ /\|5UTR/;
  my @b = split /\|/, $a[3];
  my $chr = $a[0];
  my $strand = $a[5]; 
  #my $pos = int( ($a[1] + $a[2]) / 2);
  my $pos = ($strand eq "+")?$a[1]:$a[2];
  my $sub_key = join "_", ($chr,$pos,$strand);
  $tss{$b[0]}{$sub_key} = $a[3];
}
close IN;

open(IN, $tis_file) || die "Can't open $tis_file for reading!\n";
while(<IN>) { 
  chomp; 
  my @a = split; 
  my @b = split /\|/, $a[3]; 
  my $sub_key = join "_", ($a[0],$a[1],$a[2],$a[5]); 
  $tis{$b[0]}{$sub_key} ++ ;
}
close IN; 

open(IN, $splice_file) || die "Can't open $splice_file for reading!\n";
while(<IN>) { 
  chomp;
  my @a = split;
  my @b = split "_", $a[3];
  if($b[4] eq "L") { 
    $ss_L{$a[3]} = $a[4]; 
  } else { 
    $ss_R{$a[3]} = $a[4]; 
  }
}
close IN; 

my ($ss_L_key, $ss_R_key, $ss_L_value, $ss_R_value);
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";
foreach my $tss_gene (keys %tss) {
  my @tis = keys %{$tis{$tss_gene}};
  foreach my $tss_pos (keys %{$tss{$tss_gene}}) { 
    my @a = split "_", $tss_pos; 
    for(my $i=0; $i<@tis; $i++) { 
      my (@ss_L, @ss_L_value, @ss_R, @ss_R_value); 
      my @b = split  "_", $tis[$i];
      #print STDERR "TSS and TIS chr or strand doesn't match ".$tss_gene."\n" unless ($a[0] eq $b[0] && $a[2] eq $b[3]);
      next unless ($a[0] eq $b[0] && $a[2] eq $b[3]);
      my $id = $tss{$tss_gene}{$tss_pos}."|".$tis[$i];
      my $chr = $a[0];
      my $strand = $a[2];
      my $start = ($strand eq "+")?$a[1]:$b[2];
      my $end = ($strand eq "+")?$b[1]:$a[1];
      my $score = 100; # Score == 100 indicates CS; others indicates AS *** for AS, only the major AS output with the percentage in Score
      my $no_path = 0; # no_path == 1 if the CS 5UTR or AS major 5UTR cannot be reconstructed, ie TIS not in this isoform 
      my $nblock = 0; # block number 
      my ($starts, $lens) = (".",".");
      my $start_p = $a[1];
      my (@s, @e); 
      if($a[2] eq "+") { 
        $no_path = 1 unless($a[1] < $b[1]);
        for(my $p=$a[1]; $p<$b[1]; $p++) { 
          ###
          $ss_L_key = join "_", ($a[0], $p, $p+2, $a[2], "L");
          if (exists $ss_L{$ss_L_key}) { 
            $ss_L_value = $ss_L{$ss_L_key};
            my @d = split ";",$ss_L_value;
            my (@d_pos, @d_cnt);
            for(my $i=0; $i<@d; $i++) {
              my @tmp = split ":", $d[$i];
              $d_pos[$i] = $tmp[0];
              $d_cnt[$i] = $tmp[1];
            }
            next if(sum(@d_cnt) < $min_cnt);
            my @d_frac = &percentage(\@d_cnt); # may need to check the percentage at the acceptor site
            my ($max_v, $max_i) = &max_pos(\@d_frac);
            $score *= $max_v unless($max_v >= $min_major_frac || ($max_v >= $min_major_frac_below_min_cnt && sum(@d_cnt) < $next_min_cnt));
            next if($d_pos[$max_i] == 0); 
            if(($d_pos[$max_i]-1) > $b[1]) {
              $no_path = 1;
              last;
            }
            if(($d_pos[$max_i]-1) == $b[1]) { # start codon is at junction boundary 
              $end = $p + 1; 
              $b[1] = $p + 1;
              last;
            } else {
              push @s, $start_p;
              push @e, $p + 1; 
              $start_p = $d_pos[$max_i] - 1; 
              $p = $start_p - 1; 
            }
          }

          #$ss_R_key = join "_", ($a[0], $p, $p+2, $a[2], "R");
          #push @ss_L, $ss_L_key if exists $ss_L{$ss_L_key};
          #push @ss_L_value, $ss_L{$ss_L_key} if exists $ss_L{$ss_L_key};
          #push @ss_R, $ss_R_key if exists $ss_R{$ss_R_key};
          #push @ss_R_value, $ss_R{$ss_R_key} if exists $ss_R{$ss_R_key};
        }
        push @s, $start_p;
        push @e, $b[1];
        unless($no_path) { 
          $nblock = scalar @s; 
          for(my $ii=0; $ii<$nblock; $ii++) { 
            $e[$ii] -= $s[$ii];
            $s[$ii] -= $start;
          }
          $starts = join ",", @s;
          $lens = join ",", @e; 
        }
      } else { 
        $no_path = 1 unless($a[1] > $b[2]);
        for(my $p=$a[1]; $p>$b[2]; $p--) { 
          ###
          $ss_R_key = join "_", ($a[0], $p-2, $p, $a[2], "R");
          if (exists $ss_R{$ss_R_key}) { 
            $ss_R_value = $ss_R{$ss_R_key};
            my @d = split ";",$ss_R_value;
            my (@d_pos, @d_cnt);
            for(my $i=0; $i<@d; $i++) {
              my @tmp = split ":", $d[$i];
              $d_pos[$i] = $tmp[0];
              $d_cnt[$i] = $tmp[1];
            }
            next if(sum(@d_cnt) < $min_cnt);
            my @d_frac = &percentage(\@d_cnt); # may need to check the percentage at the acceptor site
            my ($max_v, $max_i) = &max_pos(\@d_frac);
            $score *= $max_v unless($max_v >= $min_major_frac || ($max_v >= $min_major_frac_below_min_cnt && sum(@d_cnt) < $next_min_cnt));
            next if($d_pos[$max_i] == 0); 
            if(($d_pos[$max_i]+1) < $b[2]) {
              $no_path = 1;
              last;
            }
            if(($d_pos[$max_i]+1) == $b[2]) { # start codon is at junction boundary 
              $start = $p - 1; 
              $b[2] = $p - 1;
              last;
            } else {
              push @e, $start_p;
              push @s, $p - 1; 
              $start_p = $d_pos[$max_i] + 1; 
              $p = $start_p + 1; 
            }
          }

          #$ss_L_key = join "_", ($a[0], $p, $p+2, $a[2], "L");
          #push @ss_L, $ss_L_key if exists $ss_L{$ss_L_key};
          #push @ss_L_value, $ss_L{$ss_L_key} if exists $ss_L{$ss_L_key};
          #push @ss_R, $ss_R_key if exists $ss_R{$ss_R_key};
          #push @ss_R_value, $ss_R{$ss_R_key} if exists $ss_R{$ss_R_key};
        }
        push @e, $start_p;
        push @s, $b[2];
        unless($no_path) { 
          $nblock = scalar @s; 
          for(my $ii=0; $ii<$nblock; $ii++) { 
            $e[$ii] -= $s[$ii];
            $s[$ii] -= $start;
          }
          $starts = join ",", reverse @s;
          $lens = join ",", reverse @e; 
        }
      }
      #print $tss_gene."\t".$tss_pos."\t$tis[$i]\n";
      #print "@ss_L\n@ss_L_value\n";
      #print "@ss_R\n@ss_R_value\n";
      #print "\n";
      print OUT join "\t", ($chr, $start, $end, $id, $score, $strand, $start, $start, 0, $nblock, $lens, $starts), "\n"; 
    }
  }
}
close OUT;

sub percentage {
  my @a = @{$_[0]};
  my @ret;
  my $s = sum(@a);
  for(my $i=0; $i<@a; $i++) {
    $ret[$i] = $a[$i]/$s;
  }
  return(@ret);
}

sub max_pos {
  my @a = @{$_[0]};
  my ($m, $m_i) = 0;
  for(my $i=0; $i<@a; $i++) {
    if($m < $a[$i]) {
      $m = $a[$i];
      $m_i = $i;
    }
  }
  return(($m, $m_i));
}

