####################################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, xi dot wang at mdc hyphen berlin dot de
####################################################

# 1. Reconstruct 5'UTRs
perl fmt.TSS2TIS.pl 3T3.peaks.summit.bed mm10.refseq.TIS.bed 3T3_total.splice_site.bed 3T3_total.5UTR.bed

# 2. Select 5'UTRs with only alternative splicing or without splicing
awk '$5==100 && $10>0 && $4!~/INTRON/' 3T3_total.5UTR.bed | perl -ne '{@a=split; @b=split /\|/,$a[3]; $d=abs($a[1]-$a[2]); $key=join "_",@b[0..3]; if((exists $dist{$key} && $d < $dist{$key}) || (! exists $dist{$key})) {$dist{$key}=$d; $lines{$key}=$_;}} END {foreach $k (keys %lines) {print $lines{$k}}}' > 3T3_total.5UTR_CS.bed

