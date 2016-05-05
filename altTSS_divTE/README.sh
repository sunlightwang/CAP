####################################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, xi dot wang at mdc hyphen berlin dot de
####################################################

# 1. bootstrap tags
mkdir tag_boostrap
for s in Frac1R1 Frac2R1 Frac3R1 Frac4R1 Frac5R1 Frac6R1 Frac7R1 Frac1R2 Frac2R2 Frac3R2 Frac4R2 Frac5R2 Frac6R2 Frac7R2; do
  for i in `seq 1 1000`; do
    perl bootstrap_tags.pl data/$s.tag.cnt tag_boostrap/$s.tag.cnt.bt$i
  done
done

# 2. count tags in clusters & calculate alt-TSS TE-fc
mkdir bt_cnt
mkdir bt_alt_TSS_TE_fc
for i in `seq 1 1000`; do
  cnt.core.sh $i
done

# 3. get bootstraping mean and sd
Rscript bt.pval.R

# 4. analysis the real data
Rscript diffTSS_TE.R 3T3.peaks.gene_assign.cnt.matrix 3T3.diff_TSS_TE_fc.txt 

# 5. generate plots 
Rscript gen_plot.R

# 6. FDR estimation 
Rscript esti.FDR.R
