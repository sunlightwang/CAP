############################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
############################################

#!/bin/bash
if [ $# != 1 ]; then
  echo "Para error."
  exit 1;
fi
j=${1}

for i in Frac1R1 Frac2R1 Frac3R1 Frac4R1 Frac5R1 Frac6R1 Frac7R1 Frac1R2 Frac2R2 Frac3R2 Frac4R2 Frac5R2 Frac6R2 Frac7R2; do
  awk -vOFS='\t' '{split($1, a, /_/); print a[1],a[2]-1,a[2],"NA",$2,a[3]}' tag_boostrap/$i.tag.cnt.bt${j} | intersectBed -wo -s -a 3T3.peaks.gene_assign.bed -b - | awk -vOFS='\t' '{cnt[$4]+=$11} END{for(i in cnt){print i,cnt[i]}}' >  bt_cnt/$i.peaks.gene.cnt.bt${j}
done
perl genGeneCountMatrix.pl bt_cnt/ .peaks.gene.cnt.bt${j} ../sample.mapping 3T3.peaks.gene_assign.ID bt_cnt/3T3.peaks.gene_assign.cnt.matrix.bt${j}
Rscript diffTSS_TE.bootstrap.R $j

