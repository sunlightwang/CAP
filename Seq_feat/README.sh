####################################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, xi dot wang at mdc hyphen berlin dot de
####################################################

# 1. run seqFeat.R with sequence features except for TOP
## for example, uORFs and uAUGs
Rscript seqFeat.R 3T3_total.5UTR_CS.uORF_uAUG.cnt uORF

# 2. run TOPseq.R with TOP sequences
Rscript TOPseq.R 3T3_total.5UTR_CS.TOP.cnt TOP

