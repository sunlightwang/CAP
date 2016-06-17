############################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
############################################

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1) {
  cat("Usage: Rscript diffTSS_TE.bootstrap.R <bt_no>\n")
  quit("no")
}

bt_i <- args[1]
infile <- paste0("bt_cnt/3T3.peaks.gene_assign.cnt.matrix.bt", bt_i)
outfile <- paste0("bt_alt_TSS_TE_fc/3T3.alt_TSS_TE_fc.bt", bt_i)

#######################################
cnt.matrix <- read.table(infile)
norm.fac <- c(180.5, 113, 172.5, 79, 107, 90, 120.5, 248.5, 168.5, 290.5, 80, 128.5, 188, 276) # 7th col, upper-quantile
num.ribo <- c(0,0,1,2.5,4.5,7.5,12)

#######################################
gene.names <- function(data) {sapply(data, function(x) unlist(strsplit(x, '[|]'))[1])}
type.names <- function(data) {sapply(data, function(x) unlist(strsplit(x, '[|]'))[2])}

genes <- gene.names(rownames(cnt.matrix))
types <- type.names(rownames(cnt.matrix))

core.idx <- types %in% c("5UTR", "5UTR_INTRON", "5UTR_UP1k")

mRNA.matrix <- t(t(cnt.matrix[core.idx,]) / norm.fac)
ribo.matrix <- t(t(mRNA.matrix) * num.ribo)

mRNA.overall <- do.call(rbind, sapply(1:nrow(mRNA.matrix), function(x) {
  c(sum(mRNA.matrix[x, 1:7]), sum(mRNA.matrix[x, 8:14]))
}, simplify=F ))

ribo.overall <- do.call(rbind, sapply(1:nrow(ribo.matrix), function(x) {
  c(sum(ribo.matrix[x, 1:7]), sum(ribo.matrix[x, 8:14]))
}, simplify=F ))

nRibo <- ribo.overall / mRNA.overall + 1e-6
rownames(nRibo) <- rownames(cnt.matrix)[core.idx]

#### TSS cluster in the same genes pairwise-cmp idx
idx <- do.call(rbind, tapply(c(1:length(genes))[core.idx], as.factor(genes[core.idx]), function(x){
  if(length(x) > 1) { 
    ret <- t(combn(x,2))
    rownames(ret) <- paste0(genes[x[1]], "|", 1:nrow(ret))
    ret
  }
}))

diff_TSS.TE_fc.R1 <- apply(log2(cbind(nRibo[idx[,1],1],nRibo[idx[,2],1])),1,diff)
diff_TSS.TE_fc.R2 <- apply(log2(cbind(nRibo[idx[,1],2],nRibo[idx[,2],2])),1,diff)

diff_TSS.TE_fc <- data.frame(R1=diff_TSS.TE_fc.R1, R2=diff_TSS.TE_fc.R2)
rownames(diff_TSS.TE_fc) <- paste(rownames(cnt.matrix)[idx[,1]], rownames(cnt.matrix)[idx[,2]], sep="^")

#write.table(diff_TSS.TE_fc.R1, outfileR1, sep="\t", quote=F, row.names=F, col.names=F)
#write.table(diff_TSS.TE_fc.R2, outfileR2, sep="\t", quote=F, row.names=F, col.names=F)
write.table(diff_TSS.TE_fc, outfile, sep="\t", quote=F)

