############################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
############################################

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  cat("Usage: Rscript diffTSS_TE.R <cnt.matrix.filename> <out.filename>\n")
  quit("no")
}

infile <- args[1]
outfile <- args[2]

#######################################
cnt.matrix <- as.matrix(read.table(infile))
cnt.matrix.log2 <- log2(cnt.matrix+0.1)
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
colnames(mRNA.overall) <- c("R1", "R2")

ribo.overall <- do.call(rbind, sapply(1:nrow(ribo.matrix), function(x) {
  c(sum(ribo.matrix[x, 1:7]), sum(ribo.matrix[x, 8:14]))
}, simplify=F ))
colnames(ribo.overall) <- c("R1", "R2")

### filter mRNA.overall relative < 2%
mRNA.overall.percent.R1 <- do.call(rbind, tapply(1:sum(core.idx), as.factor(genes[core.idx]), function(x){
  ret <- mRNA.overall[x,1,drop=F] / sum(mRNA.overall[x,1,drop=F]) * 100
  rownames(ret) <-x 
  ret } ))
mRNA.overall.percent.R2 <- do.call(rbind, tapply(1:sum(core.idx), as.factor(genes[core.idx]), function(x){
  ret <- mRNA.overall[x,2,drop=F] / sum(mRNA.overall[x,2,drop=F]) * 100
  rownames(ret) <-x 
  ret } ))
### filter ribo.overall == 0

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
min_cnt.R1 <- sapply(1:nrow(idx), function(x) min(c(cnt.matrix[idx[x,1],1:7], cnt.matrix[idx[x,2],1:7])))
min_cnt.R2 <- sapply(1:nrow(idx), function(x) min(c(cnt.matrix[idx[x,1],7+1:7], cnt.matrix[idx[x,2],7+1:7])))
mean_cnt.R1 <- sapply(1:nrow(idx), function(x) 2 ^ mean(c(cnt.matrix.log2[idx[x,1],1:7], cnt.matrix.log2[idx[x,2],1:7])))
mean_cnt.R2 <- sapply(1:nrow(idx), function(x) 2 ^ mean(c(cnt.matrix.log2[idx[x,1],7+1:7], cnt.matrix.log2[idx[x,2],7+1:7])))

diff_TSS.TE_fc <- data.frame(FC.R1=diff_TSS.TE_fc.R1, FC.R2=diff_TSS.TE_fc.R2, min_cnt.R1=min_cnt.R1, min_cnt.R2=min_cnt.R2, mean_cnt.R1=mean_cnt.R1, mean_cnt.R2=mean_cnt.R2, 
                             mRNA.percent.L.R1=mRNA.overall.percent.R1[as.character(idx[,1]),], mRNA.percent.R.R1=mRNA.overall.percent.R1[as.character(idx[,2]),], 
                             mRNA.percent.L.R2=mRNA.overall.percent.R2[as.character(idx[,1]),], mRNA.percent.R.R2=mRNA.overall.percent.R2[as.character(idx[,2]),], 
                             mRNA.L=mRNA.overall[idx[,1],], mRNA.R=mRNA.overall[idx[,2],],
                             ribo.L=ribo.overall[idx[,1],], ribo.R=ribo.overall[idx[,2],])
rownames(diff_TSS.TE_fc) <- paste(rownames(cnt.matrix)[idx[,1]], rownames(cnt.matrix)[idx[,2]], sep="^")

write.table(diff_TSS.TE_fc, outfile, sep="\t", quote=F)


