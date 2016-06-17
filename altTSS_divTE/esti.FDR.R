############################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
############################################

keynames <- function(data) {sapply(data, function(x) unlist(strsplit(x, '[|]'))[1])}

nonbt.raw <- read.table("3T3.diff_TSS_TE_fc.txt")
R1.raw <- read.table("alt_TSS.TE.fc.R1.pval")
R2.raw <- read.table("alt_TSS.TE.fc.R2.pval")

flt <- nonbt.raw$mRNA.percent.L.R1 > 1 & nonbt.raw$mRNA.percent.L.R2 > 1 & nonbt.raw$mRNA.percent.R.R1 > 1 & nonbt.raw$mRNA.percent.R.R2 > 1 & 
        nonbt.raw$ribo.L.R1 > 0 & nonbt.raw$ribo.L.R2 > 0 & nonbt.raw$ribo.R.R1 > 0 & nonbt.raw$ribo.R.R2 > 0 &
        R1.raw$std < 1 & R2.raw$std < 1

nonbt <- nonbt.raw[flt, ] 
R1 <- R1.raw[flt, ]
R2 <- R2.raw[flt, ]


##########################
#FDR
getFDR.shuffle_gene <- function(Rep1, Rep2, p.cutoff=sort(c(10^(-5:0), 5*10^(-5:-1))), fc.cutoff=c(0:10)*0.2) {
    # Rep1, Rep2 are matrices after get_pval
  Reps.ID <- intersect(row.names(Rep1), row.names(Rep2))
  Rep1.sf <- Rep1[Reps.ID, ]
  Rep2.sf <- Rep2[Reps.ID, ]
  Rep1.sf[,"mean"] <- sample(Rep1.sf[,"mean"], nrow(Rep1.sf))
  Rep1.sf[,"padj"] <- sample(Rep1.sf[,"padj"], nrow(Rep1.sf))
  Rep2.sf[,"mean"] <- sample(Rep2.sf[,"mean"], nrow(Rep2.sf))
  Rep2.sf[,"padj"] <- sample(Rep2.sf[,"padj"], nrow(Rep2.sf))
  #   row.names(Rep2.sf) <- sample(row.names(Rep2.sf), length(row.names(Rep2.sf)))
  
  FDR <- sapply(p.cutoff, function(x) { 
    sapply(fc.cutoff, function(y) {
      D1.p <- row.names(Rep1)[Rep1$padj < x & Rep1$mean > y]
      D2.p <- row.names(Rep2)[Rep2$padj < x & Rep2$mean > y]
      D1.n <- row.names(Rep1)[Rep1$padj < x & Rep1$mean < -y]
      D2.n <- row.names(Rep2)[Rep2$padj < x & Rep2$mean < -y]
      D <- length(intersect(D1.p, D2.p)) + length(intersect(D1.n, D2.n))
      
      F1.p <- row.names(Rep1.sf)[Rep1.sf$padj < x & Rep1.sf$mean > y]
      F2.p <- row.names(Rep2.sf)[Rep2.sf$padj < x & Rep2.sf$mean > y]
      F1.n <- row.names(Rep1.sf)[Rep1.sf$padj < x & Rep1.sf$mean < -y]
      F2.n <- row.names(Rep2.sf)[Rep2.sf$padj < x & Rep2.sf$mean < -y]
      FP <- length(intersect(F1.p,F2.p)) + length(intersect(F1.n, F2.n))
      
      FDR <- FP / D * 100
      #FDR <- D 
    })
  })
  colnames(FDR)<- paste0("p=",p.cutoff)
  row.names(FDR)<- paste0("FC=",2^fc.cutoff)
  return(FDR)
}

get_sig_N <- function(Rep1, Rep2, p.cutoff=sort(c(10^(-5:0), 5*10^(-5:-1))), fc.cutoff=c(0:10)*0.2) {
  Reps.ID <- intersect(row.names(Rep1), row.names(Rep2))
  Rep1.sf <- Rep1[Reps.ID, ]
  Rep2.sf <- Rep2[Reps.ID, ]
  N <- sapply(p.cutoff, function(x) {
     sapply(fc.cutoff, function(y) {
     D1.p <- row.names(Rep1)[Rep1$padj < x & Rep1$mean > y]
      D2.p <- row.names(Rep2)[Rep2$padj < x & Rep2$mean > y]
      D1.n <- row.names(Rep1)[Rep1$padj < x & Rep1$mean < -y]
      D2.n <- row.names(Rep2)[Rep2$padj < x & Rep2$mean < -y]
      D <- length(intersect(D1.p, D2.p)) + length(intersect(D1.n, D2.n))
    })
  })
  colnames(N)<- paste0("p=",p.cutoff)
  row.names(N)<- paste0("FC=",2^fc.cutoff)
  return(N)
}

get_sig_gene_N <- function(Rep1, Rep2, p.cutoff=sort(c(10^(-5:0), 5*10^(-5:-1))), fc.cutoff=c(0:10)*0.2) {
  Reps.ID <- intersect(row.names(Rep1), row.names(Rep2))
  Rep1.sf <- Rep1[Reps.ID, ]
  Rep2.sf <- Rep2[Reps.ID, ]
  N <- sapply(p.cutoff, function(x) {
     sapply(fc.cutoff, function(y) {
     D1.p <- row.names(Rep1)[Rep1$padj < x & Rep1$mean > y]
      D2.p <- row.names(Rep2)[Rep2$padj < x & Rep2$mean > y]
      D1.n <- row.names(Rep1)[Rep1$padj < x & Rep1$mean < -y]
      D2.n <- row.names(Rep2)[Rep2$padj < x & Rep2$mean < -y]
      D <- c(intersect(D1.p, D2.p), intersect(D1.n, D2.n))
      length(unique(keynames(D)))
    })
  })
  colnames(N)<- paste0("p=",p.cutoff)
  row.names(N)<- paste0("FC=",2^fc.cutoff)
  return(N)
}

Padj.cutoffs <- c(0.1, 0.05,0.01, 0.001, 0.0001) 
FC.cutoffs <- log2(seq(1, 2, 0.05))

FDR.list <- sapply(1:100, function(i) 
  getFDR.shuffle_gene(R1, R2, p.cutoff=Padj.cutoffs, fc.cutoff=FC.cutoffs), simplify=FALSE)
FDR.array <- array( do.call( "c", FDR.list ), dim=c(dim(FDR.list[[1]]),length(FDR.list)) )

FDR.mean <- apply(FDR.array, c(1,2), mean)
row.names(FDR.mean) <- row.names(FDR.list[[1]])
colnames(FDR.mean) <- colnames(FDR.list[[1]])

FDR.var <- apply(FDR.array, c(1,2), var)
row.names(FDR.var) <- row.names(FDR.list[[1]])
colnames(FDR.var) <- colnames(FDR.list[[1]])
FDR.sd <- sqrt(FDR.var)

FDR.mean
get_sig_N(R1, R2, p.cutoff=Padj.cutoffs, fc.cutoff=FC.cutoffs)
get_sig_gene_N(R1, R2, p.cutoff=Padj.cutoffs, fc.cutoff=FC.cutoffs)

library(Hmisc)
pdf("FDR.pdf", width = 6, height = 6)
#
plot(2^FC.cutoffs, FDR.mean[,"p=0.1"], type="o", col="chartreuse", lwd=2, pch=20, xlab="FC cutoff", ylab="FDR (%)", ylim=c(0,30))
errbar(2^FC.cutoffs, FDR.mean[,"p=0.1"], add=T, pch=20, FDR.mean[,"p=0.1"]+FDR.sd[,"p=0.1"], FDR.mean[,"p=0.1"]-FDR.sd[,"p=0.1"], errbar.col="chartreuse", col="chartreuse")
#
lines(2^FC.cutoffs, FDR.mean[,"p=0.05"], type="o", col="cyan3", lwd=2, pch=20, xlab="FC cutoff", ylab="FDR (%)")
errbar(2^FC.cutoffs, FDR.mean[,"p=0.05"], add=T, pch=20, FDR.mean[,"p=0.05"]+FDR.sd[,"p=0.05"], FDR.mean[,"p=0.05"]-FDR.sd[,"p=0.05"], errbar.col="cyan3", col="cyan3")
#
lines(2^FC.cutoffs, FDR.mean[,"p=0.01"], type="o", col="slateblue3", lwd=2, pch=20, xlab="FC cutoff", ylab="FDR (%)")
errbar(2^FC.cutoffs, FDR.mean[,"p=0.01"], add=T, pch=20, FDR.mean[,"p=0.01"]+FDR.sd[,"p=0.01"], FDR.mean[,"p=0.01"]-FDR.sd[,"p=0.01"], errbar.col="slateblue3", col="slateblue3")
#
lines(2^FC.cutoffs, FDR.mean[,"p=0.001"], type="o", col="darkorchid1", lwd=2, pch=20, xlab="FC cutoff", ylab="FDR (%)")
errbar(2^FC.cutoffs, FDR.mean[,"p=0.001"], add=T, pch=20, FDR.mean[,"p=0.001"]+FDR.sd[,"p=0.001"], FDR.mean[,"p=0.001"]-FDR.sd[,"p=0.001"], errbar.col="darkorchid1", col="darkorchid1")
#
lines(2^FC.cutoffs, FDR.mean[,"p=1e-04"], type="o", col="deeppink", lwd=2, pch=20, xlab="FC cutoff", ylab="FDR (%)")
errbar(2^FC.cutoffs, FDR.mean[,"p=1e-04"], add=T, pch=20, FDR.mean[,"p=1e-04"]+FDR.sd[,"p=1e-04"], FDR.mean[,"p=1e-04"]-FDR.sd[,"p=1e-04"], errbar.col="deeppink", col="deeppink")
#
abline(h=5, col="grey50", lwd=2, lty=2)
#
legend("topright", c("p=0.1","p=0.05","p=0.01","p=0.001","p=0.0001"), text.col=c("chartreuse","cyan3","slateblue3","darkorchid1","deeppink"), bty="n")
dev.off()
