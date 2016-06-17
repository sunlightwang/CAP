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

############################
padj_cutoff <- 0.01
fc_cutoff <- log2(1.5)
############################

pdf("alt_TSS.TE_fc.pdf")

plot(R1$mean, R2$mean, col=rgb(0.1,0.1,0.1,0.3), pch=20, xlim=c(-8,8), ylim=c(-8,8))
legend("topleft", legend=paste0("r=", signif(cor(R1$mean, R2$mean),2)), bty="n")
plot(nonbt$FC.R1, nonbt$FC.R2, col=rgb(0.1,0.1,0.1,0.3), pch=20, xlim=c(-8,8), ylim=c(-8,8))
legend("topleft", legend=paste0("r=", signif(cor(nonbt$FC.R1, nonbt$FC.R2),2)), bty="n")

plot(nonbt$FC.R1, R1$mean, col=rgb(0.1,0.1,0.1,0.3), pch=20, xlim=c(-8,8), ylim=c(-8,8))
plot(nonbt$FC.R2, R2$mean, col=rgb(0.1,0.1,0.1,0.3), pch=20, xlim=c(-8,8), ylim=c(-8,8))

plot(abs(R1$mean), R1$std, col=rgb(0.1,0.1,0.1,0.3), pch=20, xlim=c(0,8), ylim=c(0,1))
#points(abs(R1$mean[R1$padj<padj_cutoff]), R1$std[R1$padj<padj_cutoff], col="purple", pch=20)
points(abs(R1$mean[R1$padj<padj_cutoff & abs(R1$mean)>fc_cutoff]), R1$std[R1$padj<padj_cutoff & abs(R1$mean)>fc_cutoff], col="red", pch=20)
abline(v=fc_cutoff, col="cyan", lty=2)
#legend("bottomright", legend=c(length(R1$mean), sum(R1$padj<padj_cutoff), sum(R1$padj<padj_cutoff & abs(R1$mean)>fc_cutoff)), text.col=c(1,"purple",2), bty="n")
legend("bottomright", legend=c(length(R1$mean), sum(R1$padj<padj_cutoff & abs(R1$mean)>fc_cutoff)), text.col=c(1,2), bty="n")
#
plot(abs(R2$mean), R2$std, col=rgb(0.1,0.1,0.1,0.3), pch=20, xlim=c(0,8), ylim=c(0,1))
#points(abs(R2$mean[R2$padj<padj_cutoff]), R2$std[R2$padj<padj_cutoff], col="purple", pch=20)
points(abs(R2$mean[R2$padj<padj_cutoff & abs(R2$mean)>fc_cutoff]), R2$std[R2$padj<padj_cutoff & abs(R2$mean)>fc_cutoff], col="red", pch=20)
abline(v=fc_cutoff, col="cyan", lty=2)
#legend("bottomright", legend=c(length(R2$mean), sum(R2$padj<padj_cutoff), sum(R2$padj<padj_cutoff & abs(R2$mean)>fc_cutoff)), text.col=c(1,"purple",2), bty="n")
legend("bottomright", legend=c(length(R2$mean), sum(R2$padj<padj_cutoff & abs(R2$mean)>fc_cutoff)), text.col=c(1,2), bty="n")

###########
both.sig <- (R1$padj<padj_cutoff & R1$mean < -fc_cutoff & R2$padj<padj_cutoff & R2$mean < -fc_cutoff) | (R1$padj<padj_cutoff & R1$mean> fc_cutoff & R2$padj<padj_cutoff & R2$mean> fc_cutoff)
control <- abs(R1$mean) < fc_cutoff & abs(R2$mean) < fc_cutoff
both.mean <- (R1$mean + R2$mean) / 2 
both.std.max <- apply(cbind(R1$std, R2$std), 1, max)
plot(abs(both.mean), both.std.max, col=rgb(0.1,0.1,0.1,0.3), pch=20, xlim=c(0,8), ylim=c(0,1), main="both")
points(abs(both.mean[both.sig]), both.std.max[both.sig], col="red", pch=20)
abline(v=fc_cutoff, col="cyan", lty=2)
legend("bottomright", legend=c(length(both.mean), sum(both.sig)), text.col=c(1,2), bty="n")
###########
dev.off()

