############################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
############################################

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) {
  cat("Usage: Rscript <input.filename> <output.prefix>\n")
  quit("no")
}
inputfile <- args[1]
outprefix <- args[2]

keynames <- function(data) {sapply(data, function(x) unlist(strsplit(x, '[|]'))[1])}
leftnames <- function(data) {sapply(data, function(x) unlist(strsplit(x, '\\^'))[1])}
rightnames <- function(data) {sapply(data, function(x) unlist(strsplit(x, '\\^'))[2])}

seqFeat <- read.table(inputfile)

nonbt.raw <- read.table("../altTSS_divTE/3T3.diff_TSS_TE_fc.txt")
R1.raw <- read.table("../altTSS_divTE/alt_TSS.TE.fc.R1.pval")
R2.raw <- read.table("../altTSS_divTE/alt_TSS.TE.fc.R2.pval")
UTR.len <- read.table("3T3_total.5UTR_CS.length")

# filter out minor isoforms
flt <- nonbt.raw$mRNA.percent.L.R1 > 5 & nonbt.raw$mRNA.percent.L.R2 > 5 & nonbt.raw$mRNA.percent.R.R1 > 5 & nonbt.raw$mRNA.percent.R.R2 > 5 & 
        nonbt.raw$ribo.L.R1 > 0.05 & nonbt.raw$ribo.L.R2 > 0.05 & nonbt.raw$ribo.R.R1 > 0.05 & nonbt.raw$ribo.R.R2 > 0.05 &
        R1.raw$std < 1 & R2.raw$std < 1
nonbt <- nonbt.raw[flt, ] 
R1 <- R1.raw[flt, ]
R2 <- R2.raw[flt, ]

both.mean <- (R1$mean + R2$mean) / 2 
both.std.max <- apply(cbind(R1$std, R2$std), 1, max)

########################
## seqFeat enrichment 
########################

length_btw_clusters <- function( cmp_name ) {
  sapply(1:length(cmp_name), function(x) {
    y <- unlist(strsplit(cmp_name[x], "\\^"))
    diff(UTR.len[y,])
})}

seqFeat_diff_btw_clusters <- function( cmp_name ) { 
  do.call(rbind, sapply(1:length(cmp_name), function(x) {
    y <- unlist(strsplit(cmp_name[x], "\\^"))
    apply(seqFeat[y,,drop=F], 2, diff)
}, simplify=F))}
seqFeat_min_btw_clusters <- function( cmp_name ) { 
  do.call(rbind, sapply(1:length(cmp_name), function(x) {
    y <- unlist(strsplit(cmp_name[x], "\\^"))
    apply(seqFeat[y,,drop=F], 2, min)
}, simplify=F))}

UTR_len_diff <- length_btw_clusters(rownames(R1))
seqFeat.id <- rownames(seqFeat)
common.idx <- leftnames(rownames(R1)) %in% seqFeat.id & rightnames(rownames(R1)) %in% seqFeat.id & !is.na(UTR_len_diff) & abs(UTR_len_diff) < 1000

both.mean.sel <- both.mean[common.idx]
UTR_len_diff.sel <- UTR_len_diff[common.idx]
seqFeat_diff <- seqFeat_diff_btw_clusters(rownames(R1)[common.idx])
seqFeat_min <- seqFeat_min_btw_clusters(rownames(R1)[common.idx])

neg_len_diff.idx <- UTR_len_diff.sel < 0
TE_fc.adj <- both.mean.sel
TE_fc.adj[neg_len_diff.idx] <- -TE_fc.adj[neg_len_diff.idx]
log10_len_diff <- log10(abs(UTR_len_diff.sel))
seqFeat_diff.adj <- seqFeat_diff
seqFeat_diff.adj[neg_len_diff.idx] <- -seqFeat_diff.adj[neg_len_diff.idx]

match_sampling <- function(sel.idx, ctrl.idx, feat, n) { 
  feat.range <- range(feat)
  sel.density <- density(feat[sel.idx], from=feat.range[1], to=feat.range[2], n=256)
  ctrl.density <- density(feat[ctrl.idx], from=feat.range[1], to=feat.range[2], n=256)
  relative.density <- sel.density$y / (ctrl.density$y + 1e-12)
  relative.density <- relative.density / max(relative.density)
  ret <- c()
  m <- 0
  while(1) {
    m <- m + 1
    r_unif <- runif(1)
    sample.idx <- sample(ctrl.idx,1)
    r_sample <- relative.density[which.min(abs(feat[sample.idx]-sel.density$x))]
    if(r_unif <= r_sample) ret <- c(ret, sample.idx)
    if(length(unique(ret)) == n | m > 10000000) return(ret)
  }
}

pdf(paste0(outprefix,".pdf"), width=4, height=5)
op <- par(mar = c(5,5,4,3), las=2)
#####################################
### wilcox test for each seqFeat 
seqFeat.test <- do.call(rbind, sapply(1:ncol(seqFeat_diff.adj), function(x) { 
  seqFeat.idx <- which(seqFeat_diff.adj[,x] >= 1)
  if(length(seqFeat.idx) < 10) return(data.frame(case.n=length(seqFeat.idx),ctrl.n=NA, ctrl.uniq.n=NA, seqFeat.mean=NA, ctrl.mean=NA, seqFeat.median=NA, ctrl.median=NA, U.p=NA))
  ctrl.idx <- which(seqFeat_diff.adj[,x] == 0)
  ctrl.idx.match_sampling <- match_sampling(seqFeat.idx, ctrl.idx, log10_len_diff, 500)
  ctrl.n <- length(ctrl.idx.match_sampling)
  ctrl.uniq.n <- length(unique(ctrl.idx.match_sampling))
  if(ctrl.uniq.n < 10) return(data.frame(case.n=length(seqFeat.idx), ctrl.n=ctrl.n, ctrl.uniq.n=ctrl.uniq.n, seqFeat.mean=NA, ctrl.mean=NA, seqFeat.median=NA, ctrl.median=NA, U.p=NA))
  boxplot(TE_fc.adj[seqFeat.idx], TE_fc.adj[unique(ctrl.idx.match_sampling)], notch=T, ylim=c(-2.5,2.5), pch=".", 
          pars = list(boxwex = 0.8, staplewex = 0.3, outwex = 0.1, outcol=c("pink","grey70"), whiskcol=c("deeppink","grey50"), staplecol=c("deeppink","grey50")),
          names=c("presence", "controls"), col=c("deeppink","grey50"), main=colnames(seqFeat)[x], ylab="TE log2(FC)")
  U.p <- wilcox.test(TE_fc.adj[seqFeat.idx], TE_fc.adj[unique(ctrl.idx.match_sampling)])$p.val 
  seqFeat.mean <- mean(TE_fc.adj[seqFeat.idx])
  ctrl.mean <- mean(TE_fc.adj[unique(ctrl.idx.match_sampling)])
  seqFeat.median <- median(TE_fc.adj[seqFeat.idx])
  ctrl.median <- median(TE_fc.adj[unique(ctrl.idx.match_sampling)])
  abline(h=seqFeat.median, lty=2, col="deeppink")
  legend("topleft", legend=paste0("P = ",signif(U.p,2)), bty="n")
  data.frame(case.n=length(seqFeat.idx),ctrl.n=ctrl.n, ctrl.uniq.n=ctrl.uniq.n, seqFeat.mean=seqFeat.mean, ctrl.mean=ctrl.mean, seqFeat.median=seqFeat.median, ctrl.median=ctrl.median, U.p=U.p)
}, simplify=F))
dev.off()

rownames(seqFeat.test) <- colnames(seqFeat)
write.table(seqFeat.test, paste0(outprefix,".txt"), quote=F, sep="\t")

