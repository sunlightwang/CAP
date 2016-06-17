############################################
# This file is part of CAPTRE
#
# Author:
#   Xi Wang, Xi.Wang@mdc-berlin.de
############################################

file.pattern <- "3T3.alt_TSS_TE_fc.bt"
file <- dir("bt_alt_TSS_TE_fc", pattern=file.pattern, full.names = TRUE)
file.cont <- lapply(file, function(x) read.table(x,stringsAsFactors = FALSE))
alt_TSS.TE.fc.R1 <- sapply(file.cont, `[[`, "R1")
alt_TSS.TE.fc.R2 <- sapply(file.cont, `[[`, "R2")
rownames(alt_TSS.TE.fc.R1) <- rownames(file.cont[[1]])
rownames(alt_TSS.TE.fc.R2) <- rownames(file.cont[[1]])

bootstrap.pval <- function(data) { 
  mean <- rowMeans(data)
  std <- sqrt(apply(data, 1, var))
  pval <- apply(cbind(mean, std), 1, function(X) 2*(1-pnorm(abs(X[1])/X[2])) )
  padj <- p.adjust(pval, method="BH")
  data.frame(mean=mean, std=std, pval=pval, padj=padj)
  #ret <- data.frame(mean=mean, std=std, pval=pval, padj=padj)
  #rownames(ret) <- rownames(data)
  #ret
}

alt_TSS.TE.fc.R1.bt <-  bootstrap.pval(alt_TSS.TE.fc.R1)
alt_TSS.TE.fc.R2.bt <-  bootstrap.pval(alt_TSS.TE.fc.R2)

write.table(alt_TSS.TE.fc.R1.bt, "alt_TSS.TE.fc.R1.pval", sep="\t", quote=F)
write.table(alt_TSS.TE.fc.R2.bt, "alt_TSS.TE.fc.R2.pval", sep="\t", quote=F)
