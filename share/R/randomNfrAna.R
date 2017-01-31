suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library(gridExtra))

## load arguments
ARGV <- commandArgs(TRUE)

nfr <- read.table(ARGV[1])
random <- read.table(ARGV[2])

for (i in 1:length(nfr$V5)){ nfr$V16[i] <- length(which(random$V4>nfr$V5[i]))/length(random$V4); }
#nfr$V16 <- unlist(lapply(as.list(nfr$V5), function(x) t<-length(which(random$V4>x))/length(random$V4)))
nfr$V17 <- p.adjust(nfr$V16, method="BH")
nfr$V18 <- abs((nfr$V9-nfr$V11)/(nfr$V9+nfr$V11))
nfr$V16 <- sprintf("%0.4f", nfr$V16)
nfr$V17 <- sprintf("%0.4f", nfr$V17)
nfr$V18 <- sprintf("%0.4f", nfr$V18)

write.table(nfr, file=ARGV[3], sep="\t", row.names=F, col.names=F, quote=F)
q()
