#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inDir"), help="input directory containing output from gene2geneDistStat")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inDir)) {
	cat("\nProgram: gene2geneDistStat.R (R script to compute genome wide distance statistics between genes))\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))

p <- lapply(list.files(opt$inDir, pattern="chr", full.names = T), function(FILE) {
    df <- read.table(pipe(sprintf("cut -f 4,12,17 %s", FILE)))
    test <- dcast(as.data.table(df, 2000), V1 ~ V2, value.var="V3")
    row.names(test) <- test$V1
    test$V1 <- NULL
    test[is.na(test)] <- 0
    # matrix2pheatmap(test, clusterRows = T, clusterCols = T, bias = 3)
    # head(test[1:10,1:10])

    #summary(log(apply(test, 1, function(x) min(x[x>0]))))
    lowLimit <- quantile(log(apply(test, 1, function(x) min(x[x>0]))))[2]
    upLimit <- quantile(log(apply(test, 1, function(x) min(x[x>0]))))[4]

    ggplot(melt(log(apply(test, 1, function(x) min(x[x>0])))), aes(value)) + geom_density(fill="#9ebcda") + theme_bw() + xlab("") + xlim(c(0,20)) +
        geom_vline(xintercept = c(lowLimit, upLimit), lty=2) +
        geom_text(aes(lowLimit, 0.1, label = as.numeric(sprintf("%0.2f", lowLimit)), hjust = 2), size = 6) +
        geom_text(aes(upLimit, 0.1, label = as.numeric(sprintf("%0.2f", upLimit)), hjust = -1), size = 6) +
        ggtitle(sprintf("%s (N=%s)", gsub("^.*\\.", "", FILE), nrow(test))) + ylab("Density")
})

pdf(sprintf("%s/gene2geneDistStat.pdf", opt$inDir), width=20, height=15)
    ggarrange(plotlist = p)
dev.off()

q()
