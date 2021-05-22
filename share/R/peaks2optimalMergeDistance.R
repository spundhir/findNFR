#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--firstFile"), help="file containing numeric distribution (should be first column)"),
    make_option(c("-j", "--secondFile"), help="file containing numeric distribution (should be first column)"),
    make_option(c("-m", "--method"), default="median", help="method to use (intersection, median, quantile; default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$firstFile) | is.null(opt$secondFile)) {
	cat("\nProgram: peaks2optimalMergeDistance.R (R script to find optimal distance to merge closely spaced peaks)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(caTools))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

## function to remove outliers from a distribution
remove_outliers <- function(x, na.rm = TRUE, ...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- 1.5 * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    y[x > (qnt[2] + H)] <- NA
    y
}

## start analysis
peaks <- read.table(opt$firstFile)
shuffle <- read.table(opt$secondFile)
peaks$class <- "peaks"
shuffle$class <- "shuffle"
df <- rbind(peaks, shuffle)
# df[which(df$V1<=0),]$V1 <- 1 ## for example, if there is one region from ChrY
# lower.limit <- min(log(df$V1))
# upper.limit <- max(log(df$V1))

# lower.limit <- as.numeric(summary(log(df$V1))[2])
# upper.limit <- as.numeric(summary(log(df$V1))[5])

lower.limit <- as.numeric(quantile(log(df$V1[!is.na(remove_outliers(df$V1))]), probs=0.05, type=8))
upper.limit <- as.numeric(quantile(log(df$V1[!is.na(remove_outliers(df$V1))]), probs=0.95, type=8))

peaks.density <- density(log(subset(df, class=="peaks")$V1), from=lower.limit, to=upper.limit, n=2^10)
shuffle.density <- density(log(subset(df, class=="shuffle")$V1), from=lower.limit, to=upper.limit, n=2^10)

density.difference <- peaks.density$y - shuffle.density$y

if(opt$method=="intersection") {
    intersection <- round(peaks.density$x[which(diff(density.difference > 0) != 0) + 1], digits=2)
    cat(round(exp(intersection)))
} else if(opt$method=="quantile") {
    cat(round(exp(round(quantile(peaks.density$x[which(density.difference>0)], probs=0.75, type=8)))))
} else {
    cat(round(exp(median(peaks.density$x[which(density.difference>0)]))))
}

# ggplot(df, aes(x=log(V1), fill=class)) + geom_density(alpha=.3)

q()
