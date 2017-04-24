#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--firstFile"), help="file containing numeric distribution (should be first column)"),
    make_option(c("-j", "--secondFile"), help="file containing numeric distribution (should be first column)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$firstFile) | is.null(opt$secondFile)) {
	cat("\nProgram: findIntersectionPoint.R (R script to find intersection point between two numeric distributions)\n")
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

peaks <- read.table(opt$firstFile)
shuffle <- read.table(opt$secondFile)
peaks$class <- "peaks"
shuffle$class <- "shuffle"
df <- rbind(peaks, shuffle)
# lower.limit <- min(log(df$V2))
# upper.limit <- max(log(df$V2))

# lower.limit <- as.numeric(summary(log(df$V2))[2])
# upper.limit <- as.numeric(summary(log(df$V2))[5])

lower.limit <- as.numeric(round(log(quantile(df$V1, probs=0.0005, type=8)), digits=2))
upper.limit <- as.numeric(round(log(quantile(df$V1, probs=0.95, type=8)), digits=2))

peaks.density <- density(log(subset(df, class=="peaks")$V1), from=lower.limit, to=upper.limit, n=2^10)
shuffle.density <- density(log(subset(df, class=="shuffle")$V1), from=lower.limit, to=upper.limit, n=2^10)

density.difference <- peaks.density$y - shuffle.density$y
intersection <- round(peaks.density$x[which(diff(density.difference > 0) != 0) + 1], digits=2)
cat(round(exp(intersection)))
#cat(exp(intersection-1))

# summary(round(peaks.density$x[which(density.difference>0)]))

# ggplot(df, aes(x=log(V2), fill=class)) + geom_density(alpha=.3)

q()
