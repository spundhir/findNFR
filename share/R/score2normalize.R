#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing score values (can be stdin)"),
    make_option(c("-c", "--scoreCol"), default=5, help="column in input file that contains score information (default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: score2normalize.R (R script to normalize input scores (rank, scaled rank, percentile rank)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(session))

## load data
if(identical(opt$inFile, "stdin")==T) {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$inFile)
}

## function to compute percentile rank
perc.rank <- function(x, xo)  length(x[x <= xo])/length(x)*100

## function to scale between 0 and 1
scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

## normalize scores
opt$scoreCol <- as.numeric(opt$scoreCol)

## compute rank
set.seed(1)
data$rank <- rank(data[,c(opt$scoreCol)], ties.method = "random")

rank_col <- ncol(data)

## compute scaled rank
data$scaledRank <- scale01(data[,rank_col])

## compute percentile rank
data$percentileRank <- apply(data, 1, function(x) 100*((as.numeric(x[rank_col])-1)/(length(data[,rank_col])-1)))
#data$percentileRank <- apply(data, 1, function(x) perc.rank(data[,opt$scoreCol], x[opt$scoreCol]))

## write output
## sort data
#write.table(data[order(-data$percentileRank),], "", sep="\t", row.names = F, col.names = F, quote = F)
## unsorted data
write.table(data, "", sep="\t", row.names = F, col.names = F, quote = F)

q()
