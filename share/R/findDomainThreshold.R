#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing histone signal in gap regions (can be stdin)"),
	make_option(c("-q", "--quantile"), help="75% quantile point as the threshold (default: median)", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: findDomainThreshold.R  (R script to determine p-value threshold for finding histone domains)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(session))

if(identical(opt$inFile, "stdin")==T) {
    df <- read.table(file("stdin"))
} else {
    df <- read.table(opt$inFile)
}

if(is.null(opt$quantile)) {
    cat(as.vector(summary(df[which(df$V1>as.vector(exp(summary(log(df$V1)))[2])),]$V6)[3]))
} else {
    cat(as.vector(summary(df[which(df$V1>as.vector(exp(summary(log(df$V1)))[2])),]$V6)[5]))
}

q()
