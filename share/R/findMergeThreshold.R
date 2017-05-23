#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing peak statistics from bed2merge (can be stdin)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: findMergeThreshold.R  (R script to determine p-value threshold for merging peaks)\n")
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

df$density <- df$V2/df$V3
z <- df$density
slope_threshold <- z[which(abs(diff(z))==max(abs(diff(z))))]
cat(max(df[which(df$density<=slope_threshold),]$V1))

outFile <- sprintf("%s.pdf", opt$inFile)
pdf(outFile)
barplot(df$density, names.arg = df$V1, las=2, ylim=c(0,1))
plot(z)
abline( v = which(abs(diff(z))==max(abs(diff(z))) ) )
abline( h = z[which(abs(diff(z))==max(abs(diff(z))) ) ] )
plot(diff(z))
invisible(dev.off())

q()
