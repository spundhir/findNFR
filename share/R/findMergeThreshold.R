#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing peak statistics from bed2merge (can be stdin)"),
	make_option(c("-m", "--max"), help="maximum difference point as the threshold (default: 75% quantile)", action="store_true")
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

outFile <- sprintf("%s.pdf", opt$inFile)

if(is.null(opt$max)) {
    z <- z[1:which(abs(diff(z))==max(abs(diff(z))))]
    cat(max(df[which(df$density<=as.vector(summary(z)[5])),]$V1))

    pdf(outFile, width=15)
    df$color <- "black"
    df[which(df$density<as.vector(summary(z)[5])),]$color <- "red"
    barplot(df$density, names.arg = df$V1, las=2, ylim=c(0,1), col=df$color, xlab="p-value (at which to merge peaks)", ylab="# reads within peaks / # mapped reads")
    invisible(dev.off())
} else {
    slope_threshold <- z[which(abs(diff(z))==max(abs(diff(z))))]
    cat(max(df[which(df$density<slope_threshold),]$V1))

    pdf(outFile, width=15)
    df$color <- "black"
    df[which(df$density<slope_threshold),]$color <- "red"
    barplot(df$density, names.arg = df$V1, las=2, ylim=c(0,1), col=df$color, xlab="p-value (at which to merge peaks)", ylab="# reads within peaks / # mapped reads")
    #plot(z)
    #abline( v = which(abs(diff(z))==max(abs(diff(z))) ) )
    #abline( h = z[which(abs(diff(z))==max(abs(diff(z))) ) ] )
    #plot(diff(z))
    invisible(dev.off())
}

q()
