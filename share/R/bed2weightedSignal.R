#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing score values (can be stdin)"),
    make_option(c("-s", "--scoreCol"), default=8, help="column in input file that contains score information (default=%default)"),
    make_option(c("-d", "--distanceCol"), default=10, help="column in input file that contains distance information (default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: bed2weightedSignal.R (R script to computed weighted signal)\n")
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
    df <- read.table(file("stdin"))
} else {
    df <- read.table(opt$inFile)
}

## function to compute weighted signal
compute_weighted_signal <- function(signal, distance, mu) {
    # signal <- df[3,8]
    # distance <- df[3,9]
    all_signal <- as.numeric(unlist(strsplit(as.character(signal),",")))
    all_distance <- as.numeric(unlist(strsplit(as.character(distance),",")))/1000

    e <- exp(1)
    weighted_signal=0
    for(i in 1:length(all_distance)) {
        weight=(2*(e^(-mu*all_distance[i])))/(1+(e^(-mu*all_distance[i])))
        weighted_signal=weighted_signal+(weight*all_signal[i])
    }
    return(weighted_signal)
}

opt$scoreCol <- as.numeric(opt$scoreCol)
opt$distanceCol <- as.numeric(opt$distanceCol)

df$weighted_signal <- 0
df.sub <- df[which(!is.na(df[,opt$scoreCol])),]

## MU is set to 2kb
df.sub$weighted_signal <- apply(df.sub, 1, function(x) compute_weighted_signal(x[opt$scoreCol], x[opt$distanceCol], 2))
df[which(!is.na(df[,opt$scoreCol])),]$weighted_signal <- df.sub$weighted_signal
write.table(df, "", sep="\t", row.names=F, col.names=F, quote=F)
#barplot(df[order(-df$weighted_signal),]$weighted_signal)

q()
