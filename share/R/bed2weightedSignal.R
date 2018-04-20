#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing score values (can be stdin) (format: tab separated <coor><segment><scores>distance>)"),
    make_option(c("-s", "--scoreCol"), help="column in input file that contains score information. If multiple, separate them by a comma"),
    make_option(c("-d", "--distanceCol"), help="column in input file that contains distance information"),
    make_option(c("-t", "--segmentCol"), help="column in input file that contains genome segment information compted by ChromHMM")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$scoreCol) | is.null(opt$distanceCol)) {
	cat("\nProgram: bed2weightedSignal.R (R script to computed weighted signal)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(plyr))

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
    if(is.na(distance)) {
        weighted_signal <- 0
    } else {
        all_signal <- as.numeric(unlist(strsplit(as.character(signal),",")))
        all_distance <- as.numeric(unlist(strsplit(as.character(distance),",")))/1000

        e <- exp(1)
        weighted_signal=0
        for(i in 1:length(all_distance)) {
            weight=(2*(e^(-mu*all_distance[i])))/(1+(e^(-mu*all_distance[i])))
            weighted_signal=weighted_signal+(weight*all_signal[i])
        }
    }
    return(weighted_signal)
}

## function to compute total signal
compute_total_signal <- function(signal) {
    if(is.na(signal)) {
        total_signal <- 0
    } else {
        total_signal <- sum(as.numeric(unlist(strsplit(as.character(signal), ","))))
    }
    return(total_signal)
}

## function to compute signal count
compute_count_signal <- function(signal) {
    if(is.na(signal)) {
        count_signal <- 0
    } else {
        count_signal <- length(unlist(strsplit(as.character(signal), ",")))
    }
    return(count_signal)
}

#opt$scoreCol <- as.numeric(opt$scoreCol)
opt$distanceCol <- as.numeric(opt$distanceCol)

#df$weighted_signal <- 0
#df.sub <- df[which(!is.na(df[,opt$scoreCol])),]

## MU is set to 2kb
#df.sub$weighted_signal <- apply(df.sub, 1, function(x) compute_weighted_signal(x[opt$scoreCol], x[opt$distanceCol], 2))
#df[which(!is.na(df[,opt$scoreCol])),]$weighted_signal <- df.sub$weighted_signal

j=ncol(df)+1
for(i in unlist(strsplit(opt$scoreCol, ","))) {
    i <- as.numeric(i)
    df[,j] <- apply(df, 1, function(x) compute_total_signal(x[i]))
    j=j+1
}

j=ncol(df)+1
for(i in unlist(strsplit(opt$scoreCol, ","))) {
    i <- as.numeric(i)
    ## MU is set to 2kb
    df[,j] <- apply(df, 1, function(x) compute_weighted_signal(x[i], x[opt$distanceCol], 2))
    j=j+1
}

j=ncol(df)+1
for(i in unlist(strsplit(opt$scoreCol, ","))) {
    i <- as.numeric(i)
    df[,j] <- apply(df, 1, function(x) compute_count_signal(x[i]))
    j=j+1
}

if(!is.null(opt$segmentCol)) {
    opt$segmentCol <- as.numeric(opt$segmentCol)
    segments <- c("enhancer_active", "enhancer_poised", "enhancer_primed", "promoter_active", 
                "promoter_bivalent", "repressed_heterochromatin", "repressed_polycomb", "transcribed")

    df_segments <- t(apply(df, 1, function(x) merge(plyr::count(segments), plyr::count(unlist(strsplit(as.character(x[11]),","))), by="x", all=T)$freq.y[1:length(segments)]))
    df <- cbind(df, df_segments)
    write.table(df[,c(1:(opt$segmentCol-1),(opt$distanceCol+1):ncol(df))], "", sep="\t", row.names=F, col.names=F, quote=F)
} else {
    write.table(df[,c(1:(as.numeric(unlist(strsplit(opt$scoreCol, ",")))[1]-1),(opt$distanceCol+1):ncol(df))], "", sep="\t", row.names=F, col.names=F, quote=F)
}

#barplot(df[order(-df$weighted_signal),]$weighted_signal)

q()
