#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing score values (can be stdin) (format: <gene coor> score (...) distance [class]; scores can be comma separated)"),
    make_option(c("-s", "--scoreCol"), help="column in input file that contains score information. If multiple, separate them by a comma"),
    make_option(c("-d", "--distanceCol"), help="column in input file that contains distance information"),
    make_option(c("-t", "--segmentCol"), help="column in input file that contains enhancer class information (eg. computed by ChromHMM | active, primed etc)"),
    make_option(c("-w", "--weight"), default=1000, help="weight parameter, higher it is less steep would be the reduction in signal with increase in distance to TSS (default: %default)")
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

opt$distanceCol <- as.numeric(opt$distanceCol)

#############################################
## function to compute total signal
#############################################
compute_total_signal <- function(signal) {
    if(is.na(signal)) {
        total_signal <- 0
    } else {
        total_signal <- sum(as.numeric(unlist(strsplit(as.character(signal), ","))))
    }
    return(total_signal)
}

#############################################
## function to compute weighted signal
## Idea source: https://www.ncbi.nlm.nih.gov/pubmed/27626379
## e <- exp(1);
## mu <- 2;
## plot(unlist(lapply(1:50000, function(dist) { dist <- dist/(10*1000); (2*(e^(-mu*dist)))/(1+((e^(-mu*dist)))); })), xlab="dist_to_tss", ylab="weight", ylim=c(0,1))
#############################################
compute_weighted_signal <- function(signal, distance, mu) {
    # signal <- df[3,8]
    # distance <- df[3,9]
    if(is.na(distance)) {
        weighted_signal <- 0
    } else {
        all_signal <- as.numeric(unlist(strsplit(as.character(signal),",")))
        all_distance <- as.numeric(unlist(strsplit(as.character(distance),",")))/as.numeric(opt$weight)

        e <- exp(1)
        weighted_signal=0
        for(i in 1:length(all_distance)) {
            weight=(2*(e^(-mu*all_distance[i])))/(1+(e^(-mu*all_distance[i])))
            weighted_signal=weighted_signal+(weight*all_signal[i])
        }
    }
    return(weighted_signal)
}

#############################################
## compute total signal
#############################################
j=ncol(df)+1
for(i in unlist(strsplit(opt$scoreCol, ","))) {
    i <- as.numeric(i)
    df[,j] <- apply(df, 1, function(x) compute_total_signal(x[i]))
    j=j+1
}

#############################################
## compute weighted signal
#############################################
j=ncol(df)+1
for(i in unlist(strsplit(opt$scoreCol, ","))) {
    i <- as.numeric(i)
    ## MU is set to 2kb
    df[,j] <- apply(df, 1, function(x) compute_weighted_signal(x[i], x[opt$distanceCol], 2))
    j=j+1
}

#############################################
## compute enhancer count
#############################################
COUNT <- unlist(lapply(df[,opt$distanceCol], function(x) length(unlist(strsplit(as.character(x), ",")))))
df <- cbind(df, COUNT)

#############################################
## compute enhancer class information, if required
#############################################
if(!is.null(opt$segmentCol)) {
    opt$segmentCol <- as.numeric(opt$segmentCol)
    segments <- unique(unlist(lapply(df[,opt$segmentCol], function(x) unique(unlist(strsplit(as.character(x), ","))))))
    #segments <- c("enhancer_active", "enhancer_poised", "enhancer_primed", "promoter_active", 
    #            "promoter_bivalent", "repressed_heterochromatin", "repressed_polycomb", "transcribed")

    df_segments <- t(apply(df, 1, function(x) merge(plyr::count(segments), plyr::count(unlist(strsplit(as.character(x[opt$segmentCol]),","))), by="x", all=T)$freq.y[1:length(segments)]))
    colnames(df_segments) <- plyr::count(unlist(strsplit(as.character(df[,opt$segmentCol]),",")))$x

    df <- cbind(df, df_segments)
    #df <- df[,c(1:(opt$segmentCol-1),(opt$distanceCol+1):ncol(df))]
}

#############################################
## organize final output results
#############################################
#print(which(!seq(1:ncol(df)) %in% c(as.numeric(unlist(strsplit(opt$scoreCol, ","))), opt$distanceCol, opt$segmentCol)))
#df_output <- df[,c(which(!seq(1:ncol(df)) %in% c(as.numeric(unlist(strsplit(opt$scoreCol, ","))), opt$distanceCol, opt$segmentCol)))]
if(!is.null(opt$segmentCol)) {
    #print(c(1:(as.numeric(unlist(strsplit(opt$scoreCol, ",")))[1]-1),(opt$distanceCol+1):ncol(df)))
    df_output <- df[,c(1:(as.numeric(unlist(strsplit(opt$scoreCol, ",")))[1]-1),(opt$distanceCol+2):ncol(df))]
} else {
    df_output <- df[,c(1:(as.numeric(unlist(strsplit(opt$scoreCol, ",")))[1]-1),(opt$distanceCol+1):ncol(df))]
}

#############################################
## write final results
#############################################
write.table(df_output, "", sep="\t", row.names=F, col.names=T, quote=F)

#barplot(df_output[order(-df$weighted_signal),]$weighted_signal)

q()
