#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--bedFile"), help="input BED file or stdin"),
    make_option(c("-l", "--list"), help="input is a list", action="store_true")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$bedFile)) {
	cat("\nProgram: bedStat.R (R script to compute statistics on BED file)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

#if(is.na(file.info("stdin")$size)) { print("Hello"); }
if(identical(opt$bedFile, "stdin")==T) {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$bedFile)
}

## check, if file has header
if(grepl("chr.*start.*end", paste(data[1,], collapse = " ")) == T) { header <- data[1,]; data <- data[-1,]; colnames(data) <- header; }

if(!is.null(opt$list)){
    cat(sprintf("Median=%0.0f\n", median(as.numeric(data[,1]))))
    cat(sprintf("Mean=%0.0f\n", mean(as.numeric(data[,1]))))
    cat(sprintf("Sum=%0.0f\n", sum(as.numeric(data[,1]))))
    cat(sprintf("Max=%0.0f\n", max(as.numeric(data[,1]))))
    cat(sprintf("Min=%0.0f\n", min(as.numeric(data[,1]))))
    cat(sprintf("Count=%0.0f\n", length(as.numeric(data[,1]))))
}else{
    cat(sprintf("Median=%0.0f\n", median(as.numeric(data[,3])-as.numeric(data[,2]))))
    cat(sprintf("Mean=%0.0f\n", mean(as.numeric(data[,3])-as.numeric(data[,2]))))
    cat(sprintf("Max=%0.0f\n", max(as.numeric(data[,3])-as.numeric(data[,2]))))
    cat(sprintf("Min=%0.0f\n", min(as.numeric(data[,3])-as.numeric(data[,2]))))
    cat(sprintf("Count=%0.0f\n", length(as.numeric(data[,3])-as.numeric(data[,2]))))
}
