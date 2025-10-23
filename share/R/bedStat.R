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
    cat(sprintf("Median=%0.4f\n", median((data[,1]))))
    cat(sprintf("Mean=%0.4f\n", mean((data[,1]))))
    cat(sprintf("Max=%0.4f\n", max((data[,1]))))
    cat(sprintf("Min=%0.4f\n", min((data[,1]))))
    cat(sprintf("Quantile (25)=%0.4f\n", (quantile(data[,1],na.rm = T,probs = c(0.25)))))
    cat(sprintf("Quantile (75)=%0.4f\n", (quantile(data[,1],na.rm = T,probs = c(0.75)))))
    cat(sprintf("Count=%0.4f\n", length((data[,1]))))
    cat(sprintf("Sum=%0.4f\n", sum((data[,1]))))
}else{
    cat(sprintf("Median=%0.4f\n", median((data[,3])-(data[,2]))))
    cat(sprintf("Mean=%0.4f\n", mean((data[,3])-(data[,2]))))
    cat(sprintf("Max=%0.4f\n", max((data[,3])-(data[,2]))))
    cat(sprintf("Min=%0.4f\n", min((data[,3])-(data[,2]))))
    cat(sprintf("Quantile (25)=%0.4f\n", (quantile(data[,3]-data[,2],na.rm = T,probs = c(0.25)))))
    cat(sprintf("Quantile (75)=%0.4f\n", (quantile(data[,3]-data[,2],na.rm = T,probs = c(0.75)))))
    cat(sprintf("Count=%0.4f\n", length((data[,3])-(data[,2]))))
}
