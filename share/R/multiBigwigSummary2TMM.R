#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--matrixFile"), help="input file containing multiBigwigSummary signal values (deepPlot -O argument) (can be stdin)"),
    make_option(c("-l", "--normMethod"), default="TMM", help="normalization method, TMM, TMMwsp, RLE, upperquartile (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$matrixFile)) {
	cat("\nProgram: multiBigwigSummary2TMM.R (R script to compute scaling factor to normalize bigWig file)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

#if(is.na(file.info("stdin")$size)) { print("Hello"); }
if(identical(opt$matrixFile, "stdin")==T) {
    data <- read.table(file("stdin"), header=T, comment.char = "")
} else {
    data <- read.table(opt$matrixFile, header=T, comment.char = "")
}

## check, if file has header
#if(grepl("chr.*start.*end", paste(data[1,], collapse = " ")) == T) { header <- data[1,]; data <- data[-1,]; colnames(data) <- header; }

suppressPackageStartupMessages(library(edgeR))

data <- data[,c(-1,-2,-3)]
NF <- edgeR::calcNormFactors(data, method=opt$normMethod) ## calculate normalization factor
LS <- colSums(data) ## calculate library size
SF <- (LS * NF)/1000000 ## calculate scaling factor
as.data.frame(SF)

q()
