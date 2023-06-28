#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("session"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file (can be stdin)"),
    make_option(c("-f", "--filterCol"), help="column to be filtered from the input file. If multiple, separate them by a comma")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$filterCol)) {
	cat("\nProgram: tabEdit.R (R script to filter select columns)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load data
if(identical(opt$inFile, "stdin")==T) {
    df <- read.table(file("stdin"), header=T, comment.char="")
} else {
    df <- read.table(opt$inFile, header=T, comment.char="")
}

#############################################
## filter columns
#############################################
filterCol <- unlist(lapply(unlist(strsplit(opt$filterCol, ",")), function(x) {
                    if(grepl("-", x)) {
                        range <- as.numeric(unlist(strsplit(x, "-")))
                        if(length(range)==1) {
                            seq(range[1], ncol(df), by=1)
                        } else {
                            seq(range[1], range[2], by=1)
                        }
                    } else {
                        as.numeric(x)
                    }
                }))

write.table(df[,filterCol], "", sep="\t", row.names=F, col.names=T, quote=F)
q()
