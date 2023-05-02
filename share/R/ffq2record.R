#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--jsonFile"), help="input JSON file or stdin")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$jsonFile)) {
	cat("\nProgram: ffq2record.R (R script to retrieve url corresponding to GEO records)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## read json file
suppressPackageStartupMessages(library("rjson"))

#if(is.na(file.info("stdin")$size)) { print("Hello"); }
if(identical(opt$jsonFile, "stdin")==T) {
    df <- fromJSON(file=file("stdin"))
} else {
    df <- fromJSON(file=opt$jsonFile)
}

## parse json file
NAME <- names(df)
lst <- lapply(names(df[[NAME]]$study$runs), function(x) {
            name <- tolower(gsub("\"", "", gsub("^.*\\s+", "", unlist(strsplit(df[[NAME]]$study$runs[[x]]$title, ";"))[2])))
            record <- t(as.data.frame(lapply(1:length(df[[NAME]]$study$runs[[x]]$files), function (y) {
                t <- (as.data.frame(c(x, name, df[[NAME]]$study$runs[[x]]$files[[y]]$url, df[[NAME]]$study$runs[[x]]$title, df[[NAME]]$study$runs[[x]]$sample$attributes$`cell type`)))
                row.names(t) <- NULL
                t
            })))
            row.names(record) <- NULL
            record
        })
lst <- lst[lapply(lst,length)>0]
lst <- do.call(rbind.data.frame, lst)
row.names(lst) <- NULL
colnames(lst) <- c("sra_id", "name", "url", "title", "cell_type")
write.table(lst, "", col.names=T, row.names=F, quote=F, sep="\t")
q()
