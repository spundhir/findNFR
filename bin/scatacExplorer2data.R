#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
	make_option(c("-g", "--genome"), default = "mm10", help="genome for which to fetch the data (default = %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$genome)) {
	cat("\nProgram: scatacExplorer2data.R (R script to fetch scATAC-seq data using scATAC.Explorer package)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: sachin.pundhir@finsenlab.dk\n");
	print_help(parser)
	q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load packages
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scATAC.Explorer))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## start analysis
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
meta = queryATAC(metadata_only = TRUE)[[1]]
acc = meta %>% dplyr::filter(Genome_Build == opt$genome) %>% dplyr::pull(Accession) %>% unique
df <- lapply(acc, function(x) {
  res = queryATAC(accession = x);
  do.call(rbind, lapply(names(res), function(y) {
    do.call(rbind, lapply(metadata(res[[y]])$regions, function(z) c(unlist(strsplit(z, "-")), y)))
    })) %>% as.data.frame %>% dplyr::mutate(accession = x)
  })

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## output file containing tissue specificity information
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
write.table(do.call(rbind, df), file="", sep="\t", col.names = T, row.names = F, quote = F)

q()
