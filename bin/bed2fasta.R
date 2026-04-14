#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
	make_option(c("-i", "--bedFile"), help="input BED file or stdin"),
    make_option(c("-g", "--genome"), default="mm10", help="genome (mm10 or hg38; default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$bedFile)) {
	cat("\nProgram: bed2fasta.R (R script to retrieve DNA sequence for input genomic coordinates)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load packages
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library("GenomicRanges"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## initialize input files based on genome
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
#if(is.na(file.info("stdin")$size)) { print("Hello"); }
if(identical(opt$bedFile, "stdin")==T) {
    data <- read.table(file("stdin"))
} else {
    data <- read.table(opt$bedFile)
}

## check, if file has header
if(grepl("chr.*start.*end", paste(data[1,], collapse = " ")) == T) { header <- data[1,]; data <- data[-1,]; colnames(data) <- header; }

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## start analysis
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## query coordinate (format: chr start end name)
colnames(data)[1:3] <- c("chr", "start", "end")
  
query = makeGRangesFromDataFrame(data[,c("chr", "start", "end")], 
                                seqnames.field = "chr", start.field = "start", end.field = "end",
                                keep.extra.columns = T, ignore.strand = T)

if(opt$genome=="hg19") {
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    bsgenome=BSgenome.Hsapiens.UCSC.hg19
    # bm = useMart(host = "grch37.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$genome=="mm9") {
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm9))
    bsgenome=BSgenome.Mmusculus.UCSC.mm9
    # bm = useMart(host = "feb2014.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
} else if(opt$genome=="hg38") {
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    bsgenome=BSgenome.Hsapiens.UCSC.hg38
    # bm = useMart(host = "apr2020.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
} else if(opt$genome=="mm10") {
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
    bsgenome=BSgenome.Mmusculus.UCSC.mm10
    # bm = useMart(host = "apr2020.archive.ensembl.org",  biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
}

seqs <- BSgenome::getSeq(bsgenome, query)
names(seqs) <- paste0(
    seqnames(query), ":",
    start(query), "-", end(query)
)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save output files
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
#writeXStringSet(seqs, filepath = "t.fasta")
#writeLines(as.vector(rbind(names(seqs), seqs)), "")
writeLines(seqs, "")
q()
