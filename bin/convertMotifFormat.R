#!/localhome/bric/xfd783/miniconda3/bin/Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing motifs (can be stdin)"),
	make_option(c("-o", "--outFile"), help="output file containing motifs (optional, default=stdout)"),
    make_option(c("-f", "--inFormat"), default="matrix", help="format of motifs in input file (homer, jaspar, meme, transfac, uniprobe, cisbp, matrix; default=%default)"),
    make_option(c("-F", "--outFormat"), default="homer", help="format of motifs in output file (homer, jaspar, meme, transfac, uniprobe, cisbp, matrix; default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: convertMotifFormat.R (R script to convert format of motifs)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(universalmotif))


## read input motif file
if(identical(opt$inFile, "stdin")==T) {
    tmpInFile <- tempfile(fileext = ".txt")
    #tmpInFile <- "/localhome/bric/xfd783/software/homer/data/knownTFs/vertebrates/2025_version/t.txt"
    writeLines(readLines(file("stdin")), tmpInFile)
    opt$inFile <- tmpInFile
}

if(opt$inFormat=="homer") {
    df <- suppressWarnings(read_homer(opt$inFile))
} else if(opt$inFormat=="jaspar") {
    df <- suppressWarnings(read_jaspar(opt$inFile))
} else if(opt$inFormat=="meme") {
    df <- suppressWarnings(read_meme(opt$inFile, readsites=T, readsites.meta=T))
} else if(opt$inFormat=="transfac") {
    df <- suppressWarnings(read_transfac(opt$inFile))
} else if(opt$inFormat=="uniprobe") {
    df <- suppressWarnings(read_uniprobe(opt$inFile))
} else if(opt$inFormat=="cisbp") {
    df <- suppressWarnings(read_cisbp(opt$inFile))
} else {
    df <- suppressWarnings(read_matrix(opt$inFile, positions = "columns", headers = T, rownames = T))
}

## convert motifs into desired format and write to ouput file
if(is.null(opt$outFile)) {
    tmpOutFile <- tempfile(fileext = ".txt")
    opt$outFile <- tmpOutFile
} else {
    tmpOutFile <- NULL
}

if(opt$outFormat=="homer") {
    write_homer(df, opt$outFile, overwrite = T)
} else if(opt$outFormat=="jaspar") {
    write_jaspar(df, opt$outFile, overwrite = T)
} else if(opt$outFormat=="meme") {
    write_meme(df, opt$outFile, version = 5, overwrite = T)
} else if(opt$outFormat=="transfac") {
    write_transfac(df, opt$outFile, overwrite = T)
} else if(opt$outFormat=="uniprobe") {
    write_uniprobe(df, opt$outFile, overwrite = T)
} else if(opt$outFormat=="cisbp") {
    write_cisbp(df, opt$outFile, overwrite = T)
} else {
    write_matrix(df, opt$outFile, positions = "columns", headers = T, rownames = T)
}

if(!is.null(tmpOutFile)) {
    cat(readLines(tmpOutFile), sep="\n")
}
q()
