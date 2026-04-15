#!/usr/local/bin/Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing motifs (can be stdin)"),
	make_option(c("-j", "--motifName"), help="name of motifs needs to be extracted. if multiple, separate them by a comma"),
    make_option(c("-f", "--inFormat"), default="matrix", help="format of motifs in input file (homer, jaspar, meme, transfac, uniprobe, cisbp, matrix; default=%default)"),
    make_option(c("-F", "--outFormat"), default="homer", help="format of motifs in output file (homer, jaspar, meme, transfac, uniprobe, cisbp, matrix; default=%default)"),
	make_option(c("-P", "--motifPlotFile"), help="pdf file to plot the motif(s) as sequence logo"),
	make_option(c("-W", "--motifPlotWidth"), default=6, help="plot width (default=%default)"),
	make_option(c("-H", "--motifPlotHeight"), default=10, help="plot height (default=%default)"),
	make_option(c("-M", "--fullMatch"), action="store_true", help="full match of name")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$motifName)) {
	cat("\nProgram: extractMotifs.R (R script to extract motifs of interest)\n")
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

    ## remove temporary files
    invisible(file.remove(c(tmpInFile)))
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
tmpOutFile <- tempfile(fileext = ".txt")

## extract motifs of interest
#df <- df[which(grepl("^Msx3$", sapply(df, slot, "name")))]
if(!is.null(opt$fullMatch)) { 
    df <- df[grepl(paste0(paste0("^", gsub(",", "|^", opt$motifName)), "$"), sapply(df, slot, "name"))]
} else {
    df <- df[grepl(paste0("^", gsub(",", "|^", opt$motifName)), sapply(df, slot, "name"))]
}

## plot extracted motifs as sequence logo
if(!is.null(opt$motifPlotFile)) {
    pdf(opt$motifPlotFile, width=opt$motifPlotWidth, height=opt$motifPlotHeight)
    view_motifs(df)
    dev.off()
} else if(length(df) > 0) {
    if(opt$outFormat=="homer") {
        write_homer(df, tmpOutFile, overwrite = T)
    } else if(opt$outFormat=="jaspar") {
        write_jaspar(df, tmpOutFile, overwrite = T)
    } else if(opt$outFormat=="meme") {
        write_meme(df, tmpOutFile, version = 5, overwrite = T)
    } else if(opt$outFormat=="transfac") {
        write_transfac(df, tmpOutFile, overwrite = T)
    } else if(opt$outFormat=="uniprobe") {
        write_uniprobe(df, tmpOutFile, overwrite = T)
    } else if(opt$outFormat=="cisbp") {
        write_cisbp(df, tmpOutFile, overwrite = T)
    } else {
        write_matrix(df, tmpOutFile, positions = "columns", headers = T, rownames = T)
    }

    #cat(readLines(pipe(sprintf("cat %s | awk -v OFS='\t' '{ if($0 ~ /^>/) { print $1,$2\"/\"$(gsub(\">\", \"\", $1)),$3; } else { print $0; }}'", tmpOutFile))), sep="\n")
    cat(readLines(tmpOutFile), sep="\n")

    ## remove temporary files
    invisible(file.remove(c(tmpOutFile)))
} else { cat(sprintf("\nERROR: no motif found with name, %s\n\n", opt$motifName)); }

q()
