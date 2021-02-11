#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(seqPattern))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--bedFile"), help="input BED file, can be stdin\n\t\tformat: chr strand end name score strand promoter_class"),
    #make_option(c("-m", "--motifFile"), help="input PWM file"),
    make_option(c("-g", "--genome"), default="mm10", help="genome (default: %default)"),
    make_option(c("-o", "--outPDF"), help="output pdf file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$bedFile) | is.null(opt$outPDF)) {
	cat("\nProgram: bed2seqPattern.R (R script to compute sequence motif distribution across input genomic regions)\n")
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

## load genome
if(opt$genome=="mm9") {
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm9))
    GENOME=seqlengths(Mmusculus)
} else if(opt$genome=="mm10") {
    suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))
    GENOME=seqlengths(Mmusculus)
} else if(opt$genome=="hg19") {
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    GENOME=seqlengths(Hsapiens)
} else if(opt$genome=="hg38") {
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    GENOME=seqlengths(Hsapiens)
} else {
	print_help(parser)
	q()
}

## load genomic coordinates
coor <- GRanges(seqnames = data[,1],
                ranges = IRanges(start = data[,2],
                                 end = data[,2]),
                strand = data[,6],
                promoter_class = data[,7],
                seqlengths = GENOME)

## define windows
if(opt$genome=="mm9" | opt$genome=="mm10") {
    coorFlank_oligo <- promoters(coor, upstream = 400, downstream = 800)
    coorFlankSeq_oligo <- getSeq(Mmusculus, coorFlank_oligo)
    coorFlank_pwm <- promoters(coor, upstream = 200, downstream = 200)
    coorFlankSeq_pwm <- getSeq(Mmusculus, coorFlank_pwm)
} else if(opt$genome=="hg19" | opt$genome=="hg38") {
    coorFlank_oligo <- promoters(coor, upstream = 400, downstream = 800)
    coorFlankSeq_oligo <- getSeq(Hsapiens, coorFlank_oligo)
    coorFlank_pwm <- promoters(coor, upstream = 200, downstream = 200)
    coorFlankSeq_pwm <- getSeq(Hsapiens, coorFlank_pwm)
}

## identify row index of narrow, medium and broad promoters
sIdx <- coor$promoter_class=="1_narrow"
mIdx <- coor$promoter_class=="2_medium"
bIdx <- coor$promoter_class=="3_broad"

## create output pdf file
pdf(opt$outPDF)
    ## plot WW and SS frequency plots
    par(mfrow = c(1,3), mar = c(4.5,4,1,1))
    plotPatternOccurrenceAverage(regionsSeq = coorFlankSeq_oligo[sIdx],
        patterns = c("WW", "SS"), flankUp = 400, flankDown = 800,
        smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)
    plotPatternOccurrenceAverage(regionsSeq = coorFlankSeq_oligo[mIdx],
        patterns = c("WW", "SS"), flankUp = 400, flankDown = 800,
        smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)
    plotPatternOccurrenceAverage(regionsSeq = coorFlankSeq_oligo[bIdx],
        patterns = c("WW", "SS"), flankUp = 400, flankDown = 800,
        smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)

    ## plot PWM frequency plots
    NAME="BREd,BREu,CCAAT,DCE_S_I,DCE_S_II,DCE_S_III,DPE,GC,INR,MED,MTE,TATA,XCPE1";
    lapply(unlist(strsplit(NAME, ",")), function(x) {
        x <- gsub("-.*", "", x)
        PWM <- as.matrix(read.table(pipe(sprintf("less /localhome/bric/xfd783/software/homer/data/knownTFs/vertebrates/2021_version/jaspar/JASPAR2020_POLII_pfms_jaspar.txt | sed -E 's/^.*\\[\\s+//g' | sed -E 's/\\]//g' | grep -w %s -A 4 | grep -v '>'", as.character(x)))))
        row.names(PWM) <- c("A", "C", "G", "T")

        ## plot motif frequency plots
        par(mfrow=c(1,1))
        plotMotifOccurrenceAverage(regionsSeq = coorFlankSeq_pwm[sIdx],
            motifPWM = PWM, minScore = "90%", flankUp = 200, flankDown = 200,
            smoothingWindow = 3, color = c("red3"), cex.axis = 0.9, yLabel = sprintf("Relative frequency (%s)", x))
        plotMotifOccurrenceAverage(regionsSeq = coorFlankSeq_pwm[mIdx],
            motifPWM = PWM, minScore = "90%", flankUp = 200, flankDown = 200,
            smoothingWindow = 3, color = c("green3"), cex.axis = 0.9, add=T)
        plotMotifOccurrenceAverage(regionsSeq = coorFlankSeq_pwm[bIdx],
            motifPWM = PWM, minScore = "90%", flankUp = 200, flankDown = 200,
            smoothingWindow = 3, color = c("blue3"), add = T)
        legend("topright", legend = c("sharp", "medium", "broad"), col = c("red3", "green3", "blue3"), bty = "n", lwd = 1)
    })
dev.off()

q()
