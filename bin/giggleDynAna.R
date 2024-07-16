#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inFile"), help="input file containing tf enrichment dynamics (can be stdin)"),
  make_option(c("-o", "--outPdfFile"), help="output pdf image file"),
  make_option(c("-c", "--topN"), default=10, help="number of top hits to analyze from each sample (defaut=%default)"),
  make_option(c("-x", "--pVal"), default=1e-15, help="p-value cutoff (defaut=%default)"),
  make_option(c("-y", "--minOverlap"), default=50, help="minimum frequency of overlaps in each class (defaut=%default)"),
  make_option(c("-z", "--minOddsRatio"), default=8, help="minimum odds ratio (log) of at least one class (defaut=%default)"),
  make_option(c("-b", "--colorBias"), default="1", help="bias in color (default=%default)"),
  make_option(c("-l", "--mustInclude"), help="name of element that must be included in the final output (if multiple, separate them by a comma)"),
  make_option(c("-R", "--clusterRows"), default=F, help="cluster rows in the heatmap (default=%default)"),
  make_option(c("-C", "--clusterCols"), default=T, help="cluster columns in the heatmap (default=%default)"),
  make_option(c("-W", "--plotWidth"), default=5, help="width of the heatmap (default=%default)"),
  make_option(c("-H", "--plotHeight"), default=5, help="height of the heatmap (default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outPdfFile)) {
	cat("\nProgram: giggleDynAna.R (R script to plot giggle enrichment dynamics)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))

## header of DENOVO- or KNOWN-MOTIF_ENRICHMENT_DYNAMICS.TXT
# 1.  SAMPLE_ID 
# 2.  TF_NAME 
# 3.  FILE_SIZE
# 4.  OVERLAPS
# 5.  ODDS_RATIO
# 6.  FISHERS_TWO_TAIL
# 7.  FISHERS_LEFT_TAIL
# 8.  FISHERS_RIGHT_TAIL
# 9.  COMBO_SCORE
# 10. TOTAL_REGIONS

if(identical(opt$inFile, "stdin")==T) { 
    data <- read.table(file("stdin"), header=T)
} else {
    data <- read.table(opt$inFile, header=T)
}

no_rows=nrow(data)/length(unique(data$id))
data$file <- gsub(".bed.gz", "", gsub("^.*/", "", data$file))

df <- lapply(c("overlaps", "odds_ratio", "combo_score", "fishers_two_tail"), function(x) {
            t <- as.data.frame(matrix(data[,x], nrow=no_rows))
            colnames(t) <- unique(data$id)
            row.names(t) <- unique(data$file)
            t
        })
names(df) <- c("overlaps", "odds_ratio", "combo_score", "fishers_two_tail")

## identify significant overlaps
sig_rows <- which(rowMin(as.matrix(df[["fishers_two_tail"]])) < opt$pVal &
                    rowMin(as.matrix(df[["overlaps"]])) > opt$minOverlap &
                    log(rowMax(round(as.matrix(df[["odds_ratio"]])))) >= opt$minOddsRatio &
                    log(rowSds(as.matrix(df[["odds_ratio"]]))) > 4)

if(!is.null(opt$mustInclude)) {
    sig_rows = c(sig_rows, which(row.names(df$overlaps) %in% unlist(strsplit(opt$mustInclude, ","))))
}

## filter top N overlaps ordered based on combo score
sig_rows <- sig_rows[sig_rows %in% which(row.names(df$combo_score) %in% unique(unlist(lapply(colnames(df$combo_score), function(x) head(row.names(df$combo_score[order(-df$combo_score[,x]),]), opt$topN)))))]

if(length(sig_rows)>2) {
    df_sig <- cbind(df[["overlaps"]][sig_rows,] %>% mutate(name = row.names(df[["overlaps"]][sig_rows,])) %>% reshape2::melt(),
                    df[["odds_ratio"]][sig_rows,] %>% mutate(name = row.names(df[["odds_ratio"]][sig_rows,])) %>% reshape2::melt() %>% pull(value),
                    df[["fishers_two_tail"]][sig_rows,] %>% mutate(name = row.names(df[["fishers_two_tail"]][sig_rows,])) %>% reshape2::melt() %>% pull(value),
                    df[["combo_score"]][sig_rows,] %>% mutate(name = row.names(df[["combo_score"]][sig_rows,])) %>% reshape2::melt() %>% pull(value))

    colnames(df_sig) <- c("name", "class", "overlaps", "odds_ratio", "pvalue", "combo_score")
    df_sig$odds_ratio <- log(df_sig$odds_ratio)
    df_sig$pvalue <- -log10(df_sig$pvalue)

    pdf(opt$outPdfFile, height=opt$plotHeight, width=opt$plotWidth)
        ggplot(df_sig, aes(x = class, y = name)) +
            geom_point(aes(colour=odds_ratio,size=combo_score)) +
            #geom_point(aes(colour=pvalue,size=GeneDensity*100)) +
            #scale_colour_gradient(high="brown", low="yellow", name="p-value") +
            #scale_colour_gradient(low="#00441b", high="#ccece6", name="p-value") +
            scale_colour_gradient(high="#00441b", low="#99d8c9", name="log(odds_ratio") +
            #breaks=c(1e-20, 1e-15, 1e-10, 1e-5, 0.01, 0.05)) +
            #scale_size(range=c(1,10)) +
            #scale_size(range=c((min(data_sig$GeneDensity)*100),(max(data_sig$GeneDensity*100))), breaks=waiver(), labels=waiver()) +
            #scale_size_area(breaks=seq(0,100,by=10), max_size=20) +
            theme(text=element_text(size=10), axis.text.x=element_text(angle=90)) +
            theme_bw(base_size=15)

        #ggarrange(matrix2Heatmap((t(df[["odds_ratio"]][sig_rows,])), scale="col", clusterRows = F, clusterCols = T))
    dev.off()
    df[["sig"]] <- df_sig
} else {
    cat("\nNo enrichment is found\n")
    df[["sig"]] <- NA
}

#saveRDS(df, file = opt$outFile)
