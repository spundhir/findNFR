#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing tf enrichment dynamics (can be stdin)"),
	make_option(c("-o", "--outFile"), help="output .rds file"),
    make_option(c("-n", "--minFreq"), default="50", help="minimum frequency of target sequences in each class (defaut=%default)"),
    make_option(c("-y", "--minOverlap"), default="100", help="minimum frequency of overlapping tf peaks (defaut=%default)"),
    make_option(c("-p", "--minOddsRatio"), default="3", help="minimum odd ratio of at least one class (defaut=%default)"),
    make_option(c("-d", "--diffFreq"), default="0", help="minimum difference in enrichment between categories (defaut=%default)"),
    make_option(c("-f", "--onlyDiffFreq"), action="store_true", help="less stringent criteria. use only difference in frequency"),
    make_option(c("-t", "--plotRange"), action="store_true", help="even less stringent criteria. plot motifs within a range of differential enrichment"),
    make_option(c("-v", "--plotSim"), action="store_true", help="instead plot motifs similar in frequency"),
    make_option(c("-r", "--rescale"), action="store_true", help="rescale color bar between min and max value"),
    make_option(c("-k", "--breaks"), default="0.3", help="breaks for rescale color bar (default=%default)"),
    make_option(c("-c", "--color"), default="RdBu", help="color palette (default=%default)"),
    make_option(c("-b", "--colorBias"), default="1", help="bias in color (default=%default)"),
    make_option(c("-l", "--mustIncludeMotif"), help="name of motifs that must be included in the final output (if multiple, separate them by a comma)"),
    make_option(c("-R", "--clusterRows"), default=F, help="cluster rows in the heatmap (default=%default)"),
    make_option(c("-C", "--clusterCols"), default=T, help="cluster columns in the heatmap (default=%default)"),
    make_option(c("-W", "--plotWidth"), default=5, help="width of the heatmap (default=%default)"),
    make_option(c("-H", "--plotHeight"), default=5, help="height of the heatmap (default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outFile)) {
	cat("\nProgram: giggleDynAna.R (R script to plot TF enrichment dynamics)\n")
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
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(pheatmap))

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

#sig_rows <- which(rowMax(as.matrix(df[["odds_ratio"]])) > 400 & 
#                    rowSds(as.matrix(df[["odds_ratio"]])) > 30 &
#                    row.names(df[[1]]) %in% df_promoter[which(df_promoter$logCPM_expr>1),]$name)

#cat(sprintf("%d out of %d motifs passed filter criteria..", length(sig_rows), nrow(df[["overlaps"]])))
#cat("\n")

saveRDS(df, file = opt$outFile)

#if(length(sig_rows)>2) {
#    pdf(opt$outPdfFile, height=opt$plotHeight, width=opt$plotWidth)
#    #breaks <- as.vector(summary(as.vector(mat[sig_rows,])))
#    #len=2
#    #breaks1 <- seq(breaks[1], breaks[2], length=len)
#    #breaks2 <- seq(breaks[2], breaks[3], length=len)
#    #breaks3 <- seq(breaks[3], breaks[4], length=len)
#    #breaks4 <- seq(breaks[4], breaks[5], length=len)
#    #breaks5 <- seq(breaks[5], breaks[6], length=len)
#    #breaks <- c(breaks1[1:length(breaks1)], breaks2[2:length(breaks2)],
#    #            breaks3[2:length(breaks3)], breaks4[2:length(breaks4)],
#    #            breaks5[2:length(breaks5)])
#    if(!is.null(opt$rescale)) {
#        breaks <- seq(min(mat[sig_rows,]), max(mat[sig_rows,]), by=as.numeric(opt$breaks))
#        #if(length(breaks)<=12) {
#        #    myCol <- rev(brewer.pal(length(breaks)-1, "RdBu"))
#        #} else {
#        #    myCol <- colorRampPalette(c("dodgerblue4","grey97", "sienna3"))(length(breaks)-1)
#        #    #myCol <- colorpanel(n=length(breaks)-1,low="blue",mid="#f7f7f7",high="red")
#        #    #myCol <- colorpanel(n=length(breaks)-1,low="#2166ac",mid="#f7f7f7",high="#b2182b")
#        #}
#        myCol <- rev(colorRampPalette(brewer.pal(11, opt$color), bias=as.numeric(opt$colorBias))(length(breaks)-1))
#        #heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(15,25), cexCol=1, cexRow=1, breaks=breaks)
#        pheatmap(mat[sig_rows,], col=myCol, border_color = NA, cluster_rows=opt$clusterRows, cluster_cols=opt$clusterCols, breaks=breaks)
#    } else {
#        myCol <- rev(colorRampPalette(brewer.pal(11, opt$color), bias=as.numeric(opt$colorBias))(256))
#        #myCol <- rev(brewer.pal(11, opt$color))
#        #heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(15,25), cexCol=1, cexRow=1)
#        pheatmap(mat[sig_rows,], col=myCol, border_color = NA, cluster_rows=opt$clusterRows, cluster_cols=opt$clusterCols)
#        #heatmap.2(mat, trace="none", col=myCol, margins=c(15,20), cexCol=1, cexRow=1)
#    }
#    dev.off()
#
#    txtFile=sprintf("%s.txt", opt$outPdfFile)
#    mat_df <- as.data.frame(mat)
#    colnames(mat_df) <- unlist(lapply(colnames(mat_df), function(x) sprintf("%s (freq_target/freq_bkg)", x)))
#    mat_df$diff <- apply(mat_df, 1, function(x) max(x)-min(x))
#    dat_df <- as.data.frame(dat)
#    colnames(dat_df) <- unlist(lapply(colnames(dat_df), function(x) sprintf("%s (per_target)", x)))
#    dat_df$class <- 0
#    if(nrow(dat_df[which(apply(dat_df, 1, function(x) max(x)>5)==T),])>0) {
#        dat_df[which(apply(dat_df, 1, function(x) max(x)>5)==T),]$class <- 1
#    }
#    if(nrow(dat_df[which(apply(dat_df, 1, function(x) max(x)>10)==T),])>0) {
#        dat_df[which(apply(dat_df, 1, function(x) max(x)>10)==T),]$class <- 2
#    }
#    if(nrow(dat_df[which(apply(dat_df, 1, function(x) max(x)>20)==T),])>0) {
#        dat_df[which(apply(dat_df, 1, function(x) max(x)>20)==T),]$class <- 3
#    }
#    mat_df <- merge(mat_df, dat_df, by=0)
#    colnames(mat_df)[1] <- "motif"
#    write.table(mat_df[order(-mat_df[,"diff"]),], txtFile, quote=F, append=F, sep="\t", col.names=T, row.names=F)
#} else if(length(sig_rows)>=1) {
#    txtFile=sprintf("%s.txt", opt$outPdfFile)
#    cat("\nNumber of significant motifs are too low for a heatmap. Instead writing results to txtFile\n\n")
#    if(length(sig_rows)==1) {
#        write(rownames(mat)[sig_rows], file="motif_dynamics.txt", sep="\t")
#        write.table(mat[sig_rows,], txtFile, quote=F, append=TRUE, sep="\t")
#    } else {
#        write.table(mat[sig_rows,], txtFile, quote=F, append=FALSE, sep="\t")
#    }
#} else {
#    cat("\nNo significant motif is found\n")
#}
