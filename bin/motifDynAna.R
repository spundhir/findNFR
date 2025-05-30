#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing motif enrichment dynamics (can be stdin)"),
    make_option(c("-n", "--minFreq"), default="50", help="minimum frequency of target sequences in each class (defaut=%default)"),
    make_option(c("-p", "--minPer"), default="3", help="minimum percentage of target sequences that should have the motif (defaut=%default)"),
    make_option(c("-y", "--minIdentity"), default="0.8", help="minimum identity of denovo motif with the known motif (defaut=%default)"),
    make_option(c("-d", "--diffFreq"), default="0", help="minimum difference in enrichment between categories (defaut=%default)"),
	make_option(c("-o", "--outPdfFile"), help="output pdf image file"),
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
    make_option(c("-M", "--clusterColsMeasure"), default="euclidean", help="cluster columns in the heatmap using (default=%default)"),
    make_option(c("-W", "--plotWidth"), default=5, help="width of the heatmap (default=%default)"),
    make_option(c("-H", "--plotHeight"), default=5, help="height of the heatmap (default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile) | is.null(opt$outPdfFile)) {
	cat("\nProgram: motifDynAna.R (R script to plot motif dynamics)\n")
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
# 1. SAMPLE_ID 
# 2. MOTIF_NAME 
# 3. MOTIF
# 4. P-value
# 5. log_P-value
# 6. q-value 
# 7. No_of_target_sequences_with_motif 
# 8. %_of_target_sequences_with_motif
# 9. No_of_background_sequences_with_motif
# 10. %_of_background_sequences_with_motif
# 11. No_of_target_sequences
# 12. No_of_background_sequences

if(identical(opt$inFile, "stdin")==T) { 
    dataRaw <- read.table(file("stdin"))
} else {
    dataRaw <- read.table(opt$inFile)
}

dataRaw$V13 <- log2((dataRaw$V8+0.01)/(dataRaw$V10+0.01))
dataRaw$V14 <- gsub("^.*BestGuess:", "", dataRaw$V2)
dataRaw$V14 <- gsub("\\(.*", "", dataRaw$V14)
dataRaw$V15 <- gsub("^.*BestGuess:", "", dataRaw$V2)
dataRaw$V15 <- as.numeric(gsub("\\)", "", gsub("^.*\\(", "", dataRaw$V15)))

if(length(which( !is.na(dataRaw$V15), arr.ind=TRUE))>0) {
    data <- dataRaw[which(dataRaw$V11>=as.numeric(opt$minFreq) & dataRaw$V15>=as.numeric(opt$minIdentity) &
                           dataRaw$V2 %in% unique(dataRaw[which(dataRaw$V8>=as.numeric(opt$minPer)),]$V2)),]
} else {
    data <- dataRaw[which(dataRaw$V11>=as.numeric(opt$minFreq) &
                          dataRaw$V2 %in% unique(dataRaw[which(dataRaw$V8>=as.numeric(opt$minPer)),]$V2)),]
    data$V15 <- 0
}

if(!is.null(opt$mustIncludeMotif)) {
    data <- rbind(data, dataRaw[grep(paste(unlist(strsplit(as.character(opt$mustInclude),",")),collapse="\\/|"), dataRaw$V2),])
    data <- data[order(data$V1, data$V2),]
    #data[is.na(data$V15),]$V15 <- 0
    data <- unique(data)
}

cat(sprintf("%d out of %d motifs passed identity (%s) and percentage (%s) score criteria..", nrow(data)/length(unique(data$V1)), nrow(dataRaw)/length(unique(dataRaw$V1)), opt$minIdentity, opt$minPer))
cat("\n")

no_rows=nrow(data)/length(unique(data$V1))
tf_info <- as.data.frame(data[1:no_rows, c(14,15)])
dat <- matrix(data$V8, nrow=no_rows)
pat <- matrix(data$V4, nrow=no_rows)
## use log fold change
mat <- matrix(data$V13, nrow=no_rows)
## use q-value
#mat <- matrix(data$V6, nrow=no_rows)
data$V1 <- sprintf("%s_%s", data$V1, data$V11)
colnames(mat) <- as.vector(unique(data$V1))
row.names(mat) <- data[1:no_rows,2]
colnames(dat) <- as.vector(unique(data$V1))
row.names(dat) <- data[1:no_rows,2]
if(!is.null(opt$plotSim)) {
    #sig_rows <- which(apply(mat, 1, function(x) max(x) > 0.1 & min(x) > 0.1 & max(x)-min(x) <= as.numeric(opt$diffFreq)))
    #sig_rows <- which(apply(mat, 1, function(x) max(x)-min(x) <= as.numeric(opt$diffFreq)))
    sig_rows <- row(mat)[,1]
} else if(!is.null(opt$plotRange)) {
    sig_rows <- which(apply(mat, 1, function(x) max(x) > 1.5 & max(x)-min(x) > 0.3 & max(x)-min(x) <= as.numeric(opt$diffFreq)))
} else if(!is.null(opt$onlyDiffFreq)) {
    sig_rows <- which(apply(mat, 1, function(x) max(x)-min(x)>as.numeric(opt$diffFreq)))
} else {
    sig_rows <- which(apply(mat, 1, function(x) max(x) > 0 & min(x) < 0 & max(x)-min(x) > as.numeric(opt$diffFreq)))
    #sig_rows <- which(apply(mat, 1, function(x) max(x) > as.numeric(opt$diffFreq) & min(x) < -1*as.numeric(opt$diffFreq) & max(x)-min(x) > 1))
}
#sig_rows <- which(apply(dat, 1, function(x) max(x)>3))
if(!is.null(opt$mustIncludeMotif)) {
    sig_rows = sort(c(sig_rows, grep(paste(unlist(strsplit(as.character(opt$mustInclude),",")),collapse="\\/|"), row.names(dat))))
}

cat(sprintf("%d out of %d motifs passed filter criteria..", length(sig_rows), nrow(mat)))
cat("\n")

if(length(sig_rows)>2) {
    pdf(opt$outPdfFile, height=opt$plotHeight, width=opt$plotWidth)
    #breaks <- as.vector(summary(as.vector(mat[sig_rows,])))
    #len=2
    #breaks1 <- seq(breaks[1], breaks[2], length=len)
    #breaks2 <- seq(breaks[2], breaks[3], length=len)
    #breaks3 <- seq(breaks[3], breaks[4], length=len)
    #breaks4 <- seq(breaks[4], breaks[5], length=len)
    #breaks5 <- seq(breaks[5], breaks[6], length=len)
    #breaks <- c(breaks1[1:length(breaks1)], breaks2[2:length(breaks2)],
    #            breaks3[2:length(breaks3)], breaks4[2:length(breaks4)],
    #            breaks5[2:length(breaks5)])
    if(!is.null(opt$rescale)) {
        breaks <- seq(min(mat[sig_rows,]), max(mat[sig_rows,]), by=as.numeric(opt$breaks))
        #if(length(breaks)<=12) {
        #    myCol <- rev(brewer.pal(length(breaks)-1, "RdBu"))
        #} else {
        #    myCol <- colorRampPalette(c("dodgerblue4","grey97", "sienna3"))(length(breaks)-1)
        #    #myCol <- colorpanel(n=length(breaks)-1,low="blue",mid="#f7f7f7",high="red")
        #    #myCol <- colorpanel(n=length(breaks)-1,low="#2166ac",mid="#f7f7f7",high="#b2182b")
        #}
        myCol <- rev(colorRampPalette(brewer.pal(11, opt$color), bias=as.numeric(opt$colorBias))(length(breaks)-1))
        #heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(15,25), cexCol=1, cexRow=1, breaks=breaks)
        pheatmap(mat[sig_rows,], col=myCol, border_color = NA, cluster_rows=opt$clusterRows, cluster_cols=opt$clusterCols, breaks=breaks, clustering_distance_cols=opt$clusterColsMeasure)
    } else {
        myCol <- rev(colorRampPalette(brewer.pal(11, opt$color), bias=as.numeric(opt$colorBias))(256))
        #myCol <- rev(brewer.pal(11, opt$color))
        #heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(15,25), cexCol=1, cexRow=1)
        pheatmap(mat[sig_rows,], col=myCol, border_color = NA, cluster_rows=opt$clusterRows, cluster_cols=opt$clusterCols, clustering_distance_cols=opt$clusterColsMeasure)
        #heatmap.2(mat, trace="none", col=myCol, margins=c(15,20), cexCol=1, cexRow=1)
    }
    dev.off()

    txtFile=sprintf("%s.txt", opt$outPdfFile)
    mat_df <- as.data.frame(mat)
    colnames(mat_df) <- unlist(lapply(colnames(mat_df), function(x) sprintf("%s (freq_target/freq_bkg)", x)))
    mat_df$diff <- apply(mat_df, 1, function(x) max(x)-min(x))
    dat_df <- as.data.frame(dat)
    colnames(dat_df) <- unlist(lapply(colnames(dat_df), function(x) sprintf("%s (per_target)", x)))
    dat_df$class <- 0
    if(nrow(dat_df[which(apply(dat_df, 1, function(x) max(x)>5)==T),])>0) {
        dat_df[which(apply(dat_df, 1, function(x) max(x)>5)==T),]$class <- 1
    }
    if(nrow(dat_df[which(apply(dat_df, 1, function(x) max(x)>10)==T),])>0) {
        dat_df[which(apply(dat_df, 1, function(x) max(x)>10)==T),]$class <- 2
    }
    if(nrow(dat_df[which(apply(dat_df, 1, function(x) max(x)>20)==T),])>0) {
        dat_df[which(apply(dat_df, 1, function(x) max(x)>20)==T),]$class <- 3
    }
    mat_df$is_sig <- "no"
    mat_df[sig_rows,]$is_sig <- "yes"
    mat_df <- merge(mat_df, dat_df, by=0)
    colnames(mat_df)[1] <- "motif"
    write.table(mat_df[order(-mat_df[,"diff"]),], txtFile, quote=F, append=F, sep="\t", col.names=T, row.names=F)
} else if(length(sig_rows)>=1) {
    txtFile=sprintf("%s.txt", opt$outPdfFile)
    cat("\nNumber of significant motifs are too low for a heatmap. Instead writing results to txtFile\n\n")
    if(length(sig_rows)==1) {
        write(rownames(mat)[sig_rows], file="motif_dynamics.txt", sep="\t")
        write.table(mat[sig_rows,], txtFile, quote=F, append=TRUE, sep="\t")
    } else {
        write.table(mat[sig_rows,], txtFile, quote=F, append=FALSE, sep="\t")
    }
} else {
    cat("\nNo significant motif is found\n")
}
#save.session("test.session")

## old code (v2.0)
#data <- read.table(opt$inFile)
#data$V10 <- apply(data, 1, function(x) binom.test(as.numeric(x[6]), as.numeric(x[7]), as.numeric(x[8])/as.numeric(x[9]), alternative="g")$p.value)
#no_rows=nrow(data)/length(unique(data$V5))
#dat <- matrix(data$V2, nrow=no_rows)
## use log fold change
#mat <- matrix(data$V4, nrow=no_rows)
## use difference between real and background percentages
#mat <- matrix(data$V2-data$V3, nrow=no_rows)
## use p-value computed using binomial test
#mat <- matrix(data$V10, nrow=no_rows)
#colnames(mat) <- as.vector(unique(data$V5))
#row.names(mat) <- data[1:no_rows,1]
#myCol <- brewer.pal(9, "Blues")
#sig_rows <- which(apply(dat, 1, function(x) max(x))>0.33)
#pdf(opt$outPdfFile)
#heatmap.2(mat[sig_rows,], trace="none", col=myCol, margins=c(12,20), cexCol=1, cexRow=1)
#dev.off()

## old code (v1.0)
#data$count <- 1
#dat <- dcast(data, V2 ~ V1)
#dat[is.na(dat)] <- 0
#cols <- ncol(dat)
#dat$mean <- apply(dat[,c(2:cols)],1, mean)
#dat <- dat[with(dat, order(-mean)),]
#mat <- data.matrix(dat)[,c(2:cols)]
#rownames(mat) <- dat$V2
#mycol <- colorpanel(n=99,low="white",high="black")
#colnames(mat) <- gsub("granulocytes", "grn", colnames(mat))
#pdf(opt$outPdfFile)
#heatmap.2(mat, trace="none", col=mycol, margins=c(11,10), cexCol=1)
#dev.off()
