#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
  make_option(c("-i", "--inFile"), help="input file containing tf enrichment dynamics (can be stdin)"),
  make_option(c("-o", "--outPdfFile"), help="output pdf image file"),
  make_option(c("-c", "--topN"), help="number of top hits to plot from each sample (default=%default)"),
  make_option(c("-f", "--filterPval"), action="store_true", help="filter based on -x, -y and -z"),
  make_option(c("-x", "--pVal"), default=1e-15, help="p-value cutoff (default=%default)"),
  make_option(c("-y", "--minOverlap"), default=50, help="minimum frequency of overlaps in each class (default=%default)"),
  make_option(c("-z", "--minOddsRatio"), default=2, help="minimum odds ratio of at least one class (default=%default)"),
  make_option(c("-F", "--noFilter"), action="store_true", help="no filter, output all hits"),
  make_option(c("-l", "--mustInclude"), help="name of overlap(s) that must be included in the final output (if multiple, separate them by a comma)"),
  make_option(c("-q", "--quaNorm"), action="store_true", help="plot quantile normalized data (default=%default)"),
  make_option(c("-b", "--colorBias"), default=1, help="bias in color (default=%default)"),
  make_option(c("-R", "--clusterRows"), default=T, help="cluster rows in the heatmap (default=%default)"),
  make_option(c("-C", "--clusterCols"), default=T, help="cluster columns in the heatmap (default=%default)"),
  make_option(c("-W", "--plotWidth"), default=15, help="width of the heatmap (default=%default)"),
  make_option(c("-H", "--plotHeight"), default=20, help="height of the heatmap (default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$inFile) | is.null(opt$outPdfFile)) {
	cat("\nProgram: giggleDynAna.R (R script to plot giggle enrichment dynamics)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load libraries
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
#suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load custom functions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
##-----------------------------------------##
## scale df function
##-----------------------------------------##
scale_df <- function(x, y) {
  ## de (default center normalized)
  ## zo (Zero-One); 
  ## zs (Z-Score);
  ## qu (quantile normalization)
  ## rank (rank normalization)
  ## cp (contribution percentage)
  if(missing(y)) { y="de"; }
  if(y=="de") {
    res <- apply(x, 2, function(x) scale(x))
  }
  else if(y=="zo") {
    maxs <- apply(x, 2, max)
    mins <- apply(x, 2, min)
    res <- scale(x, center = mins, scale = maxs - mins)
  } else if(y=="zs") {
    maxs <- apply(x, 2, max)
    mins <- apply(x, 2, min)
    res <- scale(x, center=(maxs+mins)/2, scale=(maxs-mins)/2)
  } else if(y=="qu") {
    res <- normalize.quantiles(as.matrix(x))
  } else if(y=="rank") {
    set.seed(1)
    res <- apply(x, 2, function(x) rank(x, ties.method = "last"))
  } else if (y=="cp") {
    res <- sweep(x, 2, colSums(x), FUN = "/") * 100
  } else {
    res <- x
  }
  
  if(!is.null(colnames(x))) {
    colnames(res) <- colnames(x)
  }
  return(res)
}

##-----------------------------------------##
## heatmap function
##-----------------------------------------##
matrix2Heatmap <- function(x, y=NULL, scale=NULL, col=NULL, bias=NULL, clusterRows=FALSE, clusterCols=FALSE, displayN=FALSE, return_rows=FALSE, title=NULL,
                           columnLabels=NULL, rowAnnot=NULL, legendPos="bottom", showRowNames=T, displayN_customMat=NULL) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(ggpubr))
  suppressPackageStartupMessages(library(tiff))
  
  # x = matrix or data frame
  # y = cluster information (from cluster_matrix)
  if(is.null(y)) {
    y <- rep(1, nrow(x));
  }
  
  if(is.null(scale)) { scale="none"; }
  else if(scale=="col") {
    scale="column"
  }
  
  if(is.null(bias)) { bias=1; }
  if(is.null(displayN_customMat)) { displayN_customMat <- x }
  
  if(scale=="row") {
    x <- t(scale(t(x)))
  } else if(scale=="col" | scale=="column") {
    x <- scale(x)
  }
  if(is.null(col)) { col=colorRampPalette(c("dodgerblue4","grey97", "sienna3"), bias=bias)(250); } 
  else if(length(col) > 1) { col=col; }
  else { col=(colorBrewer2palette(name = col, count = 250, bias=bias)); }
  
  if(is.null(title)) { title=NA; }
  
  if(is.null(columnLabels)) { columnLabels = colnames(x); }
  
  mat <- as.matrix(x)
  df <- as.data.frame(mat)
  df$clusters <- y
  df <- df[order(df$clusters),]
  
  mat <- as.matrix(df[,c(1:ncol(df)-1)])
  #col=colorRamp2(c(seq(min(mat), max(mat), length.out = length(col))), col)
  # if(missing(legend)) {
  #   legend=as.numeric(sprintf("%0.2f", c(seq(min(!is.na(mat)), max(!is.na(mat)), length.out = 10))))
  # }
  if(return_rows==TRUE) {
    return(row_order(draw(Heatmap(mat, cluster_rows=clusterRows, cluster_columns=clusterCols, col=col, use_raster = T, column_title = title,
                                  show_row_names=showRowNames, raster_device="tiff", raster_quality = 2, split = df$clusters, gap = unit(2, "mm"),
                                  rect_gp = gpar(col = "black")))))
  } else if(displayN==FALSE) {
    grid.grabExpr(draw(Heatmap(mat, cluster_rows=clusterRows, cluster_columns=clusterCols, col=col, use_raster = T, column_title = title,
                               show_row_names=showRowNames, raster_device="tiff", raster_quality = 2, split = df$clusters, gap = unit(2, "mm"),
                               border=F, column_labels=columnLabels, right_annotation = rowAnnot, heatmap_legend_param = list(direction = "horizontal")),
                       heatmap_legend_side = legendPos, annotation_legend_side = legendPos))
  } else {
    grid.grabExpr(draw(Heatmap(mat, cluster_rows=clusterRows, cluster_columns=clusterCols, col=col, use_raster = T, column_title = title,
                               show_row_names=showRowNames, raster_device="tiff", raster_quality = 2, split = df$clusters, gap = unit(2, "mm"),
                               cell_fun = function(j, i, x, y, width, height, fill) { grid.text(sprintf("%s", displayN_customMat[i, j]), x, y, 
                                                                                                gp = gpar(fontsize = 8)) },
                               rect_gp = gpar(col = "black"), column_labels=columnLabels, right_annotation = rowAnnot, heatmap_legend_param = list(direction = "horizontal")),
                       heatmap_legend_side = legendPos, annotation_legend_side = legendPos))
  }
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## read input data
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(identical(opt$inFile, "stdin")==T) { 
    data <- read.table(file("stdin"), header=T)
} else {
    data <- read.table(opt$inFile, header=T)
}
# data <- read.table("~/project/misc/julia//multiClassGiggleAna/giggleDynAna/GIGGLE_ENRICHMENT_HOMER.TXT", header=T)
# opt <- NULL; opt$topN <- 10; opt$pVal <- 1e-15; opt$minOverlap <- 50; opt$minOddsRatio <- 8; 
# opt$clusterRows <- T; opt$clusterCols <- T; opt$colorBias <- 1; opt$topN <- 10

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## organize data frame containing enrichment results
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
no_rows=nrow(data)/length(unique(data$id))
data$file <- gsub(".bed.gz", "", gsub("^.*/", "", data$file))

df <- lapply(c("overlaps", "odds_ratio", "combo_score", "fishers_two_tail"), function(x) {
            t <- as.data.frame(matrix(data[,x], nrow=no_rows))
            colnames(t) <- unique(data$id)
            row.names(t) <- unique(data$file)
            t
        })
names(df) <- c("overlaps", "odds_ratio", "combo_score", "fishers_two_tail")

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## select hits that have enough variation in combo_score & odds_ratio across conditions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$noFilter)) {
    sig_names <- intersect(row.names(df$combo_score)[which(rowSds(normalize.quantiles(as.matrix(df$combo_score))) > summary(rowSds(normalize.quantiles(as.matrix(df$combo_score))))[3])],
                           row.names(df$odds_ratio)[which(rowSds(normalize.quantiles(as.matrix(df$odds_ratio))) > summary(rowSds(normalize.quantiles(as.matrix(df$odds_ratio))))[3])])
    cat(sprintf("%d out of %d hits passed filter criteria (>median sd of combo_score + odds_ratio)..", length(sig_names), nrow(df$overlaps)))
    cat("\n")
} else {
    sig_names <- row.names(df$combo_score)
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## select hits that have significant pVal + minOverlap + minOddsRatio
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(!is.null(opt$filterPval) & is.null(opt$noFilter)) {
  ## identify significant overlaps
  TMP <- length(sig_names)
  sig_names <- sig_names[sig_names %in% row.names(df$combo_score[which(rowMin(as.matrix(df[["fishers_two_tail"]])) < opt$pVal &
                                                                      rowMax(as.matrix(df[["overlaps"]])) > opt$minOverlap &
                                                                      rowMax(as.matrix(df[["odds_ratio"]])) >= opt$minOddsRatio),])]
  cat(sprintf("%d out of %d hits passed filter criteria (pVal + minOverlap + minOddsRatio)..", length(sig_names), TMP))
  cat("\n")
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## identify rows corresponding to significant hits in sig_names
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
sig_rows <- which(row.names(df$combo_score) %in% sig_names)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## add overlaps which must be included
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(!is.null(opt$mustInclude)) {
  sig_rows = c(sig_rows, which(row.names(df$overlaps) %in% unlist(strsplit(opt$mustInclude, ","))))
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## quantile normalize the data
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(!is.null(opt$quaNorm)) {
  df_sig <- cbind(df[["overlaps"]][sig_rows,] %>% mutate(name = row.names(df[["overlaps"]][sig_rows,])) %>% reshape2::melt(),
                  df[["odds_ratio"]][sig_rows,] %>% as.matrix() %>% scale_df(y="qu") %>% `colnames<-` (colnames(df[["combo_score"]][sig_rows,])) %>% `rownames<-` (row.names(df[["combo_score"]][sig_rows,])) %>% as.data.frame() %>% mutate(name = row.names(df[["odds_ratio"]][sig_rows,])) %>% reshape2::melt() %>% pull(value),
                  df[["fishers_two_tail"]][sig_rows,] %>% mutate(name = row.names(df[["fishers_two_tail"]][sig_rows,])) %>% reshape2::melt() %>% pull(value),
                  df[["combo_score"]][sig_rows,] %>% as.matrix() %>% scale_df(y="qu") %>% `colnames<-` (colnames(df[["combo_score"]][sig_rows,])) %>% `rownames<-` (row.names(df[["combo_score"]][sig_rows,])) %>% as.data.frame() %>% mutate(name = row.names(df[["combo_score"]][sig_rows,])) %>% reshape2::melt() %>% pull(value))
} else {
  df_sig <- cbind(df[["overlaps"]][sig_rows,] %>% mutate(name = row.names(df[["overlaps"]][sig_rows,])) %>% reshape2::melt(),
                  df[["odds_ratio"]][sig_rows,] %>% mutate(name = row.names(df[["odds_ratio"]][sig_rows,])) %>% reshape2::melt() %>% pull(value),
                  df[["fishers_two_tail"]][sig_rows,] %>% mutate(name = row.names(df[["fishers_two_tail"]][sig_rows,])) %>% reshape2::melt() %>% pull(value),
                  df[["combo_score"]][sig_rows,] %>% mutate(name = row.names(df[["combo_score"]][sig_rows,])) %>% reshape2::melt() %>% pull(value))
}

colnames(df_sig) <- c("name", "class", "overlaps", "odds_ratio", "pvalue", "combo_score")
#df_sig$odds_ratio <- log(df_sig$odds_ratio)
df_sig$pvalue <- -log10(df_sig$pvalue)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## filter out overlaps corresponding to which standard deviation is zero
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
# df_sig <- df_sig[which(df_sig$name %in% row.names(as.data.frame(matrix(df_sig[,"odds_ratio"], nrow=length(unique(df_sig$name))) %>% as.data.frame() %>% `colnames<-` (unique(df_sig$class)) %>% `rownames<-` (unique(df_sig$name)))[which(rowSds(matrix(df_sig[,"odds_ratio"], nrow=length(unique(df_sig$name))))>0),])),]

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## make the plot
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(nrow(df_sig)>2) {
  #ggboxplot(data, x="id", y="overlaps")
  
  p1 <- ggplot(df_sig, aes(x = class, y = name)) +
    geom_point(aes(colour=combo_score, size=odds_ratio)) +
    #geom_point(aes(colour=pvalue,size=GeneDensity*100)) +
    #scale_colour_gradient(high="brown", low="yellow", name="p-value") +
    #scale_colour_gradient(low="#00441b", high="#ccece6", name="p-value") +
    scale_colour_gradient(high="#00441b", low="#99d8c9", name="log(combo_score)") +
    #breaks=c(1e-20, 1e-15, 1e-10, 1e-5, 0.01, 0.05)) +
    #scale_size(range=c(1,10)) +
    #scale_size(range=c((min(data_sig$GeneDensity)*100),(max(data_sig$GeneDensity*100))), breaks=waiver(), labels=waiver()) +
    #scale_size_area(breaks=seq(0,100,by=10), max_size=20) +
    theme_bw(base_size=15) +
    theme(text=element_text(size=10), axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.position = "bottom")
  
  mat2plot <- matrix(df_sig[,"odds_ratio"], nrow=length(unique(df_sig$name))) %>% as.data.frame() %>% `colnames<-` (unique(df_sig$class)) %>% `rownames<-` (unique(df_sig$name)) %>% round(2)
  if(!is.null(opt$topN)) {
      # topN <- intersect(rowSds(as.matrix(mat2plot)) %>% sort(decreasing = T) %>% head(as.numeric(opt$topN)) %>% names,
      #           unique(unlist(lapply(colnames(mat2plot), function(x) mat2plot[order(-mat2plot[,x]),] %>% head(as.numeric(opt$topN)) %>% row.names))))
      topN <- unique(unlist(lapply(colnames(mat2plot), function(x) mat2plot[order(-mat2plot[,x]),] %>% head(as.numeric(opt$topN)) %>% row.names)))
      #topN <- rowSds(as.matrix(mat2plot)) %>% sort(decreasing = T) %>% head(as.numeric(opt$topN)) %>% names
      cat(sprintf("%d out of %d hits passed filter criteria (topN)..", length(topN), nrow(mat2plot)))
      cat("\n")
      mat2plot <- mat2plot[which(row.names(mat2plot) %in% topN),]
  }
  p2 <- matrix2Heatmap(mat2plot, 
                       scale="row", clusterRows = opt$clusterRows, clusterCols = opt$clusterCols, bias=opt$colorBias, displayN = T)

  p2 <- ggarrange(p2, ncol=2)
  #ggsave(opt$outPdfFile, ggarrange(plotlist = list(p1, p2), nrow=1, ncol=2, labels = c("A)", "B)")), height=opt$plotHeight, width=opt$plotWidth, device="pdf")
  ggsave(opt$outPdfFile, marrangeGrob(grobs = list(p2, p1), nrow=1, ncol=1, labels=c("A", "B")),  height=opt$plotHeight, width=opt$plotWidth, device="pdf")
  df[["sig"]] <- df_sig
} else {
  cat("\nNo enrichment is found\n")
  df[["sig"]] <- NA
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## write output file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
write.table(mat2plot, sprintf("%s.txt", opt$outPdfFile), sep="\t", quote = F, row.names = T, col.names = T)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save objects to .RDS file
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
saveRDS(list(df=df, mat2plot=mat2plot), file=sprintf("%s.Rds", opt$outPdfFile))

q()
