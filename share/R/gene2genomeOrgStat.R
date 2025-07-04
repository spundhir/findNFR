#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inDir"), help="input directory containing output from gene2genomeOrgStat"),
	make_option(c("-j", "--inFile"), help="input file containing output from gene2genomeOrgStat")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inDir)) {
	cat("\nProgram: gene2genomeOrgStat.R (R script to compute genome wide organization statistics of genes))\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tiff))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(tiff))

#################################
## functions
#################################
cluster_info_matrix <- function(x, maxX=NULL, title=NULL) {
    # x: matrix or data frame
    # maxX: maximum x limit (default: 30)
    if(is.null(maxX)) { maxX <- 30; }
    if(is.null(title)) { title <- ""; }
    set.seed(1)
    mat <- as.matrix(x)
    wss <- (nrow(mat)-1)*sum(apply(mat,2,var))
    for (i in 2:maxX)
    wss[i] <- sum(kmeans(mat, centers=i)$withinss)
    plot(1:maxX, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares", main=title)
}

#################################
## step-1 analysis on inDir (distribution of distance between genes on each chromosome)
#################################
p <- lapply(list.files(opt$inDir, pattern="chr", full.names = T), function(FILE) {
    df <- read.table(pipe(sprintf("cut -f 4,12,17 %s", FILE)))
    test <- dcast(data.table::as.data.table(df, 2000), V1 ~ V2, value.var="V3")
    row.names(test) <- test$V1
    test$V1 <- NULL
    test[is.na(test)] <- 0
    # matrix2pheatmap(test, clusterRows = T, clusterCols = T, bias = 3)
    # head(test[1:10,1:10])

    #summary(log(apply(test, 1, function(x) min(x[x>0]))))
    lowLimit <- quantile(log(apply(test, 1, function(x) min(x[x>0]))))[2]
    upLimit <- quantile(log(apply(test, 1, function(x) min(x[x>0]))))[4]

    ggplot(reshape2::melt(log(apply(test, 1, function(x) min(x[x>0])))), aes(value)) + geom_density(fill="#9ebcda") + theme_bw() + xlab("") + xlim(c(0,20)) +
        geom_vline(xintercept = c(lowLimit, upLimit), lty=2) +
        geom_text(aes(lowLimit, 0.1, label = as.numeric(sprintf("%0.2f", lowLimit)), hjust = 2), size = 6) +
        geom_text(aes(upLimit, 0.1, label = as.numeric(sprintf("%0.2f", upLimit)), hjust = -1), size = 6) +
        ggtitle(sprintf("%s (N=%s)", gsub("^.*\\.", "", FILE), nrow(test))) + ylab("Density") + xlab("Distance between genes (log)")
})

#################################
## step-2 analysis on inFile (statistics on distance between genes and DHS)
#################################
df <- data.table::fread(opt$inFile)
colnames(df) <- c("chr", "start", "end", "gene", "score", "strand", "gene2geneDist", "gene_length", "dhs", "dhs2geneDist")
df$dhs_class <- ifelse(df$dhs2geneDist==0, "intragenic", "intergenic")
df$gene2geneDistClass <- cut2((df$gene2geneDist), g=10, levels.mean=T)
levels(df$gene2geneDistClass) <- cut2((df$gene2geneDist), g=10, onlycuts = T)[2:length(cut2((df$gene2geneDist), g=10, onlycuts = T))]
# df$gene2geneDistClass <- cut(df$gene2geneDist, breaks = c(0, 1000, 10000, 20000, 50000, 100000, 200000, 500000, 100000000))
df$geneDensityClass <- df$gene2geneDistClass
levels(df$geneDensityClass) <- seq(length(levels(df$geneDensityClass)), 1)

# p1 <- ggplot(unique(df[,c("gene","gene2geneDist")]), aes(log(gene2geneDist))) + geom_density(fill="#9ebcda") + theme_bw() + 
#   xlim(c(0,15)) + xlab("Distance between genes (bp; log)")

p1 <- ggplot(df, aes(log(dhs2geneDist+1))) + geom_density(fill="#9ebcda") + theme_bw() + xlim(c(0,15)) + 
xlab("Distance between DHS and closest gene (bp; log+1)") + ylab("Density") +
annotate(geom="text", x=5, y=0.4, label=sprintf("# Genes=%d", length(unique(df$gene))), color="red") +
annotate(geom="text", x=5, y=0.35, label=sprintf("# DHS=%d", length(unique(df$dhs))), color="red")

t <- unlist(lapply(levels(df$gene2geneDistClass), function(x) length(unique(df[which(df$gene2geneDistClass==x),]$dhs))))
names(t) <- levels(df$gene2geneDistClass)
t <- reshape2::melt(t)
t$gene2geneDistClass <- factor(row.names(t), levels=levels(df$gene2geneDistClass), ordered=T)
p2 <- ggplot(t, aes(gene2geneDistClass, value)) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(size=10, angle=45)) + 
ylab("# of DHS") + xlab("Distance between genes (bp)")

t <- unlist(lapply(levels(df$gene2geneDistClass), function(x) length(unique(df[which(df$gene2geneDistClass==x),]$gene))))
names(t) <- levels(df$gene2geneDistClass)
t <- reshape2::melt(t)
t$gene2geneDistClass <- factor(row.names(t), levels=levels(df$gene2geneDistClass), ordered=T)
p3 <- ggplot(t, aes(gene2geneDistClass, value)) + geom_bar(stat="identity") + theme_bw() + theme(axis.text.x = element_text(size=10, angle=45)) + 
ylab("# of genes") + xlab("Distance between genes (bp)")

p4 <- ggplot(df, aes(gene2geneDistClass, log(dhs2geneDist+1))) + geom_boxplot() + theme_bw() + 
theme(axis.text.x = element_text(size=10, angle=45)) + ylab("Distance between DHS and closest gene (bp; log+1)") +
xlab("Distance between genes (bp)")

# t <- lapply(levels(df$class), function(x) as.vector(unlist(table(df[which(df$class==x),]$gene))))
# names(t) <- levels(df$class)
# t <- reshape2::melt(t)
# p5 <- ggplot(t, aes(L1, value)) + geom_boxplot() + theme_bw() + theme_bw() + theme(axis.title.x=element_blank(), axis.text.x = element_text(size=10, angle=45)) + ylab("# of DHS per gene")

t <- lapply(levels(df$gene2geneDistClass), function(x) (table(df[which(df$gene2geneDistClass==x),]$dhs_class)))
names(t) <- levels(df$gene2geneDistClass)
t <- reshape2::melt(t)[,c(3,1,2)]
colnames(t) <- c("xlabels", "class", "freq")
t$xlabels <- factor(t$xlabels, levels=as.character(unique(t$xlabels)), ordered = T)
total <- as.vector(unlist(ddply(t, .(xlabels), summarise, X2=sum(freq))[2]))
t$total <- rep(total, length(t[,1])/length(total))
t$density <- t$freq/t$total
p5 <- ggplot(t, aes(xlabels, density, fill=factor(class), label = freq)) + 
        geom_bar(stat="identity", position = "stack") +
        ggtitle("title") +
        geom_text(size = 3, position = position_stack(vjust = 0.5)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme_bw() + theme(axis.text.x = element_text(size=10, angle=45)) + xlab("Distance between genes (bp)")

p6 <- ggplot(unique(df[,c("gene","gene2geneDistClass","gene_length")]), aes(gene2geneDistClass, log(gene_length))) + geom_boxplot() + theme_bw() + 
theme(axis.text.x = element_text(size=10, angle=45)) + ylab("gene length (bp; log)") + xlab("Distance between genes (bp)")

# table(unique(df[,c("gene", "geneDensityClass")])$geneDensityClass)

#################################
## step-3 analysis on inDir (heatmap, clustering of genes based on their distance to each other)
#################################
h <- lapply(list.files(opt$inDir, pattern="chr", full.names = T), function(FILE) {
    mat <- read.table(pipe(sprintf("cut -f 4,12,17 %s", FILE)))
    mat <- dcast(data.table::as.data.table(mat, 2000), V1 ~ V2, value.var="V3")
    row.names(mat) <- mat$V1
    mat$V1 <- NULL
    mat[is.na(mat)] <- 0
    mat$Row.names <- row.names(mat)
    mat <- merge(mat, unique(df[,c("gene", "geneDensityClass")]), by.x="Row.names", by.y="gene", all.x="T")
    row.names(mat) <- mat$Row.names
    mat$Row.names <- NULL

    col_cluster <- rev(colorRampPalette(brewer.pal(length(unique(mat$geneDensityClass)), "RdBu"), bias=1)(length(unique(mat$geneDensityClass))))
    names(col_cluster) <- c(levels(mat$geneDensityClass)[levels(mat$geneDensityClass) %in% unique(mat$geneDensityClass)], "UA")

    ra <- rowAnnotation(na_col=NA,
        geneDensityClass = mat$geneDensityClass,
        col=list(geneDensityClass = col_cluster)
    )

    grid.grabExpr(draw(Heatmap(as.matrix(mat[,-ncol(mat)]), cluster_rows=T, cluster_columns=T,
        col=colorRampPalette(c("dodgerblue4","grey97", "sienna3"), bias=1.5)(250), use_raster=T,
        show_row_names=F, show_column_names=F, raster_device="tiff", raster_quality = 10, right_annotation = ra,
        column_title = sprintf("%s (N=%s)", gsub("^.*\\.", "", FILE), nrow(mat)))))
})

#save.session("test.session"); q()
#################################
## step-4 output plots
#################################
pdf(sprintf("%s/gene2genomeOrgStat.pdf", opt$inDir), width=20, height=15)
    ggarrange(plotlist = p)
    ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol=3)
    ggarrange(plotlist=h[1:4], nrow = 2, ncol=2)
    ggarrange(plotlist=h[5:8], nrow = 2, ncol=2)
    ggarrange(plotlist=h[9:12], nrow = 2, ncol=2)
    ggarrange(plotlist=h[13:16], nrow = 2, ncol=2)
    ggarrange(plotlist=h[17:20], nrow = 2, ncol=2)
    ggarrange(plotlist=h[21:length(h)], nrow = 2, ncol=2)
garbage <- dev.off()

#################################
## step-5 output file containing gene density information
#################################
write.table(unique(df[,c("chr", "start", "end", "gene", "score", "strand", "gene2geneDist", "gene_length", "gene2geneDistClass", "geneDensityClass")]), "", sep="\t", quote=F, row.names = F, col.names=T)

q()
