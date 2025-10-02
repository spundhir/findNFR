#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## parse command line arguments
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
option_list <- list(
  make_option(c("-i", "--inFile"), help="input file containing peak score (macs2 output) [format: chr start end name score strand signalValue]"),
  make_option(c("-o", "--outFile"), help="output pdf file"),
  make_option(c("-g", "--genome"), default="mm10", help="genome (mm10 or hg38; default=%default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## check, if all required arguments are given
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
if(is.null(opt$inFile) | is.null(opt$outFile)) {
  cat("\nProgram: peak2spatialSignal.R (plot peak signal classified based on distance to TSS)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load packages
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrastr"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## load custom functions
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
suppressPackageStartupMessages(source("~/software/myScripts/rScripts/FUNCTIONS.R"))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## initialize input files based on genome
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## read input file
if(identical(opt$inFile, "stdin")==T) {
    peaks <- read.table(file("stdin"))[,c(1:7)]
} else {
    peaks <- read.table(opt$inFile)[,c(1:7)]
}
# peaks <- read.table("~/data/00_ALL_CHIP-SEQ_RAW/MLL-AF9/six1_on_peaks.bed", header=F)[,c(1:7)]
colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand", "signalValue")

## reformat peak file to ensure correct format of signal values
if(length(which(!is.na(peaks$signalValue))) > 0 & is.numeric(peaks$signalValue)) {
  peaks$signalValue <- log2(peaks$signalValue+1)
} else if(length(which(!is.na(peaks$score))) > 0 & is.numeric(peaks$score)) {
  peaks$signalValue <- log2(peaks$score+1)
} else {
  cat("ERROR: no score column found\n");
  print_help(parser)
}

GENOME_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s", opt$genome), intern=T)
TSS_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -t", opt$genome), intern=T)
TISSUESPECIFICITY_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -S", opt$genome), intern=T)
GENEDENSITY_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -M", opt$genome), intern=T)
CPG_FILE <- system(sprintf("source ~/.bashrc && initialize_genome -g %s -c", opt$genome), intern=T)

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## start analysis
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## closest/target gene information
# df <- bed2closestCoor(peaks, read.table(pipe(sprintf("grep -w protein_coding %s | cut -f 1-6", TSS_FILE))), strandAware = T)[,c(1:7,11,14)]
# colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "closest_gene", "dist_to_closest_gene")
df <- merge(peaks, linkDHS2Genes(bed2window(peaks, win = 250, flank_to_tss = F), genome = opt$genome, useLoops = F, strandAware = T)[,c("name", "closest_gene", "dist_to_closest_gene")],
            by.x="name", by.y="name")
df <- df[,c(2:4,1,5:ncol(df))]
df <- merge(df, linkDHS2Genes(bed2window(peaks, win = 250, flank_to_tss = F), genome = opt$genome, useLoops = T, minoverlap = 100)[,c("name", "target_gene", "cInteraction_score", "dist_to_target_gene")],
  by.x="name", by.y="name")
df <- df[,c(2:4,1,5:ncol(df))]

## closest/target gene tissue specificity information
tmp <- read.table(TISSUESPECIFICITY_FILE, header=T)[,c("external_gene_name", "tau")] %>% dplyr::filter(!is.na(external_gene_name))
df$closest_gene_tau <- tmp$tau[match(df$closest_gene, tmp$external_gene_name)]
t <- df %>% separate_longer_delim(target_gene, ",")
t$target_gene_tau <- round(tmp$tau[match(t$target_gene, tmp$external_gene_name)],7)
df <- merge(df, aggregate(target_gene_tau ~ name, t, toString) %>% mutate(target_gene_tau = gsub("\\s+", "", target_gene_tau)), by.x="name", by.y="name", all.x=T)
df <- df[,c(2:4,1,5:ncol(df))]

## closest gene density information
df <- merge(df, read.table(GENEDENSITY_FILE, header=T)[,c("name", "geneDensityScore", "geneLength", "geneDensityClass")], by.x="closest_gene", by.y="name")
df$geneDensityClass <- factor(df$geneDensityClass)
df <- df[,c(2:8,1,9:ncol(df))]

## CpG island information
df$cpgOverlap <- ifelse(df$name %in% bed2overlap(df, CPG_FILE)$name_q, "CpG", "nonCpG")

df$dist_to_closest_gene <- log(abs(df$dist_to_closest_gene)+1) * ifelse(df$dist_to_closest_gene<0, -1, 1)
df$annot.type <- "promoter"
df[which(abs(df$dist_to_closest_gene) > log(1000) & abs(df$dist_to_closest_gene) <= log(50000)),]$annot.type <- "proximal"
df[which(abs(df$dist_to_closest_gene) > log(50000)),]$annot.type <- "distal"
df$annot.type <- factor(df$annot.type, levels=c("promoter", "proximal", "distal"), ordered=T)
# table(df$annot.type)

## class defined based on distance to closest gene
df$distClass <- cut2(abs(df$dist_to_closest_gene), cuts=log(c(0, 1000, 2500, 5000, 10000, 50000, 100000, 500000, 1000000, 40000000)), levels.mean = T)
levels(df$distClass) <- sprintf("%s_%s", seq(1,length(levels(df$distClass))), round(as.numeric(gsub("\\s+", "", levels(df$distClass)))))
# table(df$distClass)

# df %>% separate_longer_delim(dist_to_target_gene, ",") %>% 
#   mutate(dist_to_target_gene = log(as.numeric(dist_to_target_gene)+1)) %>% 
#   mutate(dist_to_closest_gene = (abs(dist_to_closest_gene)+1)) %>% 
#   ggscatter(x="dist_to_target_gene", y="dist_to_closest_gene")

brk <- c(min(df$dist_to_closest_gene), -log(50000), -log(1000), 0, log(1000), log(50000), max(df$dist_to_closest_gene))
lbl_abs <- (table(cut(df$dist_to_closest_gene, 
                      breaks=c(min(df$dist_to_closest_gene), -log(50000), -log(1000), 0, log(1000), log(50000), max(df$dist_to_closest_gene)), 
                      labels=c(min(df$dist_to_closest_gene), -log(50000), 0, 0, log(50000), max(df$dist_to_closest_gene))))) %>% as.vector
lbl_per <- gsub("\\s", "", paste(sprintf("%0.1f", (lbl_abs*100)/nrow(df)), "%"))

## make the plot
MAX_VAL=round(max(df$signalValue),0)
p1 <- ggplot(df, aes(dist_to_closest_gene, signalValue)) +
            geom_point_rast(aes(color=annot.type), alpha=0.2) +
            stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "black", alpha=0.5, show.legend = F) +
            scale_fill_distiller(palette = "Reds", direction = 1) + 
            #xlim(c(-15,15)) + 
            geom_vline(xintercept = brk, lty=2) +
            scale_x_continuous(breaks=brk, labels = as.numeric(sprintf("%0.0f", brk))) +
            #annotate("text", x = c(-10, -5, 0, 5, 10), y = 2, label = lbl,  parse = F, size=3) +  theme_classic2() +
            annotate("text", x = c(-12.5, -9, 0, 9, 12.5), y = MAX_VAL, label = lbl_abs,  parse = F, size=3) +  
            annotate("text", x = c(-12.5, -9, 0, 9, 12.5), y = MAX_VAL-(round((MAX_VAL * 5/100),1)), label = paste0("(", lbl_per, ")"),  parse = F, size=3) +
            theme_classic2() + scale_color_manual(values=c("#440154", "#21908c", "#fde725")) +
            theme(legend.position="top") + xlab("Distance to closest gene TSS in bp (log)") + ylab("Peak signalValue (Macs2; log2)") +
            labs(color = "Peak position")
p2 <- ggplot(df, aes(geneDensityScore, signalValue)) +
  geom_point_rast(aes(color=annot.type), alpha=0.2) +
  stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "black", alpha=0.5, show.legend = F) +
  scale_fill_distiller(palette = "Reds", direction = 1) + 
  theme_classic2() + scale_color_manual(values=c("#440154", "#21908c", "#fde725")) +
  theme(legend.position="none") + xlab("Gene density score") + ylab("Peak signalValue (Macs2; log2)") +
  labs(color = "Peak position") +
  theme(plot.margin = margin(t = 50, r = 5, b = 5, l = 5))
## peaks proximal to genes are located in gene dense regions (circular argument)
# p3 <- ggplot(df, aes(x=abs(dist_to_closest_gene), y=log(geneDensityScore))) + geom_point(aes(color=annot.type)) + theme_bw() +
#         xlab("Distance b/w peak to closest gene TSS in bp (log)") +
#         ylab("Gene Density Score (log)")
# p3 <- ggplot(reshape2::melt(table(df[,c("annot.type", "geneDensityClass")])), aes(x=as.factor(geneDensityClass), y=value, fill=annot.type)) + geom_bar(stat = "identity", position = "fill") +
#   scale_y_continuous(labels = scales::percent) + theme_classic() + scale_fill_manual(values=c("#440154", "#21908c", "#fde725")) +
#   xlab("geneDensityClass") + ylab("Density") + labs(fill = "Peak position") + theme(legend.position="top")
p3 <- ggbarplot(as.data.frame(table(df$geneDensityClass)), x="Var1", y="Freq", fill="Var1") + xlab("geneDensityClass") + ylab("# peaks") + 
  scale_fill_manual(values=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(10))) + theme(legend.position="none") +
  theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))
# p4 <- ggecdf(df, x="signalValue", color="geneDensityClass") + 
#   scale_color_manual(values=rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdYlBu"))(10))) + xlab("Peak signalValue (Macs2; log2)") + ylab("ECDF") +
#   theme(legend.position="none")
P <- ggarrange(p1, ggarrange(p2, p3, nrow=2, ncol=1, labels=c("B", "C")), nrow=1, ncol=2, labels=c("A", ""), widths = c(3,1))

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
## save output files
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
ggsave(filename = opt$outFile, plot = P, width = 9, height = 6)
write.table(df, "", sep="\t", col.names = T, row.names = F, quote = F)
q()