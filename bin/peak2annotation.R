#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inFile"), help="input peak annotation file"),
  make_option(c("-o", "--outFile"), help="output pdf file")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$inFile) | is.null(opt$outFile))) {
  cat("\nProgram: peak2annotation.R (plot peak signal to distance to tss plot)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ggrastr"))
suppressPackageStartupMessages(library("ggpubr"))

## read input file
if(identical(opt$inFile, "stdin")==T) {
    df <- read.table(file("stdin"))
} else {
    df <- read.table(opt$inFile)
}
colnames(df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "dist_to_geneTSS")

## reorganize values to plot
df$dist_to_geneTSS <- log(abs(df$dist_to_geneTSS)+1) * ifelse(df$dist_to_geneTSS<0, -1, 1)
df$annot.type <- "TSS"
df[which(df$dist_to_geneTSS!=0 & abs(df$dist_to_geneTSS) <= log(1000)),]$annot.type <- "proximal"
df[which(df$dist_to_geneTSS!=0 & abs(df$dist_to_geneTSS) > log(1000)),]$annot.type <- "distal"
df$annot.type <- factor(df$annot.type, levels=c("TSS", "proximal", "distal"), ordered=T)

lbl <- paste(sprintf("%0.0f", (c(table(cut(df$dist_to_geneTSS, breaks=c(min(df$dist_to_geneTSS), -log(1000), -1, 0, 1, log(1000), max(df$dist_to_geneTSS))))[c(1:2)],
                                sum(table(cut(df$dist_to_geneTSS, breaks=c(min(df$dist_to_geneTSS), -log(1000), -1, 0, 1, log(1000), max(df$dist_to_geneTSS))))[c(3:4)]),
                                table(cut(df$dist_to_geneTSS, breaks=c(min(df$dist_to_geneTSS), -log(1000), -1, 0, 1, log(1000), max(df$dist_to_geneTSS))))[c(5:6)])*100)/nrow(df)), "%")

## make the plot
if(length(which(!is.na(df$signalValue))) > 0) {
    p <- ggplot(df, aes(dist_to_geneTSS, signalValue)) 
} else {
    p <- ggplot(df, aes(dist_to_geneTSS, score)) 
}

p <- p + geom_point_rast(aes(color=annot.type), alpha=0.2) +
        stat_density_2d(geom = "polygon", contour = TRUE, aes(fill = after_stat(level)), colour = "black", alpha=0.5, show.legend = F) +
        scale_fill_distiller(palette = "Reds", direction = 1) + 
        theme_bw() + xlim(c(-15,15)) + 
        geom_vline(xintercept = c(-log(1000),log(1000)), lty=2) +
        annotate("text", x = c(-10, -5, 0, 5, 10), y = 2, label = lbl,  parse = F, size=3) +  theme_classic2() +
        theme(legend.position="top") + xlab("Distance to closest gene TSS in bp (log)") + ylab("Peak signalValue (Macs2)") +
        labs(color = "Peak position")

ggsave(filename = opt$outFile, plot = p, width = 4, height = 3)
q()
