#!/usr/bin/env Rscript
#suppressPackageStartupMessages(library("optparse", lib.loc = "/home/users/sachin/R/x86_64-redhat-linux-gnu-library/2.15"))
suppressPackageStartupMessages(library("optparse"))

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--inputDeeptoolsMatrix"), help="input deepPlot matrix file"),
  make_option(c("-o", "--outputFile"), help="output deepPlot file"),
  make_option(c("-j", "--groupAnno"), help="group annotation"),
  make_option(c("-G", "--plotHeight"), help="plot height in cm"),
  make_option(c("-W", "--plotWidth"), help="plot width in cm")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if((is.null(opt$inputDeeptoolsMatrix) | is.null(opt$outputFile))) {
  cat("\nProgram: deepPlot.R (plot deepPlot matrix using ggplot)\n")
  cat("Author: BRIC, University of Copenhagen, Denmark\n")
  cat("Version: 1.0\n")
  cat("Contact: pundhir@binf.ku.dk\n");
  print_help(parser)
  q()
}

suppressPackageStartupMessages(library("jsonlite"))
suppressPackageStartupMessages(library("ggpubr"))
suppressPackageStartupMessages(library("doBy"))
suppressPackageStartupMessages(library("ggplot2"))

## read input deepTools matrix file
df <- read.table(pipe(sprintf("zgrep -v '@' %s", opt$inputDeeptoolsMatrix)))
mat <- as.matrix(df[,-c(1:6)])

## read JSON paramters
h <- scan(opt$inputDeeptoolsMatrix, n=1, sep="\n", what=character())
params <- fromJSON(txt=gsub("@", "", h))
if(!is.null(opt$groupAnno)) {
    params$group_labels <- unlist(strsplit(opt$groupAnno, ","))
}

## plot deepPlot matrix using ggplot
p <- lapply(seq(1:(length(params$sample_boundaries)-1)), function(i) {
            df.melt <- as.data.frame(mat[,c((params$sample_boundaries[i]+1):params$sample_boundaries[i+1])])
            df.melt$class <- unlist(lapply(seq(1:(length(params$group_boundaries)-1)), function(j) {
                rep(params$group_labels[j], params$group_boundaries[j+1]-params$group_boundaries[j] )
            }))
            df.melt <- reshape2::melt(df.melt)
            #ggline(df.melt, x="variable", y="value", add="mean", color = "class", numeric.x.axis=T, title = params$sample_labels[i], plot_type = "l") + ylab("") + xlab("") +
            #        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            #        scale_x_discrete(breaks=c(1,100,200))
            df.melt <- summaryBy(value ~ variable + class, df.melt, FUN=function(x) c(mean(x), sd(x)/sqrt(length(x))))
            colnames(df.melt) <- c("variable", "class", "mean", "se")
            ggplot(df.melt, aes(as.numeric(variable), mean, color=class)) + geom_line() + 
                    geom_ribbon(aes(y = mean, ymin = mean - se, ymax = mean + se), alpha=0.1, linetype = 0) +
                    theme_bw() + 
                    scale_x_continuous(breaks=c(1,(length(unique(df.melt$variable))/2),length(unique(df.melt$variable))), 
                                        labels = c(params$upstream[i]*-1, params$`ref point`[i], params$downstream[i])) +
                    theme(axis.title.x=element_blank(), legend.position="top", legend.title = element_blank(), 
                                        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                    geom_vline(xintercept = (length(unique(df.melt$variable))/2), lty=2) +
                    ylab("") + ggtitle(params$sample_labels[i])
    })

ggarrange(plotlist = p) %>% ggexport(filename = opt$outputFile)
q()
