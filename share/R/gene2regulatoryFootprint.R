#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("session"))

## parse command line arguments
option_list <- list(
	make_option(c("-i", "--inFile"), help="input file containing score values (can be stdin) (format: gene score(s) distance(s))"),
    make_option(c("-s", "--scoreCol"), default=2, help="column in input file that contains score information. If multiple, separate them by a comma (default: %default)"),
    make_option(c("-d", "--distanceCol"), default=3, help="column in input file that contains distance information (default: %default)"),
    make_option(c("-w", "--weightModel"), default=2024, help="calculate weight using equation from 2011, 2016 or 2024 paper  (default: %default)")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)
opt <- parse_args(parser)

## check, if all required arguments are given
if(is.null(opt$inFile)) {
	cat("\nProgram: gene2regulatoryFootprint.R (R script to compute regulatory footprint of input CRE on genes)\n")
	cat("Author: BRIC, University of Copenhagen, Denmark\n")
	cat("Version: 1.0\n")
	cat("Contact: pundhir@binf.ku.dk\n");
	print_help(parser)
	q()
}

## load libraries
suppressPackageStartupMessages(library(session))
suppressPackageStartupMessages(library(plyr))

## load data
if(identical(opt$inFile, "stdin")==T) {
    df <- read.table(file("stdin"))
} else {
    df <- read.table(opt$inFile)
}

opt$scoreCol <- as.numeric(opt$scoreCol)
opt$distanceCol <- as.numeric(opt$distanceCol)

#############################################
## function to compute total signal
#############################################
compute_total_signal <- function(signal) {
    if(is.na(signal)) {
        total_signal <- 0
    } else {
        total_signal <- sum(as.numeric(unlist(strsplit(as.character(signal), ","))))
    }
    return(total_signal)
}

#############################################
## function to compute weighted signal (2011 paper)
## Idea source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3610570/
## Based on the characterization of higher order chromatin interactions and our preliminary analysis, we calculated the regulatory potential for a given gene, Sg, as the sum of the nearby binding sites weighted by the distance from each site to the TSS of the gene: 𝑆𝑔=∑𝑘𝑖=1𝑒−(0.5+4Δ𝑖), where k is the number of binding sites within 100 kb of gene g and Δi is the distance between site i and the TSS of gene g normalized to 100 kb (e.g., 0.5 for a 50 kb distance). This equation models the influence of each binding site on gene regulation as a function that decreases monotonically with increasing distance from the TSS. The shape of this function approximates empirical observations of the distance between binding sites and differentially expressed genes in multiple ChIP-seq experiments. The constant in the equation enables the exponential function to adopt more flexible shapes, and 0.5 was derived to better fit ChIA-PET and Hi-C data.
#############################################
compute_weighted_signal_2011 <- function(signal, distance) {
    if(is.na(distance)) {
        weighted_signal <- 0
    } else {
        all_signal <- as.numeric(unlist(strsplit(as.character(signal),",")))
        all_distance <- abs(as.numeric(unlist(strsplit(as.character(distance),","))))/100000

        e <- exp(1)
        weighted_signal=0
        for(i in 1:length(all_distance)) {
            weight=(e^-(0.5+(4*all_distance[i])))
            weighted_signal=weighted_signal+(weight*all_signal[i])
        }
    }
    return(weighted_signal)
}

#############################################
## function to compute weighted signal (2016 paper)
## Idea source: https://www.nature.com/articles/nature19362
## We defined a distance scaled signal, pi, to present the ChIP–seq signal intensity of a certain promoter. ρ_i =∑ w_ks_k, p_i is the weighted sum of ChIP–seq reads sk at genomic positions k, where their weights decrease with distance from the TSS of transcription i. In this definition, w_k = 2e^-µ|k-t_i| / 1+e^-u|k-t_i| and t_i is the genomic position of the TSS for transcript i. The parameter μ determines the decay rate, which is a function of the distance from the TSS. For H3K4me3 and H3K27me3 markers, we set μ to be 2 kb from the TSS and the contribution of the corresponding signal decayed to 1/2 of that at the TSS.
## e <- exp(1);
## mu <- 2;
## plot(unlist(lapply(1:50000, function(dist) { dist <- dist/(10*1000); (2*(e^(-mu*dist)))/(1+((e^(-mu*dist)))); })), xlab="dist_to_tss", ylab="weight", ylim=c(0,1))
#############################################
compute_weighted_signal_2016 <- function(signal, distance) {
    mu <- 2
    if(is.na(distance)) {
        weighted_signal <- 0
    } else {
        all_signal <- as.numeric(unlist(strsplit(as.character(signal),",")))
        all_distance <- abs(as.numeric(unlist(strsplit(as.character(distance),","))))/100000

        e <- exp(1)
        weighted_signal=0
        for(i in 1:length(all_distance)) {
            weight=(2*(e^(-mu*all_distance[i])))/(1+(e^(-mu*all_distance[i])))
            weighted_signal=weighted_signal+(weight*all_signal[i])
        }
    }
    return(weighted_signal)
}

#############################################
## function to compute weighted signal (2024 paper)
## Idea source: https://academic.oup.com/nar/article/52/D1/D61/7424438
## For the gene regulator search, assignment of TFs to genes is based on regulatory potential (RP) scores that reflect the collective influence of the binding sites of a given TF on genes nearby these sites and assume that TF binding sites near the TSS are more likely to regulate the gene than those further away. The RP score for gene i and transcription factor j is 2 −x ijk / , where defined as: R ij =  is the decay peak k of TF j distance and x ijk is the genomic distance between peak k of TF j and the TSS of gene i . As different TFs might regulate genes over different ranges of genomic influence, and differ- ent genes can be influenced by enhancers over different ranges ( 34 ), RP scores are calculated for each TF and gene using short (1 kb), mid-range (10 kb) and long-range (100 kb) decay distances.
#############################################
compute_weighted_signal_2024 <- function(signal, distance) {
    if(is.na(distance)) {
        weighted_signal <- 0
    } else {
        all_signal <- as.numeric(unlist(strsplit(as.character(signal),",")))
        all_distance <- abs(as.numeric(unlist(strsplit(as.character(distance),","))))/100000

        e <- exp(1)
        weighted_signal=0
        for(i in 1:length(all_distance)) {
            weight=(2^-all_distance[i])
            weighted_signal=weighted_signal+(weight*all_signal[i])
        }
    }
    return(weighted_signal)
}

#############################################
## compute weighted signal
#############################################
if(opt$weightModel=="2011") {
    WS <-  unlist(lapply(seq(1:length(df[,c(as.numeric(opt$scoreCol))])), function(y) compute_weighted_signal_2011(df[y,c(as.numeric(opt$scoreCol))], df[y,opt$distanceCol])))
} else if(opt$weightModel=="2016") {
    WS <-  unlist(lapply(seq(1:length(df[,c(as.numeric(opt$scoreCol))])), function(y) compute_weighted_signal_2016(df[y,c(as.numeric(opt$scoreCol))], df[y,opt$distanceCol])))
} else {
    WS <-  unlist(lapply(seq(1:length(df[,c(as.numeric(opt$scoreCol))])), function(y) compute_weighted_signal_2024(df[y,c(as.numeric(opt$scoreCol))], df[y,opt$distanceCol])))
}

#############################################
## compute total signal
#############################################
#j=ncol(df)+1
#for(i in unlist(strsplit(opt$scoreCol, ","))) {
#    i <- as.numeric(i)
#    df[,j] <- apply(df, 1, function(x) compute_total_signal(x[i]))
#    j=j+1
#}
TOTAL <- unlist(lapply(df[,c(as.numeric(opt$scoreCol))], function(x) sum(as.numeric(unlist(strsplit(as.character(x),","))))))

#############################################
## compute enhancer count
#############################################
COUNT <- as.data.frame(unlist(lapply(df[,opt$distanceCol], function(x) length(unlist(strsplit(as.character(x), ","))[!is.na(unlist(strsplit(as.character(x), ",")))]))))

#############################################
## organize final output results
#############################################
#print(which(!seq(1:ncol(df)) %in% c(as.numeric(unlist(strsplit(opt$scoreCol, ","))), opt$distanceCol, opt$segmentCol)))
#df_output <- df[,c(which(!seq(1:ncol(df)) %in% c(as.numeric(unlist(strsplit(opt$scoreCol, ","))), opt$distanceCol, opt$segmentCol)))]
df_output <- cbind(df[,1], TOTAL, WS, COUNT)
colnames(df_output) <- c("gene", "total_signal", "weighted_signal", "cre_count")

#############################################
## write final results
#############################################
write.table(df_output, "", sep="\t", row.names=F, col.names=T, quote=F)

#save.session("test.session")
q()
