
====================================================================
README file for findNFR (v0.01)
====================================================================

Description
===========
    findNFR is a computational method to predict Nucleosome Free Regions corresponding to input genomic regions in BED format.

Version
=======
    0.01

Citation
========

Programs and datasets
=====================

Installation
============

    1. To install findNFR, download findNFR.tar.gz and unpack it. A directory, findNFRSuite will be created

        tar -zxvf findNFR.tar.gz

    2. Now compile and create executable blockbuster

        make or make all

    3. Export environment variable 'FINDNFRPATH' containing path to findNFR installation directory

        export FINDFRPATH=<path to findNFR installation directory>

    4. Add 'FINDNFRPATH' to your 'PATH' environment variable

        export PATH=$PATH:$FINDNFRPATH/bin

    5. Add 'FINDNFRPATH' to your 'PERL5LIB' environment variable

        export PERL5LIB=$PERL5LIB:$FINDNFRPATH/share/perl/

    To permanently add or update the environment variable(s), add the last three export commands in your ~/.bashrc file

Dependency
==========

    We assume that the following programming platforms are installed and working: perl, R, and gcc. Besides, following packages should be installed.

    1. Install the needed perl modules

        sudo cpan Tie::IxHash Statistics::Basic

    2. R modules are installed by entering R (type R on the cmdline) and then enter the following three commands (follow the instructions on the screen):

        install.packages(c("ggplot2", "gridExtra", "optparse", "randomForest", "e1071"))

        source("http://bioconductor.org/biocLite.R")

        biocLite(c("DESeq"))

    3. download samtools from http://sourceforge.net/projects/samtools/files/samtools/1.2/samtools-1.2.tar.bz2/download, go to the download location and do

        tar xjf samtools-1.2.tar.bz2

        cd samtools-1.2

        make -j10 prefix=$HOME install

    4. download bedtools from https://github.com/arq5x/bedtools2/releases/download/v2.23.0/bedtools-2.23.0.tar.gz, go to the download location and do

        tar xzf BEDTools.v2.23

        cd bedtools-2.23.0/

        make -j 10

        cp bin/* $HOME/bin

    5. download featureCounts (subread) from http://sourceforge.net/projects/subread/files/subread-1.4.6-p4/, go to the download location and do

        tar xzf subread-1.4.6-p4-Linux-x86_64.tar.gz
        
        cd subread-1.4.6-p3-Linux-x86_64
        
        cp bin/featureCounts $HOME/bin

    6. download bedGraphToBigWig from http://hgdownload.soe.ucsc.edu/admin/exe/ for your operating system, go to the download location and do

        cp bedGraphToBigWig $HOME/bin

        chmod 755 $HOME/bin/bedGraphToBigWig

    7. download macs2 version 2.1.0 from https://github.com/taoliu/MACS/, go to the download location and install as mentioned in INSTALL.rst file

Usage
=====

    findNFR is called with the following parameters

    pare -i <BAM file(s)> -k <BED file> [OPTIONS]

Example
=======

    findNFR -i data/histone_Rep1.bam,data/histone_Rep2.bam -k input.bed -o results -m hg19 -p 10 &>findNFR.log

Input
=====

    - BAM file(s) corresponding to each replicate(s)

    - BED file containing genomic coordinates of regions of interest

    The chromosome identifier in the input BAM files should start with chr, for example as chrY and not like Y.

Output
======

    The results from the findNFR are presented in two text files:

    a) RESULTS.TXT: main result file in BED format 

    For easy access, the html version of this file (RESULTS.HTML) is also available within the output directory

    b) RESULTS.UCSC: file to view the enhancer and promoter regions in UCSC browser

More info
=========

License
=======

    findNFR: a computational method to Predict Nucleosome Free Regions using histone marks

    Copyright (C) 2017  Sachin Pundhir (pundhir@binf.ku.dk)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

