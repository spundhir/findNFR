#!/bin/bash

GENOME_SIZE=3000000000

#### usage ####
usage() {
	echo Program: "binomTest.sh (perform binomial test)"
	echo Author: BRIC, University of Copenhagen, Denmark
	echo Version: 1.0
	echo Contact: pundhir@binf.ku.dk
	echo "Usage: binomTest.sh -i <file> -j <file> -k <string>"
	echo "Options:"
    echo " -i <file>   [input file containing features of interest in BED format (can be stdin)]"
    echo " -j <file>   [input file containing reference features in BED format, eg. genome segementation]"
	echo "[optional]"
    echo " -I <string> [list of specific features to filter regions of interest (-i file)]"
    echo " -J <string> [list of specific reference features to search for eg. E, TSS etc (-j file)]"
    echo "             [if multiple, seperate them by a comma]"
    echo " -s <int>    [genome size (default: 3000000000)]"
	echo " -h          [help]"
	echo
	exit 0
}

MAPPING_FREQUENCY=1

#### parse options ####
while getopts i:j:I:J:s:h ARG; do
	case "$ARG" in
        i) FILE_INT=$OPTARG;;
        j) FILE_REF=$OPTARG;;
        I) FILTER_INT=$OPTARG;;
        J) FILTER_REF=$OPTARG;;
        s) GENOME_SIZE=$OPTARG;;
		h) HELP=1;;
	esac
done

## usage, if necessary file and directories are given/exist
if [ -z "$FILE_INT" -o ! -f "$FILE_REF" -o "$HELP" ]; then
	usage
fi

## create temporary BED file if input is from stdin
TMP=$(date | md5sum | cut -f 1 -d " ")
#TMP=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 32 | head -n 1)
TMP="8297ab2cc40b06008ace7b732b89b1ff"
if [ "$FILE_INT" == "stdin" ]; then
    while read LINE; do
        echo -e "${LINE}"
    done > $TMP
    FILE_INT=$TMP
else
    scp $FILE_INT $TMP
    FILE_INT=$TMP
fi

if [ ! -z "$FILTER_INT" ]; then
    perl -ane 'if($_=~/\s+'$FILTER_INT'\s+/) { print $_; }' $FILE_INT > $FILE_INT.tmp;
    mv $FILE_INT.tmp $FILE_INT
fi

echo -e "#file\tfeatures\ttotal\toverlap\tmean\tstddev\tp-value\tpercentage\texp_overlap\texp_per"
#echo -e "$FILE_REF\t$FILE_INT\t$FILTER_REF"; exit
if [ -z "$FILTER_REF" ]; then
    N=`zless $FILE_INT | wc -l | cut -f 1 -d " "`;
    Pr=`tabEdit -i $FILE_REF -D | perl -ane '$cov+=($F[2]-$F[1])+1; END { printf("%0.6f", $cov/'$GENOME_SIZE'); }'`;
    mean=`echo $Pr | perl -ane '$mean='$N'*$_; printf("%0.6f", $mean);'`;
    stdev=`echo $Pr | perl -ane '$stdev=sqrt('$N'*$_*(1-$_)); printf("%0.6f", $stdev);'`;

    overlap=`intersectBed -a $FILE_INT -b <(tabEdit -i $FILE_REF -D) -u | wc -l`;
    pvalue=`Rscript ${FINDNFRPATH}/share/R/pnorm.R $overlap $mean $stdev | cut -f 2 -d " "`;
    per=`perl -e '$per=('$overlap'*100)/'$N'; printf("%0.2f", $per);'`;

    exp_overlap=`Rscript ${FINDNFRPATH}/share/R/qnorm.R $mean $stdev | cut -f 2 -d " "`;
    exp_per=`perl -e '$exp_per=('$exp_overlap'*100)/'$N'; printf("%0.2f", $exp_per);'`;

    #echo -e "$entity\t$N\t$Pr\t$mean\t$stdev\t$overlap\t$per\t$exp_overlap\t$exp_per"; exit;
    file=`echo $FILE_REF | sed 's/^.*\///g'`;
    echo -e "$file\tNA\t$N\t$overlap\t$mean\t$stdev\t$pvalue\t$per\t$exp_overlap\t$exp_per";
else
    IFS=","
    FEATURES=($FILTER_REF)
    IFS=""
    FIELD=$(intersectBed -a $FILE_INT -b <(tabEdit -i $FILE_REF -D) -wo | head -n 1 | perl -ane '$field=scalar(@F); printf("%drn,%d", $field, $field);')
    NCOL=$(intersectBed -a $FILE_INT -b <(tabEdit -i $FILE_REF -D) -wo | head -n 1 | perl -ane 'print scalar(@F);')
    #echo -e "$FIELD\t$NCOL"; exit

    for (( i=0; i<${#FEATURES[@]}; i++ )); do
        entity=${FEATURES[$i]};
        ENTITY_COL=$(intersectBed -a $FILE_INT -b <(tabEdit -i $FILE_REF -D) -wo | grep -w $entity | head -n 1 | perl -ane 'BEGIN { $col=1; } foreach(@F) { if($_=~/^'$entity'$/) { print $col; last; } $col++; }')
        N=`zless $FILE_INT | wc -l | cut -f 1 -d " "`;
        Pr=`zgrep -w $entity <(tabEdit -i $FILE_REF -D) | perl -ane '$cov+=($F[2]-$F[1])+1; END { printf("%0.6f", $cov/'$GENOME_SIZE'); }'`;
        mean=`echo $Pr | perl -ane '$mean='$N'*$_; printf("%0.6f", $mean);'`;
        stdev=`echo $Pr | perl -ane '$stdev=sqrt('$N'*$_*(1-$_)); printf("%0.6f", $stdev);'`;

        if [ ! -z "$ENTITY_COL" ]; then
            #echo "intersectBed -a $FILE_INT -b <(tabEdit -i $FILE_REF -D) -wao | sort -k 1,1 -k 2n,2 -k 3n,3 -k $FIELD | perl -ane '\$key=\"\$F[0]_\$F[1]_\$F[2]\"; if(!defined(\$seen{\$key})) { print \$_; \$seen{\$key}=1;}' | cut -f $ENTITY_COL | sort | uniq -c | sed -E 's/^\s+//g' | grep -w $entity | cut -f 1 -d \" \""; exit;
            overlap=`intersectBed -a $FILE_INT -b <(tabEdit -i $FILE_REF -D) -wao | sort -k 1,1 -k 2n,2 -k 3n,3 -k $FIELD | perl -ane '$key="$F[0]_$F[1]_$F[2]"; if(!defined($seen{$key})) { print $_; $seen{$key}=1;}' | cut -f $ENTITY_COL | sort | uniq -c | sed -E 's/^\s+//g' | perl -ane 'if($F[1]=~/^'$entity'$/) { print $_; }' | cut -f 1 -d " "`;
            if [ -z "$overlap" ]; then
                overlap=0
            fi
        else 
            overlap=0
        fi
        #echo -e "Entity: $entity; N: $N; Pr: $Pr; Mean: $mean; Stdev: $stdev; Overlap: $overlap; Expected overlap: $exp_overlap"; exit;

        pvalue=`Rscript ${FINDNFRPATH}/share/R/pnorm.R $overlap $mean $stdev | cut -f 2 -d " "`;
        per=`perl -e '$per=('$overlap'*100)/'$N'; printf("%0.2f", $per);'`;

        exp_overlap=`Rscript ${FINDNFRPATH}/share/R/qnorm.R $mean $stdev | cut -f 2 -d " "`;
        exp_per=`perl -e '$exp_per=('$exp_overlap'*100)/'$N'; printf("%0.2f", $exp_per);'`;
        
        #echo -e "$entity\t$N\t$Pr\t$mean\t$stdev\t$overlap\t$pvalue\t$exp_overlap";
        file=`echo $FILE_REF | sed 's/^.*\///g'`;
        echo -e "$file\t${FEATURES[$i]}\t$N\t$overlap\t$mean\t$stdev\t$pvalue\t$per\t$exp_overlap\t$exp_per";
    done

    ## performing enrichment analysis for features that do not overlap with reference
    N=`zless $FILE_INT | wc -l | cut -f 1 -d " "`;
    Pr=`zless <(tabEdit -i $FILE_REF -D) | perl -ane '$cov+=($F[2]-$F[1])+1; END { printf("%0.6f", ('$GENOME_SIZE'-$cov)/'$GENOME_SIZE'); }'`;
    mean=`echo $Pr | perl -ane '$mean='$N'*$_; printf("%0.6f", $mean);'`;
    stdev=`echo $Pr | perl -ane '$stdev=sqrt('$N'*$_*(1-$_)); printf("%0.6f", $stdev);'`;
    overlap=`intersectBed -a $FILE_INT -b <(tabEdit -i $FILE_REF -D) -v | wc -l`;
    pvalue=`Rscript ${FINDNFRPATH}/share/R/pnorm.R $overlap $mean $stdev | cut -f 2 -d " "`;
    per=`perl -e '$per=('$overlap'*100)/'$N'; printf("%0.2f", $per);'`;

    exp_overlap=`Rscript ${FINDNFRPATH}/share/R/qnorm.R $mean $stdev | cut -f 2 -d " "`;
    exp_per=`perl -e '$exp_per=('$exp_overlap'*100)/'$N'; printf("%0.2f", $exp_per);'`;

    #echo -e "$entity\t$N\t$Pr\t$mean\t$stdev\t$overlap\t$pvalue\t$exp_overlap";
    file=`echo $FILE_REF | sed 's/^.*\///g'`;
    echo -e "$file\tother\t$N\t$overlap\t$mean\t$stdev\t$pvalue\t$per\t$exp_overlap\t$exp_per";
fi

if [ ! -z "$TMP" ]; then
    rm $TMP
fi

## regioneR package
#library(regioneR)
#library("BSgenome.Mmusculus.UCSC.mm9")
#peaks_cebpe <- toGRanges(pipe("grep -w promoter /home/pundhir/project/chip-seq-analysis/analysis_cebpe/analysis/peaks_cebpe/peaks_cebpe_class.bed"), header=F)
#peaks_e2f <- toGRanges(pipe("grep -v chr8_random /home/pundhir/project/chip-seq-analysis/analysis_cebpe/analysis/E2f1_mm9_lifted.bed"))

#cebpe_e2f <- permTest(A=peaks_e2f, B=peaks_cebpe, ntimes=1000,
#                       randomize.function=circularRandomizeRegions,
#                       evaluate.function=numOverlaps, count.once=TRUE,
#                       genome="mm9", mc.set.seed=FALSE, mc.cores=1)

#plot(cebpe_e2f)
