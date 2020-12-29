#!/bin/bash
#
# This script will find the consensus peak regions from peak files (in
# BED format) of multiple samples by:
#
# 1. Converting the peak file of each sample into non-overlapping 3
# cols BED file and concatenating them;
#
# 2. Sorting the concatenated file and Building a genome coverage
# track in BedGraph, of which the value (the 3rd col) indicates the
# number of samples with a peak at a particular location;
# 
# 3. Using MACS to call regions being covered by peaks from more than
# a certain number of samples.

# ------------------------------------------------------------------

if [ $# -lt 2 ];then
    echo "Need at least 3 parameters! <consensus prefix> <cutoff> <all the peak files in BED format>"
    exit
fi

OUTPUT_PREFIX=$1
CUTOFF=$2
SAMPLE_PEAKS=${@:3}

# Modify the following parameters
# define the minlen and maxgap for peak calling (arbitrary)
MINLEN=200
MAXGAP=30

# ------------------------------------------------------------------

# 1 pileup 
A=${OUTPUT_PREFIX}.all.bed
P=${OUTPUT_PREFIX}.pileup.bdg

rm -f $A $P
touch $A $P

cat $SAMPLE_PEAKS | cut -f 1,2,3 > $A
macs3 pileup -i $A -f BEDPE -o $P

#2 
O=${OUTPUT_PREFIX}.consensus.bed

macs3 bdgpeakcall -i $P -o $O --no-trackline -c $CUTOFF -g $MAXGAP -l $MINLEN

# end
echo "All done. Check ${O}"

