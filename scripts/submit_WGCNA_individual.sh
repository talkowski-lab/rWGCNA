#!/bin/bash

###############################################
# This script is used to submit WGCNA_individual.R to the server
# by Heather Wick hwick@broad institute.org
# usage:
# bash submit_rWGCNA_args_blockwise.sh project_name /path/to/Data_counts_combat_norm_log2.txt /path/to/output/base/dir/ /path/to/traitData.txt [block_size]
# example:
# bash submit_rWGCNA_args_blockwise.sh project_name /path/to/Data_counts_combat_norm_log2.txt /path/to/output/base/dir/ /path/to/traitData.txt 17000
###############################################

module purge
module load gcc/6.3.0
module load glibc/2.14
module load R/4.0.2

name=$1 #project name
data=$2 #log2 normalized expression data
baseDir=$3 #base directory to contain the output
traitData=$4 #temporary arg while testing
maxBlockSize=$5 #I recommend using a number slightly larger than the number of genes, if you have the computational power
WORKDIR=${name}_bs${iter}_${test_info}
PROJECTDIR=${baseDir}${WORKDIR}
mkdir -p ${baseDir}${WORKDIR}/Plots ${baseDir}${WORKDIR}/Jobs
echo "WD: "$WORKDIR
echo "PD: "$PROJECTDIR


bsub -sla miket_sc -q big-multi -n 4 -M 20000 -v 20000 -R 'rusage[mem=20000]' \
      -o ${PROJECTDIR}/Jobs/%J.out -e ${PROJECTDIR}/Jobs/%J.err \
      "Rscript --vanilla --verbose WGCNA_individual.R ${PROJECTDIR} $name $data $test_info $maxBlockSize"


