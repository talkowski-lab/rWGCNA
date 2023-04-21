#!/bin/bash

###############################################
# This script is used to submit rWGCNA_consensus.R to the server
# by Heather Wick hwick@broadinstitute.org heather.c.wick@gmail.com
# usage:
# bash submit_rWGCNA_args_blockwise.sh project_name [num_iter] /path/to/Data_counts_combat_norm_log2.txt /path/to/output/base/dir/ /path/to/traitData.txt [block_size] [bool_use_default_soft_power_threshold]
# example:
# bash submit_rWGCNA_args_blockwise.sh project_name 100 /path/to/Data_counts_combat_norm_log2.txt /path/to/output/base/dir/ /path/to/traitData.txt 17000 F
###############################################

module purge
module load gcc/6.3.0
module load glibc/2.14
module load R/4.0.2

name=$1 #project name
iter=$2 #number of times to resample/bootstrap data set
data=$3 #log2 normalized expression data
baseDir=$4 #base directory to contain the output
traitData=$5 #trait data with Samples column corresponding to colnames of data
maxBlockSize=$6 #I recommend using a number slightly larger than the number of genes, if you have the computational power
defaultPower=$7 #should be T or F; F will calculate soft power threshold for each resampled data set and choose the median
WORKDIR=${name}_NumIter${iter}
PROJECTDIR=${baseDir}${WORKDIR}
mkdir -p ${PROJECTDIR}/Plots ${PROJECTDIR}/Jobs
echo "PD: "$PROJECTDIR

bsub -sla miket_sc -q big-multi -n 4 -M 20000 -v 20000 -R 'rusage[mem=20000]' \
   -o ${PROJECTDIR}/Jobs/%J.out -e ${PROJECTDIR}/Jobs/%J.err \
   "Rscript --vanilla --verbose rWGCNA_args_blockwise.R ${PROJECTDIR} $name $iter $data $traitData $maxBlockSize $defaultPower"
