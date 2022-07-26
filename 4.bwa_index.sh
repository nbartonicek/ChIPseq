#!/bin/bash

#load all modules
#number of cores
ncore=16

#project directory
#homedir="/share/ClusterScratch/nenbar"

#genome directory
genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/bwa/hg38_ercc"
mkdir -p $genomeDir

bwa index -a bwtsw hg38_ercc.fa -p hg38_bwa &


