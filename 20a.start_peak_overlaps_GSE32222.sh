#!/bin/bash

module load briglo/R/3.4.2
#module load gi/bcftools/1.0
#directory hierarchy
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Chris"
resultsDir="$projectDir/project_results/"
numcores=8
tag="-P DSGClinicalGenomics"
#extension of the files to be used
inExt="gz"

#scripts directory
scriptsPath="$homedir/projects/Chris/scripts"
logDir=$scriptsPath"/logs"

mkdir -p $logDir
        #out directory
types=( "resistant" "sensitive" )
for i in {1..11}; do
	for type in ${types[@]}; do
        qsub -b y -j y -N pt$i$type -wd $logDir -pe smp $numcores $tag -V "R --vanilla <$scriptsPath/20a.peak_overlaps_GSE75372.R --id=$i --type=$type"
	done;
done;

