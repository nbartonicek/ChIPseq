#!/bin/bash

module load briglo/R/3.4.2
#module load gi/bcftools/1.0
#directory hierarchy
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Chris"
resultsDir="$projectDir/project_results/"
numcores=3
#tag="-P DSGClinicalGenomics"
tag=""
#extension of the files to be used
inExt="gz"

#scripts directory
scriptsPath="$homedir/projects/Chris/scripts"
logDir=$scriptsPath"/logs"

mkdir -p $logDir
        #out directory
classes=( "enhancer" "superenhancer" "DE" "resistant" "sensitive" )
classes=( "DE" )
for typeID in {67..134}; do
	for class in ${classes[@]}; do
        qsub -q short.q -b y -j y -N pt$typeID$class -wd $logDir -pe smp $numcores $tag -V "R --vanilla <$scriptsPath/30b.grammar_pvalues_parallel.R --typeID=$typeID --class=$class"
	done;
done;

