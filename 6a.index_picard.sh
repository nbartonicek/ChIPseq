#!/bin/bash


#directory hierarchy
#raw files directory
numcores=3
tag="-P DSGClinicalGenomics"
tag=""

homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Chris"
resultsDir="$projectDir/project_results/"
logDir="/share/ScratchGeneral/nenbar/projects/Chris/scripts/logs"
inPath=$resultsDir"ELF5.picard"


files=`ls $inPath/*.bam`
for file in ${files[@]};do
        sample=`basename $file`
        qsub -N $sample -b y -wd $logDir -j y -R y -pe smp $numcores $tag -q short.q -V samtools index $file

done;
