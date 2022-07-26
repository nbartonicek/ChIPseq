#!/bin/bash

module load gi/repeatmasker/4.0.5
#module load gi/bcftools/1.0
#directory hierarchy
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Chris"
resultsDir="$projectDir/project_results/"
numcores=3
tag="-P DSGClinicalGenomics"
tag=""
#extension of the files to be used
inExt="gz"

#scripts directory
scriptsPath="$homedir/projects/Chris/scripts"
logDir=$scriptsPath"/logs"

inFile="/share/ScratchGeneral/nenbar/projects/Chris/scripts/repeats/ELF5_1kb.fasta"
mkdir -p $logDir
        #out directory

qsub -b y -j y -N repeatMasker -wd $logDir -pe smp $numcores $tag -V "RepeatMasker --species human $inFile -pa $numcores -qq -a -xm -gff -dir /share/ScratchGeneral/nenbar/projects/Chris/project_results/repeats"


