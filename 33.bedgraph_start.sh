#!/bin/bash

module load gi/igvtools/2.3.23
#module load gi/bcftools/1.0
#directory hierarchy
#raw files directory
homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Chris"
resultsDir="$projectDir/project_results/"
numcores=6
#tag="-P DSGClinicalGenomics"
tag=""
#extension of the files to be used
inExt="_pileup.bedgraph"

#scripts directory
scriptsPath="$homedir/projects/Chris/scripts"
logDir=$scriptsPath"/logs"
inDir=$resultsDir"ELF5.macs"
mkdir -p $logDir
        #out directory
for file in `ls $inDir/*$inExt | grep H3K4me3`;do
	sampleName=`basename $file`
	inFile=`echo $sampleName | sed s/bdg/bedgraph/`
	outFile=`echo $sampleName | sed s/_treat_pileup.bdg/.tdf/`

	#mv $file $inDir/$inFile
	

    qsub -q short.q -b y -j y -N inFile -wd $logDir -pe smp $numcores $tag -V "igvtools toTDF $inDir/$inFile $inDir/$outFile $inDir/hg38.chrom.sizes"
done;

