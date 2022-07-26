#!/bin/bash
module load briglo/R/3.4.2

#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="/share/ScratchGeneral/nenbar/"
projectDir=$homedir"projects/Chris/"
resultsDir=$projectDir"project_results/"


numcores=8

#extension of the files to be used
inExt="dedup.bam$"

#scripts directory
scriptsPath=$projectDir"/scripts"
logDir=$scriptsPath"/logs"
mkdir -p $logPath

projectnames=( "ELF5" )

for projectname in ${projectnames[@]}; do


	inPath="$resultsDir$projectname.picard/"
        outDir=$resultsDir$projectname".bigwig"
        mkdir -p $outPath

        echo $inPath
	files=`ls $inPath | grep $inExt | grep -v NoDox | grep Pool`
	for file in ${files[@]};do
                uniqueID=`basename $file | sed s/.dedup.bam//`
                echo $uniqueID
                bwJobName="bw."$uniqueID
                bw_line="R --vanilla <$scriptsPath/11.bam2bw.R --inFile=$inPath$file"
                qsub -N $uniqueID -b y -wd $logDir -j y -R y -pe smp $numcores -q short.q -V $bw_line

        done 
done;



