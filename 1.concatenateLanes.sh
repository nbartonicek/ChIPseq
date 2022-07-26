#!/bin/bash

sample="L5"
tag="-P DSGClinicalGenomics"
homeDir="/share/ScratchGeneral/nenbar/projects/Chris"
scriptDir="/share/ScratchGeneral/nenbar/projects/Chris/scripts/"
logDir=$scriptDir"logs"
mkdir -p $logDir
numcoresSmall=1
queue="short"


qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 1 $tag -V"

declare -A used
outdir="$homeDir/raw_files/$sample/merged"
mkdir -p $outdir
for file in `ls $homeDir/raw_files/$sample/*.fastq.gz`;do
#for file in `dx ls /projects/$sample/trimgalore/$sample/*.fq.gz`;do 
        SUBMIT=()
        name=`echo $file | sed 's/.*_180522_CIX--/180522_CIX--/' | sed 's/_Human_.*M00/_M00/'`
        #echo $file
        echo $name
        SUBMIT+=($file)
        for secondFile in `ls $homeDir/raw_files/$sample/*.fastq.gz`;do
                secondName=`echo $secondFile | sed 's/.*_180522_CIX--/180522_CIX--/' | sed 's/_Human_.*M00/_M00/'`
                if [[ "$file" != "$secondFile" ]] && [[ $name == $secondName ]]; then
                        SUBMIT+=($secondFile)
                        #echo $secondFile
                fi;
        done;

        #run the app
        if [[ ${used[$name]} != 1 ]] ; then
                $qsubLine -N concat_$name "cat  ${SUBMIT[0]} ${SUBMIT[1]} > $outdir/$name;"
                #echo ${SUBMIT[0]}
                #echo ${SUBMIT[1]}              
        fi
        used[$name]=1

        #dx run MRTA:/concatenate_lanes --dest /raw_files/$sample -ireads1=/raw_files/$sample/150626_D00119_0128_AC75H9ANXX/lanes/C75H9ANXX_7_$name -ireads2=/raw_files/$sample/150626_D00119_0128_AC75H9ANXX/lanes/C75H9ANXX_8_$name -isample_name=$name --yes;
done;




