#!/bin/bash

gi/homer/4.7

#directory hierarchy
#raw files directory
numcores=8
tag="-P DSGClinicalGenomics"
tag=""

homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Chris"
resultsDir="$projectDir/project_results/"

genomeFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/bwa/hg38_bwa/hg38_bwa"
#extension of the files to be used
inExt=".bam"


#scripts directory
scriptsPath="/share/ScratchGeneral/nenbar/projects/Chris/scripts"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

i=0
        
inPath="/share/ScratchGeneral/nenbar/projects/Chris/project_results/ELF5.picard/"
inputInPath="/share/ScratchGeneral/nenbar/projects/Chris/project_results/ELF5.picard.merged/"
outPath="/share/ScratchGeneral/nenbar/projects/Chris/project_results/ELF5.homer_separate/"
#log and command files for bsub
commandPath="commands"
#make the directory structure   
mkdir -p $outPath

echo $inPath
files=`ls $inPath/*.bam`

tfs=( "ELF5" "ER" "FOXA1" )
#tfs=( "ELF5" )
treatments=( "Dox" "NoDox")
#treatments=( "Dox" )

for tf in ${tfs[@]};do
        for treatment in ${treatments[@]}; do
                for i in {1..4};do
                       homerJobName=$tf$treatment$i
                       mkdir -p $outPath$homerJobName
                       homerDir=$outPath$homerJobName
                       cp 
                       inputDir=
                        macsLine="findPeaks $homerDir \
                                -style factor \
                                
                                -i 
                                -c "$inPath"Pool1Input$treatment""2.dedup.bam "$inPath"Pool1Input$treatment""4.dedup.bam "$inPath"Pool2Input$treatment""2.dedup.bam "$inPath"Pool2Input$treatment""4.dedup.bam "$inPath"Pool3Input$treatment""2.dedup.bam "$inPath"Pool3Input$treatment""4.dedup.bam \
                                --keep-dup all -f BAMPE -n $tf$treatment$i -B -q 0.01 --outdir $outPath"
                        #echo $macsLine
                        qsub -N $macsJobName -b y -wd $logDir -j y -R y -pe smp $numcores $tag -q short.q -V $macsLine 
                done;
        done;
done;

#tfs=( "H3K4me3" )
#treatments=( "Dox" )
#
#
#for tf in ${tfs[@]};do
#        for treatment in ${treatments[@]}; do
#                macsJobName=$tf$treatment
#                is=( "2" "4" )
#                for i in ${is[@]};do
#                
#                        macsLine="macs2 callpeak \
#                                -t $inPath$tf$treatment""$i.dedup.bam  \
#                                -c "$inPath"Pool1Input$treatment""2.dedup.bam "$inPath"Pool1Input$treatment""4.dedup.bam "$inPath"Pool2Input$treatment""2.dedup.bam "$inPath"Pool2Input$treatment""4.dedup.bam "$inPath"Pool3Input$treatment""2.dedup.bam "$inPath"Pool3Input$treatment""4.dedup.bam \
#                                --keep-dup all --broad -f BAMPE -g hs -n $tf$treatment$i -B -q 0.01 --outdir $outPath"
#        #                #echo $macsLine
#                        qsub -N $macsJobName$i -b y -wd $logDir -j y -R y -pe smp $numcores $tag -q short.q -V $macsLine 
#                done;
#        done;
#done;
#
#
#tfs=( "H3K4me3" )
#treatments=( "NoDox" )
#
#
#for tf in ${tfs[@]};do
#        for treatment in ${treatments[@]}; do
#                macsJobName=$tf$treatment
#                is=( "1" "3" )
#                for i in ${is[@]};do
#                        macsLine="macs2 callpeak \
#                                -t $inPath$tf$treatment""$i.dedup.bam $inPath$tf$treatment""3.dedup.bam \
#                                -c "$inPath"Pool1Input$treatment""2.dedup.bam "$inPath"Pool1Input$treatment""4.dedup.bam "$inPath"Pool2Input$treatment""2.dedup.bam "$inPath"Pool2Input$treatment""4.dedup.bam "$inPath"Pool3Input$treatment""2.dedup.bam "$inPath"Pool3Input$treatment""4.dedup.bam \
#                                --keep-dup all --broad -f BAMPE -g hs -n $tf$treatment$i -B -q 0.01 --outdir $outPath"
#                        #echo $macsLine
#                        qsub -N $macsJobName$i -b y -wd $logDir -j y -R y -pe smp $numcores $tag -q short.q -V $macsLine 
#                done;
#        done;
#done;
#
#
#
#