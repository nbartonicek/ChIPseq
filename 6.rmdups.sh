#!/bin/bash


#module load gi/bwa/0.7.9a
#cd $GENOMES/bwa/mm10_BAC
#bwa index -p mm10_BAC -a bwtsw mm10_BAC.fa 

module load gi/bwa/0.7.8
module load gi/gcc/4.8.2
module load gi/picard-tools/1.138
module load gi/novosort/precompiled/1.03.08

#directory hierarchy
#raw files directory
numcores=12
tag="-P DSGClinicalGenomics"
#tag=""

homedir="/share/ScratchGeneral/nenbar"
projectDir="$homedir/projects/Chris"
resultsDir="$projectDir/project_results/"

genomeFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/bwa/hg38_bwa/hg38_bwa"
#extension of the files to be used
inExt="fastq.gz"


#scripts directory
scriptsPath="/share/ScratchGeneral/nenbar/projects/Chris/scripts"
logDir=$scriptsPath"/logs"

projectnames=( "ELF5" )
i=0
for projectname in ${projectnames[@]}; do


        
        inPath="/share/ScratchGeneral/nenbar/projects/Chris/raw_files/$projectname/renamed/"

	outPath="$projectDir/project_results/$projectname.bwa"
        #log and command files for bsub
        logPath="logs"
        commandPath="commands"
        #make the directory structure   
        mkdir -p $outPath
        mkdir -p $logPath

        subs=0

        #get the name of the script for the logs
        scriptName=`basename $0`
        
        echo $inPath
	files=`ls $inPath/*.fastq.gz`
	for file in ${files[@]};do
                        echo The file used is: $file
                        filesTotal[i]=$file;
                        let i++;
        done 
done;

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

        #dir=`echo ${filesTotal[$j]}`
        #files=`ls $inPath/$dir/*.$inExt`
	
        inFile1=${filesTotal[$j]}
        inFile2=${filesTotal[$(($j+1))]}
        uniqueID=`basename $inFile1 | sed s/_R1.fastq.gz//`
        name=$uniqueID
        outDir=$outPath/
	mkdir -p $outDir
        echo $name
	#echo $command_line


	bwaJobName="dedup."$name
	samSortJobName="samSort"$name
	bamJobName="dedup."$name
	sortJobName="sort."$name
	
	indexJobName="index."$name
	indexStatsJobName="indexstats."$name
	outSam=$outDir"Aligned.out.sam"
	outSortedSam=$outDir"Aligned.sorted.sam"
	outBam=$outDir"$name.bam"
	outSortedBam=$outDir"Aligned.sortedByCoord.out.bam"

        dupDir=$outPath"dedup/"
        mkdir -p $dupDir
	dup_line="samtools fixmate $outDir$name.sorted.bam $dupDir$name.dedup.bam"
        
        # The first sort can be omitted if the file is already name ordered
        dup_line="java -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.138/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$outDir$name.sorted.bam O=$dupDir$name.dedup.bam M=$dupDir$name.metrics.txt"
        dup_line="java -Xmx2g -jar /share/ClusterShare/software/contrib/gi/picard-tools/1.95/CollectInsertSizeMetrics.jar INPUT=$outDir$name.sorted.bam HISTOGRAM_FILE=$dupDir$name.metrix.jpg OUTPUT=$dupDir$name.metrix.txt"
        #echo $dup_line
        qsub -N $bwaJobName -b y -wd $logDir -j y -R y -pe smp $numcores $tag -q short.q -V $dup_line 
	
        j=$(($j+2))


done;
