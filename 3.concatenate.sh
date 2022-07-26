#!/bin/bash

outDir="../ATT1_combined"
dx mkdir -p $outDir

pairs=( "R1" "R2" )
for pair in ${pairs[@]};do 
	for file in `dx ls *L001_$pair*`; do 
		name=`echo $file | sed 's/_L00.*//'`;file2=$name"_L002_"$pair"_001.fastq.gz"; 
		echo $name
		dx run concatenate -ireads1=$file -ireads2=$file2 --yes;
	done;
done


