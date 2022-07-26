#data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66733

library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2018)
library(Biostrings)
library(ChIPpeakAnno)
library(UpSetR) 
library(grid)
library(plyr)
library(gplots)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(org.Hs.eg.db)
library(circlize)

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
chrGR<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))
#homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs/")
inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate/")

inBams=paste0(homedir,"/projects/Chris/project_results/ELF5.picard/")
outPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate_consensus/")
system(paste0("mkdir -p ",outPath))



######## directory structure #######
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/TADdifferent/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
tadDir=paste0(resultsDir,"/GSE66733/Hi-C_MCF7_MCF10A_processed_HiCfiles/TAD_boundaries/")

system(paste0("mkdir -p ",imageDir))

liftover<-function(gr){
  gr.hg19=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg19ToHg38.over.chain"))
  seqlevelsStyle(gr.hg19) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg19, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
load(paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))

for(file in list.files(peakRobjectsDir,pattern="Rdata",full.names=T)){
  load(file)
  sampleName<-gsub(".Rdata","",basename(file))
  cat(sampleName)
  cat("\n")
  gr=experiment.DB
  values(gr)=NULL
  gr$score=experiment.DB$Fold
  gr$pval=-log10(experiment.DB$FDR)
  grPos<-gr[gr$score>0]
  grNeg<-gr[gr$score<0]
  cleanGRsPeaks[[paste0(sampleName,"_Dox_enriched")]]<-grPos
  cleanGRsPeaks[[paste0(sampleName,"_Dox_depleted")]]<-grNeg

}

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


allPeaks<-reduce(unlist(cleanGRsPeaks[1:5]))
for(sampleName in names(cleanGRsPeaks[1:5])){
  cat(sampleName)
  cat("\n")
  values(allPeaks)[sampleName]<-countOverlaps(allPeaks,cleanGRsPeaks[sampleName])
}

df<-as.data.frame(values(allPeaks))

#find the FOXA1 clean ones
uniqueFD<-df[(df$FOXA1Dox==1)&(df$FOXA1NoDox!=1),]
foxA1DoxUnique<-allPeaks[(df$FOXA1Dox==1)&(df$FOXA1NoDox!=1)]
mat<-findOverlaps(foxA1DoxUnique,cleanGRsPeaks[["ELF5Dox"]])
foxA1DoxUnique=foxA1DoxUnique[-unique(queryHits(mat))]

#first, how many of them overlap hotspots
data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)))
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
cat(sum(countOverlaps(foxA1DoxUnique,hotGR)>0))

#1857, 1241 that don't overlap ELF5 - possibly other TFs are there. but lets check the distance for 
#hotspots, not hotspots


#TAD analysis 
tadFiles<-list.files(tadDir,full.names=T,pattern="MCF7")
tadL<-GRangesList()
for(tadFile in tadFiles){
  sampleName=gsub("HiCStein-MCF7-WT__hg19__","",basename(tadFile))
  sampleName=gsub("_.*","",sampleName)
  cat(sampleName)
  cat("\n")
  tf<-fread(tadFile)
  tf=as.data.frame(tf)
  gr<-GRanges(seqnames=sampleName,IRanges(start=as.integer(apply(tf[,c("start","end")],1,mean)),width=1),seqlengths=seqlengths(cleanGRsPeaks[["FOXA1Dox"]]))
  gr<-liftover(gr)

  boundaries<-setdiff(chrGR,gr)
  tadL[[sampleName]]=boundaries
}

tads<-unlist(tadL)

#Ok, now that we have TADS

#What percentage of TADS have FOXA1Dox: 3699 out of 3841, 96%
sum(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0)
sum(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0)/length(tads)

#What percentage of TADS have ELF5: 3385 out of 3841, 88%
sum(countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)
sum(countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)/length(tads)


#What percentage have FOXA1Dox and ELF5: 3353 of 3841: 87%
sum((countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)&(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0))
sum((countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)&(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0))/length(tads)

#What percentage of TADS have FOXA1Dox: 2561 out of 3841, 67%
sum(countOverlaps(tads,foxA1DoxUnique)>0)
sum(countOverlaps(tads,foxA1DoxUnique)>0)/length(tads)


#What percentage of FOXA1UniqueDox without ELF5 have ELF5 in TAD: 2464 of 2561, 96%
tadUniq=tads[countOverlaps(tads,foxA1DoxUnique)>0]
sum(countOverlaps(tadUniq,cleanGRsPeaks[["ELF5Dox"]])>0)
sum(countOverlaps(tadUniq,cleanGRsPeaks[["ELF5Dox"]])>0)/length(tadUniq)


#check with a different Tad
tadFile2="/share/ScratchGeneral/nenbar/projects/Chris/project_results/GSE68858/GSM1684569_MCF7_intra_loci_pairs.txt"

tad2<-read.table(tadFile2,header=F)
gr1<-GRanges(seqnames=paste0("chr",tad2$V1),IRanges(start=tad2$V3,end=tad2$V4))
gr2<-GRanges(seqnames=paste0("chr",tad2$V2),IRanges(start=tad2$V5,end=tad2$V6))
gr1<-liftover(gr1)
gr2<-liftover(gr2)
#find those FOXA1 that overlap gr1, are there more than expected to have ELF5
mat1<-findOverlaps(resize(foxA1DoxUnique,10000,fix="center"),gr1)
mat2<-findOverlaps(resize(cleanGRsPeaks[["ELF5Dox"]],10000,fix="center"),gr2)

#out of 4601 peaks, 773 (or 17% have an interaction read)
foxa1Partner=gr2[unique(subjectHits(mat1))]
#how many of those overlap ELF5?
#ELF5 overlaps: 14% overlap the second region
#the second region has ELF5 peaks in 5636 cases (4.5%)
#cat(length(cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat2))]))
cat(length(gr2[unique(subjectHits(mat2))])/length(gr2))

#elf5OverlapOnGr2=gr2[unique(subjectHits(mat2))]
mat3<-findOverlaps(resize(cleanGRsPeaks[["ELF5Dox"]],10000,fix="center"),foxa1Partner)
#cat(length(cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat3))]))
cat(length(foxa1Partner[unique(subjectHits(mat3))])/length(foxa1Partner))


#find those FOXA1 that overlap gr1, are there more than expected to have ELF5
mat1<-findOverlaps(resize(foxA1DoxUnique,10000,fix="center"),gr2)
mat2<-findOverlaps(resize(cleanGRsPeaks[["ELF5Dox"]],10000,fix="center"),gr1)

#out of 4601 peaks, 773 (or 17% have an interaction read)
foxa1Partner=gr1[unique(subjectHits(mat1))]
#how many of those overlap ELF5?
#ELF5 overlaps: 14% overlap the second region
#the second region has ELF5 peaks in 5636 cases (4.5%)
#cat(length(cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat2))]))
cat(length(gr1[unique(subjectHits(mat2))])/length(gr1))

#elf5OverlapOnGr2=gr2[unique(subjectHits(mat2))]
mat3<-findOverlaps(resize(cleanGRsPeaks[["ELF5Dox"]],10000,fix="center"),foxa1Partner)
#cat(length(cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat3))]))
cat(length(foxa1Partner[unique(subjectHits(mat3))])/length(foxa1Partner))

























































