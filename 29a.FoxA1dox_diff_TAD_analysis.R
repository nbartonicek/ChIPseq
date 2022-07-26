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
tadDir=paste0(resultsDir,"/GSE121443/")

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

#first, how many of them overlap hotspots
data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)))
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
cat(sum(countOverlaps(foxA1DoxUnique,hotGR)>0))

#1857, 1241 that don't overlap ELF5 - possibly other TFs are there. but lets check the distance for 
#hotspots, not hotspots


#TAD analysis 
tadFile<-list.files(tadDir,full.names=T,pattern="GSE121443_enhancer")
tf<-read.table(tadFile,header=T)
coords=strsplit(as.character(tf[,1]),"_")
chromosomes=sapply(coords,function(x){x[1]})
starts=as.numeric(sapply(coords,function(x){x[2]}))
ends=as.numeric(sapply(coords,function(x){x[3]}))

gr<-GRanges(seqnames=chromosomes,IRanges(start=starts,end=ends),seqlengths=seqlengths(cleanGRsPeaks[["FOXA1Dox"]]))
gr<-liftover(gr)


tads=gr

#739 of enhancer tads
#out of those, 


#What percentage of TADS have FOXA1Dox: 42 out of 739, 5.6%
sum(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0)
sum(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0)/length(tads)

#What percentage of TADS have ELF5: 35 out of 739, 4.7%
sum(countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)
sum(countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)/length(tads)


#What percentage have FOXA1Dox and ELF5: 13 of 3841: 1.7%
sum((countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)&(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0))
sum((countOverlaps(tads,cleanGRsPeaks[["ELF5Dox"]])>0)&(countOverlaps(tads,cleanGRsPeaks[["FOXA1Dox"]])>0))/length(tads)

#What percentage of TADS have FOXA1Dox: 2 out of 739, 0.2%
sum(countOverlaps(tads,cleanGRsPeaks[["FOXA1_diff_Dox_enriched"]])>0)
sum(countOverlaps(tads,cleanGRsPeaks[["FOXA1_diff_Dox_enriched"]])>0)/length(tads)
