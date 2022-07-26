#enhancers were downloaded from this paper: https://www.nature.com/articles/s41598-017-02257-3#Sec26


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

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(org.Hs.eg.db)
library(circlize)

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]

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


load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["FOXA1Dox"]]
mat<-findOverlaps(tssDE,foxA1DoxUnique)
tssDE<-tssDE[unique(queryHits(mat))]

tssDE$hotspot<-"noHotspot"
tssDE$hotspot[countOverlaps(tssDE,hotGR)>0]="hotspot"

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)


pdf(paste0(imageDir,"dist_FOXA1DoxUnique_ELF5_ER_2kb.pdf"),width=6,height=6)
col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1, use_raster = TRUE,col = col_fun, column_title = "ELF5Dox",name ="ELF5Dox")+
EnrichedHeatmap(mat2, use_raster = TRUE,col = col_fun, column_title = "FOXA1Dox",name ="FOXA1Dox")+
EnrichedHeatmap(mat3, use_raster = TRUE,col = col_fun, column_title = "ERDox",name ="ERDox")
dev.off()


load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["FOXA1Dox"]]
mat<-findOverlaps(tssDE,foxA1DoxUnique)
tssDE<-tssDE[unique(queryHits(mat))]

tssDE$hotspot<-"noHotspot"
tssDE$hotspot[countOverlaps(tssDE,hotGR)>0]="hotspot"

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 10000, mean_mode = "w0", w = 50)

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 10000, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 10000, mean_mode = "w0", w = 50)


pdf(paste0(imageDir,"dist_FOXA1DoxUnique_ELF5_ER_20kb.pdf"),width=10,height=6)
col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1, use_raster = TRUE,col = col_fun, column_title = "ELF5Dox", name ="ELF5Dox")+
EnrichedHeatmap(mat2, use_raster = TRUE,col = col_fun, column_title = "FOXA1Dox", name ="FOXA1Dox")+
EnrichedHeatmap(mat3, use_raster = TRUE,col = col_fun, column_title = "ERDox", name ="ERDox")
dev.off()


#TAD analysis 







