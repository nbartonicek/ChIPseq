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
imageDir=paste0(resultsDir,"/figures/venn/")
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
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

cleaning="all"

pdf(paste0(imageDir,"upset_",cleaning,".pdf"),width=10,height=8)
upset(df, order.by = "freq")
dev.off()

pdf(paste0(imageDir,"venn_",cleaning,".pdf"),width=8,height=8)
venn(df)
dev.off()


##################

combined<-GRangesList()
combined[["ELF5"]]<-reduce(cleanGRsPeaks[[1]])
combined[["ER"]]<-reduce(unlist(cleanGRsPeaks[2:3]))
combined[["FOXA1"]]<-reduce(unlist(cleanGRsPeaks[4:5]))


allPeaks<-reduce(unlist(combined))
for(sampleName in names(combined)){
  cat(sampleName)
  cat("\n")
  values(allPeaks)[sampleName]<-countOverlaps(allPeaks,combined[sampleName])
}

df<-as.data.frame(values(allPeaks))

cleaning="all"
pdf(paste0(imageDir,"upset_short_",cleaning,".pdf"),width=10,height=8)
upset(df, order.by = "freq")
dev.off()


pdf(paste0(imageDir,"venn_short_",cleaning,".pdf"),width=8,height=8)
venn(df)
dev.off()

ELF5=cleanGRs[["ELF5Dox"]]
gr1kb<-resize(ELF5,fix="center",width=peakWidth)

if(cleaning=="clean"){
  data(HOT.spots)
  data(wgEncodeTfbsV3)
  hotGR <- liftover(reduce(unlist(HOT.spots)))
  wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
  removeOl <- function(.ele){
          ol <- findOverlaps(.ele, hotGR)
          if(length(ol)>0) .ele <- .ele[-unique(queryHits(ol))]
          .ele
  }
  gr1kb<-removeOl(gr1kb)
}



pdf(paste0(imageDir,"correlations_genic_Doxpeaks_binary_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)

dev.off()

####### DOX and diff peaks - separately enriched and depleted #######


allPeaks<-reduce(unlist(cleanGRsPeaks[c(1,2,4,8,12)]))
for(sampleName in names(cleanGRsPeaks[c(1,2,4,8,12)])){
  cat(sampleName)
  cat("\n")
  values(allPeaks)[sampleName]<-countOverlaps(allPeaks,cleanGRsPeaks[sampleName])
}

df<-as.data.frame(values(allPeaks))

cleaning="all"

pdf(paste0(imageDir,"upset_diff_enriched_categories_",cleaning,".pdf"),width=10,height=8)
upset(df)
dev.off()

pdf(paste0(imageDir,"venn_diff_enriched_categories_",cleaning,".pdf"),width=8,height=8)
venn(df)
dev.off()




allPeaks<-reduce(unlist(cleanGRsPeaks[c(1,3,5,9,13)]))
for(sampleName in names(cleanGRsPeaks[c(1,3,5,9,13)])){
  cat(sampleName)
  cat("\n")
  values(allPeaks)[sampleName]<-countOverlaps(allPeaks,cleanGRsPeaks[sampleName])
}

df<-as.data.frame(values(allPeaks))

cleaning="all"

pdf(paste0(imageDir,"upset_diff_depleted_categories_",cleaning,".pdf"),width=10,height=8)
upset(df)
dev.off()

pdf(paste0(imageDir,"venn_diff_depleted_categories_",cleaning,".pdf"),width=8,height=8)
venn(df)
dev.off()

####### only diff peaks - separately enriched and depleted #######


allPeaks<-reduce(unlist(cleanGRsPeaks[c(1,8,9,12,13)]))
for(sampleName in names(cleanGRsPeaks[c(1,8,9,12,13)])){
  cat(sampleName)
  cat("\n")
  values(allPeaks)[sampleName]<-countOverlaps(allPeaks,cleanGRsPeaks[sampleName])
}

df<-as.data.frame(values(allPeaks))

cleaning="all"

pdf(paste0(imageDir,"upset_diffOnly_categories_",cleaning,".pdf"),width=10,height=8)
upset(df)
dev.off()

pdf(paste0(imageDir,"venn_diffOnly_categories_",cleaning,".pdf"),width=8,height=8)
venn(df)
dev.off()








