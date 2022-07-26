library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)

homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/GSE32222/")
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")

ChIPs<-read.table("GSE32222_ids.txt")
ChIPs$type<-c(rep("sensitive",7),rep("resistant",4))

peaks<-GRangesList()
for(i in 1:length(ChIPs$V2)){
  file=ChIPs$V2[i]
  sampleName=as.character(ChIPs$V1[i])
  cat(sampleName)
  cat("\n")
  df<-read.table(paste0(inPath,file),header=T)
  gr<-GRanges(seqnames=df$chr,IRanges(df$start,df$end))
  peaks[[sampleName]]<-gr
}

allPeaks<-unlist(peaks)
allPeaks<-reduce(allPeaks)
peakCount<-list()
for(peak in names(peaks)){
  peakCount[[peak]]<-1*(countOverlaps(allPeaks,peaks[[peak]])>0)
}
df<-do.call("cbind",peakCount)

#find those are present in resistant or absent in resistant
countsRes<-apply(df[,1:7],1,sum)

#there are 144571 peaks
#the sensitive (present in >3 out of 7)
cat(sum(countsRes>3))
#38867
sensitive=allPeaks[countsRes>3]
#the resistant (present in >2 out of 4)
countsRes<-apply(df[,8:11],1,sum)
cat(sum(countsRes>2))
resistant=allPeaks[countsRes>2]
#27227

#liftover
liftover<-function(gr){
  gr.hg18=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg18ToHg38.over.chain"))
  seqlevelsStyle(gr.hg18) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg18, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}

sensitiveLift<-liftover(sensitive)
resistantLift<-liftover(resistant)

mat<-findOverlaps(sensitiveLift,resistantLift)
sensitiveUniq<-sensitiveLift[-unique(queryHits(mat))]
resistantUniq<-resistantLift[-unique(subjectHits(mat))]
#then find whether they are enriched or depleted in our peaks

#load the peaks and differentially expressed genes
#load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
load(paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))

for(file in list.files(robjectsDir,pattern="Rdata",full.names=T)){
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

ERpeaks<-cleanGRsPeaks[c(1,2,3,6,7,8,9,10,11,12,13)]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

sensitive<-enrichPeakOverlap(queryPeak     = sensitiveUniq,
                  targetPeak    = ERpeaks,
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 5000,
                  chainFile     = NULL,
                  verbose       = FALSE)


resistant<-enrichPeakOverlap(queryPeak     = resistantUniq,
                  targetPeak    = ERpeaks,
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 5000,
                  chainFile     = NULL,
                  verbose       = FALSE)










