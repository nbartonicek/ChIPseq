library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(R.utils)
library(ChIPpeakAnno)

homedir="/share/ScratchGeneral/nenbar"
inPath=paste0(homedir,"/projects/Chris/project_results/GSE32222/")
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste0(projectDir,"/scripts/")
  args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["id"]])){id = args$id} 
if (!is.null(args[["type"]])){type = args$type} 
sampleID="GSE32222"
ChIPs<-read.table(paste0(scriptsPath,sampleID,"_ids.txt"))
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
countsRes<-apply(df[,1:3],1,sum)

#there are 144571 peaks
#the sensitive (present in >1 out of 3)
cat(sum(countsRes>1))
#57283
sensitive=allPeaks[countsRes>1]
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

ERpeaks<-cleanGRsPeaks[c(1,2,3,8,9,10,11,12,13)]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#sensitive<-enrichPeakOverlap(queryPeak     = sensitiveUniq,
#                  targetPeak    = ERpeaks,
#                  TxDb          = txdb,
#                  pAdjustMethod = "BH",
#                  nShuffle      = 10000,
#                  chainFile     = NULL,
#                  verbose       = FALSE)
#
#sensitive
#
#resistant<-enrichPeakOverlap(queryPeak     = resistantUniq,
#                  targetPeak    = ERpeaks,
#                  TxDb          = txdb,
#                  pAdjustMethod = "BH",
#                  nShuffle      = 10000,
#                  chainFile     = NULL,
#                  verbose       = FALSE)
#
#resistant



data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)))
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
removeOl <- function(.ele){
        ol <- findOverlaps(.ele, hotGR)
        if(length(ol)>0) .ele <- .ele[-unique(queryHits(ol))]
        .ele
}
sensitiveLiftClean<-removeOl(sensitiveUniq)
resistantLiftClean<-removeOl(resistantUniq)
values(sensitiveLiftClean)<-NULL
values(resistantLiftClean)<-NULL


poolL<-list()
poolL[["sensitive"]] <- new("permPool", grs=GRangesList(wgEncodeTfbsV3), N=length(sensitiveLiftClean))
poolL[["resistant"]] <- new("permPool", grs=GRangesList(wgEncodeTfbsV3), N=length(resistantLiftClean))

id=as.numeric(id)

resPeaks<-GRangesList()
resPeaks[["sensitive"]]<-sensitiveLiftClean
resPeaks[["resistant"]]<-resistantLiftClean


pt <- peakPermTest(resPeaks[[type]], ERpeaks[[id]], pool=poolL[[type]], ntimes=10000)
save(pt,file=paste0(robjectsDir,sampleID,"_peakOverlap_",type,"_",names(ERpeaks)[id],".Rdata"))








