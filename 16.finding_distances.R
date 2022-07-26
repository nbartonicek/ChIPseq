

library(GenomicRanges)
library(ggplot2)
library(rtracklayer)

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
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")

system(paste0("mkdir -p ",cleanRobjectsDir))

#1. clean up from blacklist (186 in merged peaks, 50 in macs, 156 in Carroll data)
#2. plot bar plot for distribution and percentage overlap
#3. plot intensity/p-value/fdr vs presence in consensus/macs2/macs2 strict
#4. what is macs2 strict?
#liftover of blacklist sections
#data(blacklist_hg19)
#path = system.file(package="rtracklayer", "extdata", "../annotation/hg38ToHg19.over.chain")
#ch = import.chain( "../annotation/hg19ToHg38.over.chain")
#seqlevelsStyle(blacklist.hg19) = "UCSC"  # necessary
#blacklist = liftOver(blacklist.hg19, ch)
#blacklist = unlist(blacklist)
#genome(blacklist) = "hg38"

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
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

ELF5=cleanGRs[["ELF5Dox"]]

sampleName="FOXA1"

cleanGRs=cleanGRs[grepl(sampleName,names(cleanGRs))]
sapply(cleanGRs,length)

sampleName=names(cleanGRs)[1]
temp=ELF5
ELF5=resize(ELF5,width(ELF5)+1000,fix="center")
mat<-findOverlaps(ELF5,cleanGRs[[sampleName]])
ELF5S<-temp[unique(queryHits(mat))]
gr<-cleanGRs[[sampleName]][unique(subjectHits(mat))]

results<-list()
for(i in 0:1000){
  cat(".")
  results[[i+1]]<-sum(countOverlaps(shift(ELF5S,i),gr))
}


sampleName=names(cleanGRs)[2]
temp=ELF5
#ELF5=resize(ELF5,width(ELF5)+1000,fix="center")
mat<-findOverlaps(ELF5,cleanGRs[[sampleName]])
ELF5S<-temp[unique(queryHits(mat))]
gr<-cleanGRs[[sampleName]][unique(subjectHits(mat))]

results2<-list()
for(i in 0:1000){
  cat(".")
  results2[[i+1]]<-sum(countOverlaps(shift(ELF5S,i),gr))
}


pdf("testDox.pdf") 
plot(unlist(results),type="l",ylab="Frequency",xlab="distance from ELF5 site",main=paste0(sampleName," distance from ELF5 peaks"))
#points(unlist(results2),type="l",col="red")
dev.off()

results<-list()
for(i in 0:1000){
  cat(".")
  results[[i+1]]<-sum(countOverlaps(shift(ELF5S,i),gr)>0)
}
pdf("testDox_unique.pdf") 
plot(unlist(results),type="l",ylab="Frequency",xlab="distance from ELF5 site",main=paste0(sampleName," distance from ELF5 peaks"))
dev.off()
