library(GenomicRanges)
library(rtracklayer)

homedir="/share/ScratchGeneral/nenbar"
inPath=paste0(homedir,"/projects/Chris/project_results/ELF5_RNAseq.star_singlemap_chimeras/")
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste0(projectDir,"/scripts/")
 

chimeras<-GRangesList()
sampleNames<-list.files(inPath,full.names=T)
for(sampleName in sampleNames){
  sampleID=basename(sampleName)
  sampleID=gsub("MCF7_2b_V5_","",sampleID)
  sampleID=gsub("_.*","",sampleID)
  cat(sampleID)
  cat("\n")
  df=read.table(paste0(sampleName,"/Chimeric.out.junction"))
  gr1<-GRanges(seqnames=df$V1,IRanges(df$V2,width=1))
  gr2<-GRanges(seqnames=df$V4,IRanges(df$V5,width=1))
  gr<-reduce(c(gr1,gr2))
  chimeras[[sampleID]]<-gr
}



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

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

res<-list()
for(sampleName in names(cleanGRsPeaks)){
  cat(sampleName)
  cat("\n")
  gr=cleanGRsPeaks[[sampleName]]
  gr=resize(gr,1000,fix="center")
  res[[sampleName]]<-countOverlaps(chimeras,cleanGRsPeaks[[sampleName]])
}
res=res[1:15]
df<-do.call("rbind",res)
lengths<-sapply(cleanGRsPeaks[1:15],length)
normDf<-apply(df,2,function(x){100*x/lengths})
write.table(df,"../project_results/tables/chimeras.xls",row.names=T,quote=F,sep="\t")

