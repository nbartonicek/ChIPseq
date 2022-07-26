

#library(ChIPQC)
library(GenomicRanges)
library(VennDiagram)
#library(ChIPpeakAnno)
library(ggplot2)
library(rtracklayer)
library(ggfortify)



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
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")

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


macsFiles<-list.files(inPath,pattern="summits",full.names=T)
#macsFiles<-macsFiles[!grepl("gapped",macsFiles)]

#First take all combined peaks
results<-GRangesList()
for(macsFile in macsFiles){

  mypeaks <- import(macsFile)
  sampleName=gsub("_summits.bed","",basename(macsFile))
  cat(sampleName)
  cat("\n")
  mypeaks<-mypeaks[!grepl("_",as.character(seqnames(mypeaks))),]
  gr<-GRanges(seqnames=seqnames(mypeaks),IRanges(start=start(mypeaks),end=end(mypeaks)),strand="*",name=values(mypeaks)$name,score=values(mypeaks)$score)
  names(gr)<-1:length(gr)
  results[[sampleName]]<-gr
}

macsFiles<-list.files(inPath_separate,pattern="Peak",full.names=T)
macsFiles<-macsFiles[!grepl("gapped",macsFiles)]


#Then take all separate peaks
results_separate<-GRangesList()
for(macsFile in macsFiles){

  mypeaks <- read.delim(macsFile,header=F,stringsAsFactors=F)
  sampleName=gsub("_peaks.*Peak","",basename(macsFile))
  cat(sampleName)
  cat("\n")
  colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
  mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
  gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*",score=mypeaks$score,pval=mypeaks$log10.qvalue)
  names(gr)<-1:length(gr)
  results_separate[[sampleName]]<-gr
  cat(mean(width(gr)))
  cat("\n")

}

#then take summits of all combined peaks for the individual ones that have 3/4

cleanGRs<-GRangesList()
for(sampleName in names(results)){
  cat(sampleName)
  cat("\n")
  samples<-results_separate[grepl(sampleName,names(results_separate))]
  #fix Reduce - all or "inclusive"
  grAll<-reduce(unlist(samples))
  sampleCounts<-list()
  for(i in 1:length(samples)){
    sampleCounts[[names(samples)[i]]]<-1*(countOverlaps(grAll,samples[[i]])>0)
  }
  df<-do.call("cbind",sampleCounts)
  consensus<-apply(df,1,sum)
  grAll$consensus=consensus

  cutoff=3
  if(grepl("H3K4",sampleName)){cutoff=2}
  grAll=grAll[grAll$consensus>=cutoff]
  gr=grAll
  values(gr)="null"
  
  summits<-results[[sampleName]]
  mat<-findOverlaps(summits,grAll)
  counts<-table(countOverlaps(summits,grAll))
  cat(as.numeric(counts))

  summits<-summits[unique(queryHits(mat))]
  cleanGRs[[sampleName]]=summits
  save(summits,file=paste0(cleanRobjectsDir,sampleName,".Rdata"))
}
save(cleanGRs,file=paste0(cleanRobjectsDir,"all_peaks_summits.Rdata"))

#create the same but not for summits but peak regions


macsFiles<-list.files(inPath,pattern="Peak",full.names=T)
macsFiles<-macsFiles[!grepl("gapped",macsFiles)]

hists<-list()
results<-GRangesList()
for(macsFile in macsFiles){

  mypeaks <- read.delim(macsFile,header=F,stringsAsFactors=F)
  sampleName=gsub("_peaks.*Peak","",basename(macsFile))
  cat(sampleName)
  cat("\n")
  colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
  mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
  gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*",score=mypeaks$score,pval=mypeaks$log10.qvalue)
  names(gr)<-1:length(gr)
  results[[sampleName]]<-gr
 }

#ERNoDox3 elimination as an outlier
#results_separate=results_separate[-11]

cleanGRsPeaks<-GRangesList()
for(sampleName in names(results)){
  cat(sampleName)
  cat("\n")
  samples<-results_separate[grepl(sampleName,names(results_separate))]
  #fix Reduce - all or "inclusive"
  grAll<-reduce(unlist(samples))
  sampleCounts<-list()
  for(i in 1:length(samples)){
    sampleCounts[[names(samples)[i]]]<-1*(countOverlaps(grAll,samples[[i]])>0)
  }
  df<-do.call("cbind",sampleCounts)
  consensus<-apply(df,1,sum)
  grAll$consensus=consensus

  cutoff=3
  if(grepl("H3K4",sampleName)){cutoff=2}
  #if(grepl("ER",sampleName)){cutoff=2}

  grAll=grAll[grAll$consensus>=cutoff]
  gr=grAll
  values(gr)="null"
  
  summits<-results[[sampleName]]
  mat<-findOverlaps(summits,grAll)
  counts<-table(countOverlaps(summits,grAll))
  cat(as.numeric(counts))

  summits<-summits[unique(queryHits(mat))]
  cleanGRsPeaks[[sampleName]]=summits

}
save(cleanGRsPeaks,file=paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))
#save(cleanGRsPeaks,file=paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))



