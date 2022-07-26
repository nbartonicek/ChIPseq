

library(ChIPQC)
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
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
projectDir=paste0(homedir,"/projects/claudia")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")


#1. clean up from blacklist (186 in merged peaks, 50 in macs, 156 in Carroll data)
#2. plot bar plot for distribution and percentage overlap
#3. plot intensity/p-value/fdr vs presence in consensus/macs2/macs2 strict
#4. what is macs2 strict?
#liftover of blacklist sections
data(blacklist_hg19)
path = system.file(package="rtracklayer", "extdata", "../annotation/hg38ToHg19.over.chain")
ch = import.chain( "../annotation/hg19ToHg38.over.chain")
seqlevelsStyle(blacklist.hg19) = "UCSC"  # necessary
blacklist = liftOver(blacklist.hg19, ch)
blacklist = unlist(blacklist)
genome(blacklist) = "hg38"



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


macsFiles<-list.files(inPath_separate,pattern="Peak",full.names=T)
macsFiles<-macsFiles[!grepl("gapped",macsFiles)]

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

}

grTotal<-reduce(unlist(results_separate))

pcaTotal<-list()
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

  #for each of the samples write the clean results
  for(sampleOld in names(samples)){
    cat(sampleOld)
    cat("\n")
    fileOld<-macsFiles[grepl(sampleOld,macsFiles)]
    mypeaks <- read.delim(fileOld,header=F,stringsAsFactors=F)
    colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
    mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
    gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*",score=mypeaks$score,pval=mypeaks$log10.qvalue)
    names(gr)<-1:length(gr)
    mat<-findOverlaps(gr,grAll)
    mypeaksClean<-mypeaks[unique(queryHits(mat)),]
    pcaTotal[[sampleOld]]<-countOverlaps(grTotal,gr)>0
    write.table(mypeaksClean,paste0(outPath,basename(fileOld)),row.names=F,quote=F,col.names=F,sep="\t")
  }
}

df<-data.frame(samples=basename(macsFiles),type=gsub("Dox.*","Dox",basename(macsFiles)))

totMat<-do.call("cbind",pcaTotal)
#pc<-princomp(totMat)
pdf("pca_total_binary.pdf",width=8,height=6)
autoplot(prcomp(t(totMat)), data = df, colour = 'type',label = TRUE, label.size = 3)
dev.off()


pcaTotal<-list()
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

  #for each of the samples write the clean results
  for(sampleOld in names(samples)){
    cat(sampleOld)
    cat("\n")
    fileOld<-macsFiles[grepl(sampleOld,macsFiles)]
    mypeaks <- read.delim(fileOld,header=F,stringsAsFactors=F)
    colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
    mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
    gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*",score=mypeaks$score,pval=mypeaks$log10.qvalue)
    names(gr)<-1:length(gr)
    mat<-findOverlaps(gr,grAll)
    mypeaksClean<-mypeaks[unique(queryHits(mat)),]
    gr=gr[unique(queryHits(mat))]
    pcaTotal[[sampleOld]]<-countOverlaps(grTotal,gr)>0
  }
}

df<-data.frame(samples=basename(macsFiles),type=gsub("Dox.*","Dox",basename(macsFiles)))
totMat<-do.call("cbind",pcaTotal)
#pc<-princomp(totMat)
pdf(paste0("pca_total_cutoff",cutoff,".pdf"),width=8,height=6)
autoplot(prcomp(t(totMat)), data = df, colour = 'type',label = TRUE, label.size = 3)
dev.off()

