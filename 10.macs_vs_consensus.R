

library(ChIPQC)
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
library(ggplot2)
library(rtracklayer)



#homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs/")
inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate/")

inBams=paste0(homedir,"/projects/Chris/project_results/ELF5.picard/")


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
#macsFiles<-macsFiles[grepl("FOXA1",macsFiles)]

results_separate<-GRangesList()
results_promoters<-list()
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
  results_promoters[[sampleName]]<-sum(countOverlaps(gr,tss)>0)
}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tss <- promoters(genes(txdb,columns=c("gene_id")))
#promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#Test everything on Foxa1

imageDirTest=paste0(resultsDir,"/figures/test/")
system(paste0("mkdir -p ",imageDirTest))


sampleName<-names(results)[5]
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
for(i in 1:length(samples)){
  percentFromTotal<-sum(countOverlaps(grAll[consensus==i],results[[5]])>0)/length(grAll[consensus==i])
  cat(percentFromTotal)
  cat("\n")
}


#add foxa1 data
foxa1file=paste0(homedir,"/projects/Chris/project_results/GSE81714/GSE81714_foxa1_untreated-mcf7_input_peaks.narrowPeak")
foxa1<-read.delim(foxa1file,header=F,stringsAsFactors=F)
colnames(foxa1) <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue","peak")
foxa1<-foxa1[!grepl("_",foxa1$chrom),]
gr<-GRanges(seqnames=foxa1$chrom,IRanges(start=foxa1$chromStart,end=foxa1$chromEnd),strand="*",score=foxa1$score,pval=foxa1$log10.qvalue)
gr=gr[gr$pval>=5]
#48669 peaks

ch = import.chain( "../annotation/hg19ToHg38.over.chain")
seqlevelsStyle(gr) = "UCSC"  # necessary
foxa1gr = liftOver(gr, ch)
foxa1gr = unlist(foxa1gr)
genome(foxa1gr) = "hg38"
  

sampleNames<-names(results)[4]
for(sampleName in sampleNames){
  macs=results[[sampleName]]
  macs_e5=macs[macs$pval>=5]
  cat(sampleName)
  cat("\n")
  samples<-results_separate[grepl(sampleName,names(results_separate))]
  cat(names(samples))
  cat("\n")
  grAll<-reduce(unlist(samples))
  sampleCounts<-list()
  for(i in 1:length(samples)){
    sampleCounts[[names(samples)[i]]]<-1*(countOverlaps(grAll,samples[[i]])>0)
  }
  df<-do.call("cbind",sampleCounts)
  consensus<-apply(df,1,sum)
  resL<-list()
  for(i in 1:length(samples)){
    percentFromTotal<-sum(countOverlaps(grAll[consensus==i],macs)>0)/length(grAll[consensus==i])
    percentFromMacs<-sum(countOverlaps(grAll[consensus==i],macs)>0)/length(macs)
    mat<-findOverlaps(grAll[consensus==i],macs)
    #scoreTotal<-mean(grAll[consensus==i][unique(queryHits(mat))]$score)
    scoreMacs<-mean(macs[unique(subjectHits(mat))]$score)
    
    percentFromTotal_e5<-sum(countOverlaps(grAll[consensus==i],macs_e5)>0)/length(grAll[consensus==i])
    percentFromMacs_e5<-sum(countOverlaps(grAll[consensus==i],macs_e5)>0)/length(macs_e5)
    mat<-findOverlaps(grAll[consensus==i],macs_e5)
    #scoreTotal_e5<-mean(grAll[consensus==i][unique(queryHits(mat))]$score)
    scoreMacs_e5<-mean(macs_e5[unique(subjectHits(mat))]$score)
    
    percentFromTotal_Carroll<-sum(countOverlaps(grAll[consensus==i],foxa1gr)>0)/length(grAll[consensus==i])
    percentFrom_Carroll<-sum(countOverlaps(grAll[consensus==i],foxa1gr)>0)/length(foxa1gr)
    mat<-findOverlaps(grAll[consensus==i],foxa1gr)
    #scoreTotal_Carroll<-mean(foxa1gr[consensus==i][unique(queryHits(mat))]$score)
    ScoreCarroll<-mean(foxa1gr[unique(subjectHits(mat))]$score)

    resL[[i]]<-data.frame(type=paste0(i," out of ",length(samples)),counts=length(grAll[consensus==i]),
      percentWithMacs2=sprintf("%0.3f",percentFromTotal),
      percentOfMacs2=sprintf("%0.3f",percentFromMacs),
      #scoreTotal=scoreTotal,
      scoreMacs=scoreMacs,
      percentFromTotal_e5=sprintf("%0.3f",percentFromTotal_e5),
      percentOfMAcs_e5=sprintf("%0.3f",percentFromMacs_e5),
      #scoreTotal_e5=scoreTotal_e5,
      scoreMacs_e5=scoreMacs_e5,
      percentWithCarroll=sprintf("%0.3f",percentFromTotal_Carroll),
      percentOfCarroll=sprintf("%0.3f",percentFrom_Carroll),
      #scoreTotal_Carroll=scoreTotal_Carroll,
      ScoreCarroll=ScoreCarroll
     )
    #cat(percentFromTotal)
    #cat("\n")
  }
  print(do.call("rbind",resL))
}

#add er data 1
ERfiles<-list.files(paste0(homedir,"/projects/Chris/project_results/E-GEOD-41561"),full.names=T)
erList<-GRangesList()
for(file in ERfiles){
  sampleName=basename(file)
  #er<-read.delim(foxa1file,header=F,stringsAsFactors=F)
  #colnames(er) <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue","peak")
  #er<-er[!grepl("_",er$chrom),]
  #gr<-GRanges(seqnames=foxa1$chrom,IRanges(start=er$chromStart,end=er$chromEnd),strand="*",score=er$score,pval=er$log10.qvalue)
  gr<-import(file)
  gr=gr[gr$score>=50]
  erList[[sampleName]]<-gr
}
common<-reduce(unlist(erList))
co1<-1*(countOverlaps(common,erList[[1]])>0)
co2<-1*(countOverlaps(common,erList[[2]])>0)
er=common[(co1+co2)>1]

ch = import.chain( "../annotation/hg18ToHg38.over.chain")
seqlevelsStyle(er) = "UCSC"  # necessary
er = liftOver(er, ch)
er = unlist(er)
genome(er) = "hg38"
temp=er
#add er data 1

#ERfiles<-list.files(paste0(homedir,"/projects/Chris/project_results/E-MTAB-740"),full.names=T)
#erList<-GRangesList()
#for(file in ERfiles){
#  sampleName=basename(file)
#  #er<-read.delim(foxa1file,header=F,stringsAsFactors=F)
#  #colnames(er) <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue","peak")
#  #er<-er[!grepl("_",er$chrom),]
#  #gr<-GRanges(seqnames=foxa1$chrom,IRanges(start=er$chromStart,end=er$chromEnd),strand="*",score=er$score,pval=er$log10.qvalue)
#  gr<-import(file)
#  values(gr)<-NULL
#  #gr=gr[gr$score>=50]
#  erList[[sampleName]]<-gr
#}
#common<-reduce(unlist(erList))
#co1<-1*(countOverlaps(common,erList[[1]])>0)
#co2<-1*(countOverlaps(common,erList[[2]])>0)
#er=common[(co1+co2)>1]
#
#
#
#ch = import.chain( "../annotation/hg19ToHg38.over.chain")
#seqlevelsStyle(er) = "UCSC"  # necessary
#er = liftOver(er, ch)
#er = unlist(er)
#genome(er) = "hg38"



sampleNames<-names(results)[3]
for(sampleName in sampleNames){
  macs=results[[sampleName]]
  macs_e5=macs[macs$pval>=5]
  cat(sampleName)
  cat("\n")
  samples<-results_separate[grepl(sampleName,names(results_separate))]
  cat(names(samples))
  cat("\n")
  grAll<-reduce(unlist(samples))
  sampleCounts<-list()
  for(i in 1:length(samples)){
    sampleCounts[[names(samples)[i]]]<-1*(countOverlaps(grAll,samples[[i]])>0)
  }
  df<-do.call("cbind",sampleCounts)
  consensus<-apply(df,1,sum)
  resL<-list()
  for(i in 1:length(samples)){
    percentFromTotal<-sum(countOverlaps(grAll[consensus==i],macs)>0)/length(grAll[consensus==i])
    percentFromMacs<-sum(countOverlaps(grAll[consensus==i],macs)>0)/length(macs)
    mat<-findOverlaps(grAll[consensus==i],macs)
    #scoreTotal<-mean(grAll[consensus==i][unique(queryHits(mat))]$score)
    scoreMacs<-mean(macs[unique(subjectHits(mat))]$score)
    
    percentFromTotal_e5<-sum(countOverlaps(grAll[consensus==i],macs_e5)>0)/length(grAll[consensus==i])
    percentFromMacs_e5<-sum(countOverlaps(grAll[consensus==i],macs_e5)>0)/length(macs_e5)
    mat<-findOverlaps(grAll[consensus==i],macs_e5)
    #scoreTotal_e5<-mean(grAll[consensus==i][unique(queryHits(mat))]$score)
    scoreMacs_e5<-mean(macs_e5[unique(subjectHits(mat))]$score)
    
    percentFromTotal_Carroll<-sum(countOverlaps(grAll[consensus==i],er)>0)/length(grAll[consensus==i])
    percentFrom_Carroll<-sum(countOverlaps(grAll[consensus==i],er)>0)/length(er)
    mat<-findOverlaps(grAll[consensus==i],er)
    #scoreTotal_Carroll<-mean(foxa1gr[consensus==i][unique(queryHits(mat))]$score)
    ScoreCarroll<-mean(er[unique(subjectHits(mat))]$score)

    resL[[i]]<-data.frame(type=paste0(i," out of ",length(samples)),counts=length(grAll[consensus==i]),
      percentWithMacs2=sprintf("%0.3f",percentFromTotal),
      percentOfMacs2=sprintf("%0.3f",percentFromMacs),
      #scoreTotal=scoreTotal,
      scoreMacs=scoreMacs,
      percentFromTotal_e5=sprintf("%0.3f",percentFromTotal_e5),
      percentOfMAcs_e5=sprintf("%0.3f",percentFromMacs_e5),
      #scoreTotal_e5=scoreTotal_e5,
      scoreMacs_e5=scoreMacs_e5,
      percentWithCarroll=sprintf("%0.3f",percentFromTotal_Carroll),
      percentOfCarroll=sprintf("%0.3f",percentFrom_Carroll),
      #scoreTotal_Carroll=scoreTotal_Carroll,
      ScoreCarroll=ScoreCarroll
     )
    #cat(percentFromTotal)
    #cat("\n")
  }
  print(do.call("rbind",resL))
}

