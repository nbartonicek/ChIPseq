#data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66733

library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(corrplot)
library(ChIPpeakAnno)

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
imageDir=paste0(resultsDir,"/figures/grammar/")
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
grammarDir=paste0(resultsDir,"/grammar/")
system(paste0("mkdir -p ",imageDir))
permRobjectsDir=paste(resultsDir,"/Robjects/grammarPermutations/",sep="")
system(paste0("mkdir -p ",permRobjectsDir))

liftover<-function(gr,genome){
  if(genome=="hg38"){return(gr)} else {
    gr.old=gr
    ch = import.chain(paste0(projectDir,"/annotation/",genome,"ToHg38.over.chain"))
    seqlevelsStyle(gr.old) = "UCSC"  # necessary
    gr.hg38 = liftOver(gr.old, ch)
    gr.hg38  = unlist(gr.hg38 )
    genome(gr.hg38) = "hg38"
    return(gr.hg38)
  }
}

load(paste0(grammarDir,"all.Rdata"))

#find consensus of the same peaks
sampleNames<-names(grL)
replicates<-gsub("\\.r.*","",sampleNames)

replicateL<-split(grL,replicates)

consensus<-function(gr1,gr2){
  total<-reduce(c(gr1,gr2))
  mat<-findOverlaps(total,gr1)
  total<-total[unique(queryHits(mat))]  
  mat<-findOverlaps(total,gr2)
  total<-total[unique(queryHits(mat))]
  return(total)
}

consensusReps<-lapply(replicateL,function(x){
  if(length(x)>1){
    a=consensus(x[[1]],x[[2]]);return(a)
  } else {
    return(unlist(x))
  }
}
)
consensusReps=GRangesList(consensusReps)

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
  #values(grPos)<-NULL
  #values(grNeg)<-NULL
  cleanGRsPeaks[[paste0(sampleName,"_Dox_enriched")]]<-grPos
  cleanGRsPeaks[[paste0(sampleName,"_Dox_depleted")]]<-grNeg
}

for(sampleName in names(cleanGRsPeaks)[c(2:5,8:9,12:13)]){
  gr=cleanGRsPeaks[[sampleName]]
  values(gr)<-NULL
  consensusReps[[sampleName]]=gr
}

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

#first, how many of them overlap hotspots
data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftver(reduce(unlist(HOT.spots)),"hg19")
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3,"hg19")   


#now for each in enhancer, superenhancer, DE, HOT, resistant, sensitive
#get p-values with ELF5

#plot profile

targets<-GRangesList()
for(class in c("enhancer","superenhancer","HOT","DE","resistant","sensitive")){
      cat(".")
      if(class=="enhancer"){
        #first import the enhancers
        enh<-read.table(paste0(annotationDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
        enh=enh[enh$type %in% "Distal_TE",]
        enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end))
        grTarget<-liftover(enh.hg19,"hg19")
        targets[[class]]<-grTarget
      } else if (grepl("superenhancer",class)){
        enh<-read.table(paste0(annotationDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
        enh=enh[enh$type %in% "Distal_SE",]
        enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end))
        grTarget<-liftover(enh.hg19,"hg19")
        targets[[class]]<-grTarget


      } else if (grepl("HOT",class)){
        #TSS distribution
        targets[[class]]<-hotGR

      } else if (grepl("DE",class)){
        #TSS distribution
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        tss <- promoters(genes(txdb,columns=c("gene_id")), upstream=5000, downstream=5000)
        #promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
        geneid <- mapIds(org.Hs.eg.db, names(tss), "SYMBOL","ENTREZID")
        tss=tss[!duplicated(geneid)]
        names(tss)<-geneid[!duplicated(geneid)]

        #DE genes
        expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
        DEgenes<-unique(expression$hgnc_symbol[expression$adj.P.Val<0.05&abs(expression$logFC)>1])
        tssDE<-tss[names(tss) %in% DEgenes]
        grTarget <- tssDE
        values(grTarget)<-NULL
        targets[[class]]<-grTarget
      } else if (class=="resistant") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
        values(grTarget)<-NULL
        targets[[class]]<-grTarget
      } else if (class=="sensitive") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
        values(grTarget)<-NULL
        targets[[class]]<-grTarget
      } 
}



#Now for each class get p-values
#for(class in names(targets)){
for(class in c("DE")){
  cat(class)
  cat("\n")
  for(type in names(consensusReps)){
    cat(".")
    pool <- new("permPool", grs=GRangesList(targets[[class]]), N=length(consensusReps[[type]]))
    pt <- peakPermTest(resize(cleanGRs[["ELF5Dox"]],1000,fix="center"), consensusReps[[type]], pool=pool, ntimes=100)
    save(pt,file=paste0(permRobjectsDir,"peakOverlap_",class,"_",type,".Rdata"))
  }
  cat("\n")
}

res<-list()
for(file in list.files(permRobjectsDir,full.names=T)){
  sampleName<-gsub(".Rdata","",basename(file))
  sampleName<-gsub("peakOverlap_DE_","",sampleName)
  load(file)
  res[[sampleName]]<- pt$cntOverlaps$pval
}




























































