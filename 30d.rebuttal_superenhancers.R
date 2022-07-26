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
library(reshape2)
library(pheatmap)
library(UpSetR) 

#BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38","data.table","corrplot","ChIPpeakAnno","reshape2","pheatmap","UpSetR"))

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
chrGR<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))
homedir="/share/ScratchGeneral/nenbar"
#homedir="../../.."
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
permRobjectsDir=paste(resultsDir,"/Robjects/grammarPermutations/",sep="")
grammarBedDir<-paste0(resultsDir,"/grammarBeds/")

files<-list.files(grammarBedDir,full.names=T,pattern="bed")
grL<-GRangesList()
for(file in files){
	sampleName=gsub(".bed","",basename(file))
	cat(sampleName)
	cat("\n")
	grL[[sampleName]]<-import(file)
}
save(grL,file=paste0(grammarDir,"all.Rdata"))
#load(paste0(grammarDir,"all_consensus_clean.Rdata"))
load(paste0(grammarDir,"targets.Rdata"))
meds<-targets[grepl("MED1\\.|Med1\\.",names(targets))]

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))



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

save(targets,file=paste0(grammarDir,"targets.Rdata"))

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
load(paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))

mat<-findOverlaps(targets[[2]],meds)
targets[[2]][unique(queryHits(mat))]



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



targets[[2]][unique(queryHits(mat))]
ELF5SE=targets[["superenhancer"]][unique(queryHits(mat))]

consensusRepsShort=consensusReps[grepl("vehicle|estradiol|BORIS",names(consensusReps))]
seCounts<-countOverlaps(ELF5SE,consensusRepsShort)
ELF5SEover20<-ELF5SE[seCounts>20]

seELF5Counts<-countOverlaps(ELF5SEover20,meds)








