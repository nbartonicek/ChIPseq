#enhancers were downloaded from this paper: https://www.nature.com/articles/s41598-017-02257-3#Sec26


library(GenomicRanges)
library(ChIPpeakAnno)
library(ggplot2)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(rGREAT)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(hier.part)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(plyr)
library(TFBSTools)
library(JASPAR2018)
library(Biostrings)

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
imageDir=paste0(resultsDir,"/figures/enhancers/")
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")

system(paste0("mkdir -p ",imageDir))

#first import the enhancers
enh<-read.table(paste0(annotationDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
enh=enh[enh$type %in% "Distal_TE",]

enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end),type=enh$type)

liftover<-function(gr){
  gr.hg19=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg19ToHg38.over.chain"))
  seqlevelsStyle(gr.hg19) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg19, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}

enh.hg38<-liftover(enh.hg19)

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


#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")), upstream=10001, downstream=10000)


##### fetch enhancer peaks sequences, centered on ELF5Dox motif #####

mat<-findOverlaps(cleanGRs[["ELF5Dox"]],enh.hg38)
gr<-cleanGRs[["ELF5Dox"]][unique(queryHits(mat))]


gr1kb<-resize(gr,fix="center",width=1000)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/intergenic_ELF5motif_notcentered_all_1000.bed"))

################################################
######## repeat the same for promoters #########
################################################

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")), upstream=2000, downstream=500)


##### fetch enhancer peaks sequences, centered on ELF5Dox motif #####

mat<-findOverlaps(cleanGRs[["ELF5Dox"]],promoter)
gr<-cleanGRs[["ELF5Dox"]][unique(queryHits(mat))]


gr1kb<-resize(gr,fix="center",width=1000)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/genic_ELF5motif_notcentered_all_1000.bed"))



################################################
######## repeat the same for tamoxifen resistance #########
################################################



scriptsPath=paste0(projectDir,"/scripts/")
inPath=paste0(homedir,"/projects/Chris/project_results/GSE75372/")

sampleID="GSE75372"
ChIPs<-read.table(paste0(scriptsPath,sampleID,"_ids.txt"))
ChIPs$type<-c(rep("sensitive",1),rep("resistant",1))

peaks<-GRangesList()
for(i in 1:length(ChIPs$V2)){
  file=ChIPs$V2[i]
  sampleName=as.character(ChIPs$V1[i])
  cat(sampleName)
  cat("\n")
  gr<-import(paste0(inPath,file))
  #gr<-GRanges(seqnames=df$chr,IRanges(df$start,df$end))
  peaks[[sampleName]]<-gr
}

allPeaks<-unlist(peaks)
allPeaks<-reduce(allPeaks)
peakCount<-list()
for(peak in names(peaks)){
  peakCount[[peak]]<-1*(countOverlaps(allPeaks,peaks[[peak]])>0)
}
df<-do.call("cbind",peakCount)
df=as.data.frame(df)
#find those are present in resistant or absent in resistant
df$resistance=df[,1]-df[,2]
#there are 144571 peaks
#the sensitive (present in >1 out of 3)

sensitive=allPeaks[df$resistance==1]
resistant=allPeaks[df$resistance==-1]
#27227

#liftover
liftover<-function(gr){
  gr.hg19=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg19ToHg38.over.chain"))
  seqlevelsStyle(gr.hg19) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg19, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}

sensitiveLift<-liftover(sensitive)
resistantLift<-liftover(resistant)


data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)))
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
removeOl <- function(.ele){
        ol <- findOverlaps(.ele, hotGR)
        if(length(ol)>0) .ele <- .ele[-unique(queryHits(ol))]
        .ele
}
sensitiveLiftClean<-removeOl(sensitiveLift)
resistantLiftClean<-removeOl(resistantLift)
values(sensitiveLiftClean)<-NULL
values(resistantLiftClean)<-NULL



#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


##### fetch enhancer peaks sequences, centered on ELF5Dox motif #####
tam<-GRangesList()
tam[["sensitive"]]=sensitiveLiftClean
tam[["resistant"]]=resistantLiftClean

save(tam,file=paste0(robjectsDir,"tamoxifen.Rdata"))

for(type in names(tam)){

  mat<-findOverlaps(cleanGRs[["ELF5Dox"]],tam[[type]])
  gr<-cleanGRs[["ELF5Dox"]][unique(queryHits(mat))]


  gr1kb<-resize(gr,fix="center",width=1000)

  seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
  names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
  Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/tam",type,"_ELF5motif_notcentred_all_1000.bed"))

}



































