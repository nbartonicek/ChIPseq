library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(R.utils)
library(ChIPpeakAnno)

homedir="/share/ScratchGeneral/nenbar"
inPath=paste0(homedir,"/projects/Chris/project_results/GSE30407/")
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste0(projectDir,"/scripts/")
sampleID="GSE30407"


files<-list.files(inPath,pattern="csv",full.names=T)
files<-files[grepl("combined",files)]
TDpeaks<-read.table(files,header=T,sep=",")
TDpeaks<-GRanges(seqnames=TDpeaks$chr,IRanges(start=TDpeaks$start,end=TDpeaks$end))


all<-list.files(inPath,pattern="tsv",full.names=T)
TDpeaks<-read.table(all,header=T,sep="\t")
TDpeaksGR<-GRanges(seqnames=TDpeaks$space,IRanges(start=TDpeaks$start,end=TDpeaks$end))

TDpeaks24<-TDpeaksGR[TDpeaks$T24.FDR<=0.05]
TDpeaks48<-TDpeaksGR[TDpeaks$T48.FDR<=0.05]
TDpeaksBoth<-reduce(c(TDpeaks24,TDpeaks48))

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


tdLift<-liftover(TDpeaksBoth)
TDpeaks24<-liftover(TDpeaks24)
TDpeaks48<-liftover(TDpeaks48)
TDpeaksBoth<-liftover(TDpeaksBoth)

#load the peaks and differentially expressed genes
#load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
load(paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))

sampleName<-"ELF5Dox"
mcf7=cleanGRsPeaks[[sampleName]]


data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)))
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
removeOl <- function(.ele){
        ol <- findOverlaps(.ele, hotGR)
        if(length(ol)>0) .ele <- .ele[-unique(queryHits(ol))]
        .ele
}
TDliftClean<-removeOl(tdLift)
values(TDliftClean)<-NULL


poolL<-list()
poolL[["TD"]] <- new("permPool", grs=GRangesList(wgEncodeTfbsV3), N=length(TDpeaksBoth))

pt <- peakPermTest(TDpeaksBoth, mcf7, pool=poolL[["TD"]], ntimes=1000)
save(pt,file=paste0(robjectsDir,"TD_peakOverlap.Rdata"))


makeVennDiagram(list(peaks1, peaks2, peaks1, peaks2), 
                  NameOfPeaks=c("TF1", "TF2","TF3", "TF4"), 
                  totalTest=100, by="feature",
                  main = "Venn Diagram for 4 peak lists",
                  fill=c(1,2,3,4))


pdf(paste0(imageDir,"/T47D_MCF7_ELF5_venn.pdf"),width=6,height=5)
makeVennDiagram(list(TDpeaksBoth, mcf7), 
                  NameOfPeaks=c("T47D", "MCF7"), 
                  by="region",
                  main = "ELF5 peaks in T47D and MCF7 cell line",
                  fill=c(1,2))
dev.off()




pdf(paste0(imageDir,"/T47D_MCF7_ELF5_venn_all.pdf"),width=6,height=5)
makeVennDiagram(list(TDpeaks24,TDpeaks48, mcf7), 
                  NameOfPeaks=c("T47D_24","T47D_48", "MCF7"), 
                  by="region",
                  main = "ELF5 peaks in T47D and MCF7 cell line",
                  fill=c(1,2,3))
dev.off()


#genes that overlap:
mat<-findOverlaps(TDpeaksBoth,mcf7)
TDpeaksBothOL<-TDpeaksBoth[unique(queryHits(mat))]

library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")), upstream=10001, downstream=10000)
peakAnno <- annotatePeak(TDpeaksBothOL, tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
genesChip<-unique(peakAnno@anno$SYMBOL)
write.table(genesChip,"T479_genes.txt",sep="\t",quote=F)


expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]
deGenesAll<-unique(as.character(expression$hgnc_symbol[(abs(expression$logFC)>=log2(1.5)) & (expression$adj.P.Val<=0.05)]))
deGenesPos<-unique(as.character(expression$hgnc_symbol[expression$logFC>=log2(1.5) & expression$adj.P.Val<=0.05]))
deGenesNeg<-unique(as.character(expression$hgnc_symbol[expression$logFC<=-log2(1.5) & expression$adj.P.Val<=0.05]))

deGenesL<-list()
deGenesL[["DE"]]<-deGenesAll
deGenesL[["pos"]]<-deGenesPos
deGenesL[["neg"]]<-deGenesNeg

for(type in names(deGenesL)){
  cat(genesChip[genesChip %in% deGenesL[[type]]])
  cat("\n")
}

#reverse: overlap of mcf7 chips with DE genes from patients.
DEgenes4<-res[res$padj<0.1,]

