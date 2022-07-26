

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
library(fgsea)
library(stringr)


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

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")), upstream=10001, downstream=10000)
#promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#first get peaks for all combos needed.
ELF5<-cleanGRsPeaks[[1]]
values(ELF5)<-NULL
ER<-reduce(c(cleanGRsPeaks[[2]],cleanGRsPeaks[[3]]))
FOXA1<-reduce(c(cleanGRsPeaks[[4]],cleanGRsPeaks[[5]]))
H3K4me3<-reduce(c(cleanGRsPeaks[[6]],cleanGRsPeaks[[7]]))

cleanGRsPeaks=GRangesList()
cleanGRsPeaks[["ELF5"]]<-ELF5
cleanGRsPeaks[["FOXA1"]]<-FOXA1
cleanGRsPeaks[["ER"]]<-ER
cleanGRsPeaks[["H3K4me3"]]<-H3K4me3

peakSets<-GRangesList()
combinations<-combos(length(cleanGRsPeaks))$binary
for(i in 1:dim(combinations)[1]){
  idxs<-c(1:length(cleanGRsPeaks))[combinations[i,]*c(1:length(cleanGRsPeaks))]
  finalPeaks<-GRangesList()
  if(length(idxs)==1){
    cat(paste0(names(cleanGRsPeaks)[idxs]))
    cat("\t")
    cat(length(cleanGRsPeaks[[idxs]]))
    cat("\n")
    peakSets[[names(cleanGRsPeaks)[idxs]]]<-cleanGRsPeaks[[idxs]]
  } else {
    for(j in 1:c(length(idxs)-1)){
      if(j==1){
        finalPeaks[[1]]=cleanGRsPeaks[[idxs[j]]]
        names(finalPeaks)=names(cleanGRsPeaks)[idxs[j]]
      }
      #if the name of the j+1 sample is already present the do an Union, otherwise do an overlap
      oldName<-names(finalPeaks)[1]
      newName<-names(cleanGRsPeaks)[idxs[j+1]]
      #if(grepl(gsub("_.*","",newName),oldName)&length(finalPeaks[[1]])>0){
      #  finalPeaks[[1]]<-c(finalPeaks[[1]],cleanGRsPeaks[[idxs[j+1]]])
      #} else {
        gr<-finalPeaks[[1]]
        mat<-findOverlaps(gr,cleanGRsPeaks[[idxs[j+1]]])
        if(length(mat)>0){
          gr<-gr[unique(queryHits(mat))]
          finalPeaks[[1]]<-gr
        } else {
          finalPeaks[[1]]<-GRanges()
        }
      #}
      names(finalPeaks)<-paste0(names(finalPeaks)[1],"_",newName)
    }
    cat(paste0(names(finalPeaks)))
    cat("\t")
    cat(length(finalPeaks[[1]]))
    cat("\n")
    if(length(finalPeaks[[1]])>10){
      peakSets[[names(finalPeaks)]]<-finalPeaks[[1]]
    }
  }
}
temp=peakSets

peakL<-GRangesList()
for(i in 1:length(peakSets)){
  gr<-peakSets[[i]]
  names(gr)<-NULL
  names(gr)<-paste0(names(peakSets[i]),".",1:length(gr))
  peakL[[names(peakSets[i])]]<-gr
}
peakSets=peakL
save(peakSets,file="peakSets_noDiff.Rdata")
#second get the gene lists

peakAnno <- annotatePeak(unlist(peakSets), tssRegion=c(-2500, 2500), TxDb=txdb, annoDb="org.Hs.eg.db")
df=as.data.frame(peakAnno@anno)
dfL<-split(df,gsub("\\..*","",row.names(df)))

geneLists<-lapply(dfL,function(x){unique(x$ENSEMBL[abs(x$distanceToTSS)<=2500])})
geneLists=geneLists[sapply(geneLists,length)>=10]


#do the GSEA

#DE genes
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)

rank<-expression$t
names(rank)<-expression$ensembl_gene_id
rank=sort(rank)

fgseaRes <- fgsea(pathways = geneLists, 
                  stats = rank,
                  minSize=15,
                  nperm=100000)

#plot
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=11), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=4), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf(paste0(imageDir,"GSEA_peak_combinations_noDiff.pdf"),width=12,height=4)
plotGseaTable(geneLists[topPathways], rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#extract gene ids
fgseaResTab=as.data.frame(fgseaRes)

write.table(fgseaResTab[,1:7],"../project_results/tables/GSEA_noDiff.xls",row.names=F,quote=F,sep="\t")

#extract for all pathways

















