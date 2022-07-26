

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
tableDir<-paste0(resultsDir,"/tables/")
system(paste0("mkdir -p ",cleanRobjectsDir))
system(paste0("mkdir -p ",tableDir))

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

#DE genes
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]
deGenesAll<-unique(as.character(expression$hgnc_symbol[(abs(expression$logFC)>=log2(1.5)) & (expression$adj.P.Val<=0.05)]))
deGenesPos<-unique(as.character(expression$hgnc_symbol[expression$logFC>=log2(1.5) & expression$adj.P.Val<=0.05]))
deGenesNeg<-unique(as.character(expression$hgnc_symbol[expression$logFC<=-log2(1.5) & expression$adj.P.Val<=0.05]))

deGenesL<-list()
deGenesL[["DE"]]<-deGenesAll
deGenesL[["pos"]]<-deGenesPos
deGenesL[["neg"]]<-deGenesNeg

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")), upstream=10001, downstream=10000)
#promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)



#For each "diff" extract the list of proximal genes and then do the analysis
diff_enriched_samples<-names(cleanGRsPeaks)[grepl("diff.*enriched",names(cleanGRsPeaks))]

geneLists<-list()

for(sampleName in diff_enriched_samples){
  cat(sampleName)
  cat("\n")
#summits=cleanGRs[[sampleName]]
  peak=cleanGRsPeaks[[sampleName]]
  imageDir=paste0(resultsDir,"/figures/",sampleName,"/diff/")
  system(paste0("mkdir -p ",imageDir))
   

  for(geneType in c("all",names(deGenesL))){
    cat(geneType)
    cat("\n")
    if(geneType!="all"){
      deGenes<-deGenesL[[geneType]]
      deGenes<-deGenes[!is.na(deGenes)]
      gns<-bitr(deGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
      promoterSubset=promoter[promoter$gene_id %in% gns$ENTREZID]
      
    } else {promoterSubset=promoter}
    values(promoterSubset)=NULL
    names(promoterSubset)=1:length(promoterSubset)
    #tagMatrix <- getTagMatrix(peak, windows=promoterSubset)
    
    #Chromosome regions
    if(geneType!="all"){
      mat<-findOverlaps(peak,promoterSubset)
      peakSubset<-peak[unique(queryHits(mat))]
    } else {peakSubset=peak}
    if(length(peakSubset)>0){
      peakAnno <- annotatePeak(peakSubset, tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
      
      genes<-peakAnno@anno$SYMBOL[!grepl("Distal",peakAnno@anno$annotation)]
      geneLists[[paste0(sampleName,"_",geneType)]]<-unique(genes)

    }
  }
}

for(sampleName in names(geneLists)){
  outFile=paste0(tableDir,sampleName,".tab")
  write(geneLists[[sampleName]],outFile)
}

