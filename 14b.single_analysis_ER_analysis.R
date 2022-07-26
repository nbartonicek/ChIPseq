

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

#perform the analysis for each ChIP, normally
#sample peak

ELF5=cleanGRsPeaks[["ELF5Dox"]]
sampleName="FOXA1"
imageDir=paste0(resultsDir,"/figures/",sampleName,"/")
system(paste0("mkdir -p ",imageDir))

cleanGRsPeaks=cleanGRsPeaks[grepl(sampleName,names(cleanGRsPeaks))]
sapply(cleanGRsPeaks,length)
sum(countOverlaps(cleanGRsPeaks[[2]],cleanGRsPeaks[[1]])>0)


#make a venn of ELF5+/-
ovlp = findOverlaps( cleanGRsPeaks[[1]], cleanGRsPeaks[[2]] )
ov = min( length(unique( queryHits(ovlp) )), length(unique( subjectHits(ovlp) ) ) )
pdf(paste0(imageDir,"venn_",paste(names(cleanGRsPeaks)[1], names(cleanGRsPeaks)[2],sep="_"),".pdf"),width=4,height=4)
draw.pairwise.venn( 
 area1=length(cleanGRsPeaks[[1]]),
 area2=length(cleanGRsPeaks[[2]]), 
 cross.area=ov, 
 category=c(names(cleanGRsPeaks)[1], names(cleanGRsPeaks)[2]), 
 fill=c("red", "blue"), ext.text=F,fontfamily="sans",cat.fontfamily="sans",cat.pos=180,
 cat.cex=0.7)
dev.off()


sum(countOverlaps(cleanGRsPeaks[[1]],ELF5)>0)
sum(countOverlaps(cleanGRsPeaks[[2]],ELF5)>0)


for(sampleName in names(cleanGRsPeaks)[1:2]){
  cat(sampleName)
  cat("\n")
#summits=cleanGRs[[sampleName]]
  peak=cleanGRsPeaks[[sampleName]]
  imageDir=paste0(resultsDir,"/figures/",sampleName,"/")
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
    tagMatrix <- getTagMatrix(peak, windows=promoterSubset)
    if(!is.null(dim(tagMatrix))){
      if(dim(tagMatrix)[1]>0){
        pdf(paste0(imageDir,geneType,"_average_profile.pdf"),width=4,height=4)
        p<-plotAvgProf(tagMatrix, xlim=c(-10000, 10000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
        print(p)
        dev.off()

        #heatmap
        pdf(paste0(imageDir,geneType,"_heatmap.pdf"),width=3,height=8)
        tagHeatmap(tagMatrix, xlim=c(-10000, 10000), color="red")
        dev.off()
      }
    }
    #Chromosome regions
    if(geneType!="all"){
      mat<-findOverlaps(peak,promoterSubset)
      peakSubset<-peak[unique(queryHits(mat))]
    } else {peakSubset=peak}
    if(length(peakSubset)>0){
      peakAnno <- annotatePeak(peakSubset, tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Hs.eg.db")
      sum( peakAnno@annoStat$Frequency[c(1:3)])
      sum( peakAnno@annoStat[peakAnno@annoStat$Feature %in% c("1st Intron","Other Intron","Distal Intergenic"),"Frequency"])


      pdf(paste0(imageDir,geneType,"_annotation_pie.pdf"),width=8,height=4)
      p<-plotAnnoPie(peakAnno)
      print(p)
      dev.off()
      pdf(paste0(imageDir,geneType,"_annotation_bar.pdf"),width=8,height=5)
      p<-plotAnnoBar(peakAnno)
      print(p)
      dev.off()
      pdf(paste0(imageDir,geneType,"_annotation_upset.pdf"),width=8,height=4)
      upsetplot(peakAnno, vennpie=TRUE)
      print(p)
      dev.off()
      if(length(peakAnno)>1){
        #Distance to TSS
        pdf(paste0(imageDir,geneType,"_dist2TSS.pdf"),width=6,height=2)
        p<-plotDistToTSS(peakAnno,title="Distribution of transcription factor-binding loci\nrelative to TSS")
        print(p)
        dev.off()

        #Pathway analysis reactome/KEGG/interactions
        pathway <- enrichPathway(as.data.frame(peakAnno)$geneId)
        pdf(paste0(imageDir,geneType,"_reactome_pathway.pdf"),width=8,height=4)
        p<-dotplot(pathway)
        print(p)
        dev.off()
      }
    }
  }
}
  #pathway analysis great
  gr.hg38=peakSubset
  ch = import.chain(paste0(projectDir,"/annotation/hg38ToHg19.over.chain"))
  seqlevelsStyle(gr.hg38) = "UCSC"  # necessary
  gr.hg19 = liftOver(gr.hg38, ch)
  gr.hg19  = unlist(gr.hg19 )
  genome(gr.hg19) = "hg19"
  gr.hg19$score=abs(as.integer(gr.hg19$score))
  export(gr.hg19,paste0(projectDir,"/project_results/beds/",sampleName,"_hg19.bed"))

  bed=data.frame(chr=as.character(seqnames(peakSubset)),start=start(peakSubset),end=end(peakSubset))
  #bed=bed[bed$chr=="chr1",]
  job = submitGreatJob(bed)
  tb = getEnrichmentTables(job, ontology = c("GO Biological Process","MSigDB Oncogenic Signatures","MSigDB Perturbation","MSigDB Predicted Promoter Motifs"))
  for(pathwayName in names(tb)){
    pathwayOutName=gsub(" ","_",pathwayName)
    write.table(tb[[pathwayName]][1:20,],paste0(imageDir,geneType,"_",pathwayOutName,".txt"))
  }
#
  pdf(paste0(imageDir,geneType,"_great.pdf"),width=8,height=4)
  par(mfrow = c(1, 3))
  res = plotRegionGeneAssociationGraphs(job)
  dev.off()

#}


#redo analysis for the thesis
 kb3=promoters(promoterSubset,upstream=3000,downstream=3000)
peaks<-cleanGRsPeaks[[sampleName]]
sum(countOverlaps(peaks,kb3)>0)

#for enriched/depleted
gr=unlist(cleanGRs[c("FOXA1Dox","FOXA1NoDox")])
mat<-findOverlaps(gr,peaks)
gr=reduce(gr[unique(queryHits(mat))])
#export 500 and 1kb around the summits, submit to meme/dreme
#gr=cleanGRs[[sampleName]]
gr500<-resize(gr,fix="center",width=500)
gr1kb<-resize(gr,fix="center",width=1000)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/",sampleName,"_500bp.bed"))


seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
names(seq) = paste0("SEQUENCE_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/",sampleName,"_1kb.bed"))


for(deName in names(deGenesL)){
  deGenes<-deGenesL[[deName]]
  deGenes<-deGenes[!is.na(deGenes)]
  gns<-bitr(deGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  promoterSubset=promoter[promoter$gene_id %in% gns$ENTREZID]
  cat(sum(countOverlaps(peakSubset,promoterSubset)>0))
  cat("\n")
  cat(sum(countOverlaps(peakSubset,promoterSubset)>0)/length(promoterSubset))
  cat("\n")
}


