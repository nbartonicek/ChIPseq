


library(GenomicRanges)
#library(ChIPpeakAnno)
library(ggplot2)
library(rtracklayer)
#library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(clusterProfiler)
#library(ReactomePA)
#library(rGREAT)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(hier.part)
library(fgsea)


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


geneFiles<-list.files(".",pattern="outcome")
pathways<-list()

for(geneFile in geneFiles){
  listName=gsub(".txt","",geneFile)
  pathway<-read.table(geneFile,header=T)
  pathway=pathway[,1]
  pathways[[listName]]<-pathway

}

#DE genes
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
geneid <- mapIds(org.Hs.eg.db, expression$ENSEMBL_ID, "SYMBOL","ENSEMBL")
expression$symbol=geneid
aggregated<-aggregate(expression$t,list(as.character(expression$symbol)),mean)

rank<-aggregated$x
names(rank)<-as.character(aggregated$Group.1)
rank=sort(rank)

fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=10, maxSize=500,
                  nperm=100000)
#plot
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=1), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=9), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf(paste0(imageDir,"MDSC_ER_GSEA.pdf"),width=12,height=4)
plotGseaTable(pathways[1:2], rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#extract gene ids
fgseaRes=as.data.frame(fgseaRes)

#extract for all pathways



#DE genes
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
geneid <- mapIds(org.Hs.eg.db, expression$ENSEMBL_ID, "SYMBOL","ENSEMBL")
expression$symbol=geneid
aggregated<-aggregate(expression$t,list(as.character(expression$symbol)),mean)

rank<-aggregated$x
names(rank)<-as.character(aggregated$Group.1)
rank=sort(rank)

pathways<-list()
pathways[[1]]<-pathway
fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=10, maxSize=500,
                  nperm=100000)

for(gene in pathway){
  cat(gene)
  cat("\n")
  cat(names(rank)[grepl(gene,names(rank))])
  cat("\n")
}

#plot
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=1), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=9), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

pdf(paste0(imageDir,"ERassociated_GSEA.pdf"),width=12,height=4)
plotGseaTable(pathways, rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#extract gene ids
fgseaRes=as.data.frame(fgseaRes)

write.table(fgseaRes[,1:7],"../project_results/tables/GSEA.xls",row.names=F,quote=F,sep="\t")

#extract for all pathways
















