

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
library(ggrepel)


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
imageDir=paste0(resultsDir,"/figures/GSEA/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
system(paste0("mkdir -p ",imageDir))



pathway<-read.table("ER_associated.txt",header=T)
pathway=pathway[,1]
pathways<-list()
pathways[["ER_associated"]]<-pathway


pathway<-read.table("ER_associated_B.txt",header=T)
pathway=pathway[,1]
pathways[["ER_associated_B"]]<-pathway


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

for(gene in pathway){
  cat(gene)
  cat("\n")
  cat(names(rank)[grepl(gene,names(rank))])
  cat("\n")
}

pdf(paste0(imageDir,"ER_Alist_ER_Blist_GSEA.pdf"),width=12,height=4)
plotGseaTable(pathways[1:2], rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#extract gene ids
df=as.data.frame(fgseaRes)

write.table(df[,1:7],"../project_results/tables/GSEA_Alist_Blist.xls",row.names=F,quote=F,sep="\t")

#extract for all pathways

#leading edge

leGenes<-unlist(df[2,8])


pathway=pathways[["ER_associated_B"]]
stats=rank
gseaParam=1
rnk <- rank(-stats)
ord <- order(rnk)
statsAdj <- stats[ord]
statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
statsAdj <- statsAdj/max(abs(statsAdj))
pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
pathway <- sort(pathway)
gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
    returnAllExtremes = TRUE)
bottoms <- gseaRes$bottoms
tops <- gseaRes$tops
n <- length(statsAdj)
xs <- as.vector(rbind(pathway - 1, pathway))
ys <- as.vector(rbind(bottoms, tops))
toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
diff <- (max(tops) - min(bottoms))/8
x = y = NULL

leRanks<-data.frame(rank=rnk[names(rnk) %in% leGenes],genes=names(rnk)[names(rnk) %in% leGenes])

temp=toPlot
merged=merge(toPlot,leRanks,by.x="x",by.y="rank",all.x=T)
#merged$genes[is.na(merged$genes)]=""
toPlot=merged
dup=duplicated(toPlot$genes)
toPlot$genes[dup]=NA
pdf(paste0(imageDir,"ERassociated_GSEA_Blist_test.pdf"),width=12,height=6)
g <- ggplot(toPlot, aes(x = x, y = y)) + 
  geom_point(color = "green", size = 0.1) + 
  geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
  geom_hline(yintercept = 0, colour = "black") + 
  geom_line(color = "green") + 
  theme_bw() + 
  ylim(c(-0.5,0.6)) +
  geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
        y = -diff/2, xend = x, yend = diff/2), size = 0.2) + 
  theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "rank", y = "enrichment score") + 
  geom_text_repel(data=subset(toPlot,!is.na(genes)),
      angle        = 90,
      segment.size  = 0.5,
      segment.colour = "gray80",
      segment.alpha = 0.5,
       direction     = "y",
      ylim=c(max(tops),0.6),size=2,
#      xlim=(c(0,n+1)),
      aes(x = x, y = diff/2, 
        label = genes))
g

dev.off()

leadingEdge<-expression[expression$hgnc_symbol %in% pathway,]
leadingEdge=leadingEdge[order(leadingEdge$t),]
write.table(leadingEdge,"../project_results/tables/leadingEdge_ER_pathway.xls",row.names=F,quote=F,sep="\t")



pathway=pathways[["ER_associated_B"]]
pathwayGenes=pathway
stats=rank
gseaParam=1
rnk <- rank(-stats)
ord <- order(rnk)
statsAdj <- stats[ord]
statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
statsAdj <- statsAdj/max(abs(statsAdj))
pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
pathway <- sort(pathway)
gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
    returnAllExtremes = TRUE)
bottoms <- gseaRes$bottoms
tops <- gseaRes$tops
n <- length(statsAdj)
xs <- as.vector(rbind(pathway - 1, pathway))
ys <- as.vector(rbind(bottoms, tops))
toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
diff <- (max(tops) - min(bottoms))/8
x = y = NULL

leRanks<-data.frame(rank=rnk[names(rnk) %in% pathwayGenes],genes=names(rnk)[names(rnk) %in% pathwayGenes])

temp=toPlot
merged=merge(toPlot,leRanks,by.x="x",by.y="rank",all.x=T)
#merged$genes[is.na(merged$genes)]=""
toPlot=merged
dup=duplicated(toPlot$genes)
toPlot$genes[dup]=NA
pdf(paste0(imageDir,"ERassociated_GSEA_Blist_all_B.pdf"),width=12,height=6)
g <- ggplot(toPlot, aes(x = x, y = y)) + 
  geom_point(color = "green", size = 0.1) + 
  geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
  geom_hline(yintercept = 0, colour = "black") + 
  geom_line(color = "green") + 
  theme_bw() + 
  ylim(c(-0.5,0.6)) +
  geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
        y = -diff/2, xend = x, yend = diff/2), size = 0.2) + 
  theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "rank", y = "enrichment score") + 
  geom_text_repel(data=subset(toPlot,!is.na(genes)),
      angle        = 90,
      segment.size  = 0.5,
      segment.colour = "gray80",
      segment.alpha = 0.5,
       direction     = "y",

      ylim=c(max(tops),0.6),size=2,
#      xlim=(c(0,n+1)),
      aes(x = x, y = diff/2, 
        label = genes))
g

dev.off()

leadingEdge<-expression[expression$hgnc_symbol %in% pathwayGenes,]
leadingEdge=leadingEdge[order(leadingEdge$t),]
write.table(leadingEdge,"../project_results/tables/leadingEdge_ER_pathway.xls",row.names=F,quote=F,sep="\t")
write.table(leRanks,"../project_results/tables/leadingEdge_ER_pathway_ranks.xls",row.names=F,quote=F,sep="\t")








