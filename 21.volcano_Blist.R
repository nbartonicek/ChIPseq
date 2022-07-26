
library(ShortRead)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
#library(Go.db)
library(ggrepel)
timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "human"
library(DOSE)



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


expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
pathway<-read.table("ER_associated_genes.txt",header=T)


expression$genelabels <- ""
expression$genelabels[expression$hgnc_symbol %in% pathway[,1]] = expression$hgnc_symbol[expression$hgnc_symbol %in% pathway[,1]]
#expression$genelabels[expression$adj.P.Val>0.1]<-""
expression$color<-"gray"
expression$color[expression$hgnc_symbol %in% pathway[,1]]<-"darkred"
expression$color[(expression$hgnc_symbol %in% pathway[,1]) & (expression$adj.P.Val>0.1)]<-"black"

pdf(paste0(imageDir,"volcano_ER_B_list.pdf"),width=8,height=8)
  p<-ggplot(expression) +
    geom_point(aes(x=logFC, y=-log10(adj.P.Val)),colour=expression$color) +
    ggtitle("ER associated") +
    xlab("log2 fold change") + xlim(c(-6,6)) +
    ylab("-log10 adjusted p-value") +
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
    geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val), label = genelabels)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) 
  print(p)
dev.off()


expression$genelabels <- ""
expression$genelabels[expression$hgnc_symbol %in% pathway[,1]] = expression$hgnc_symbol[expression$hgnc_symbol %in% pathway[,1]]
expression$genelabels[expression$adj.P.Val>0.1]<-""
expression$color<-"gray70"
expression$color[expression$hgnc_symbol %in% pathway[,1]]<-"darkred"
expression$color[(expression$hgnc_symbol %in% pathway[,1]) & (expression$adj.P.Val>0.1)]<-"black"
#expression=expression[1:100,]
pdf(paste0(imageDir,"volcano_ER_B_list_selected2.pdf"),width=8,height=8)
  p<-ggplot(expression) +
    geom_point(aes(x=logFC, y=-log10(adj.P.Val)),colour=expression$color) +
    ggtitle("ER associated") +
    xlab("log2 fold change") + xlim(c(-6,6)) +
    ylab("-log10 adjusted p-value") +
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
    geom_text_repel(data=subset(expression,logFC<0),
      nudge_x = -6 - subset(expression,logFC<0)$logFC, 
      segment.size  = 0.6,
      segment.colour = "gray40",
#      direction     = "y",
#      hjust         = 1, xlim=c(-6,6),
      aes(x = logFC, y = -log10(adj.P.Val), 
        label = genelabels)) +
    geom_text_repel(data=subset(expression,logFC>0),
      nudge_x = 6 - subset(expression,logFC>0)$logFC,       
      segment.size  = 0.6,
#      direction     = "y",
      segment.colour = "gray40",

#      hjust         = 1, xlim=c(-6,6),
      aes(x = logFC, y = -log10(adj.P.Val), 
        label = genelabels)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) 
  print(p)
dev.off()



expression$genelabels <- ""
expression$genelabels[expression$hgnc_symbol %in% pathway[,1]] = expression$hgnc_symbol[expression$hgnc_symbol %in% pathway[,1]]
expression$genelabels[expression$adj.P.Val>0.1]<-""
expression$color<-"gray70"
expression$color[expression$hgnc_symbol %in% pathway[,1]]<-"darkred"
expression$color[(expression$hgnc_symbol %in% pathway[,1]) & (expression$adj.P.Val>0.1)]<-"black"
#expression=expression[1:100,]
pdf(paste0(imageDir,"volcano_ER_B_list_selected_2rows.pdf"),width=8,height=8)
  p<-ggplot(expression) +
    geom_point(aes(x=logFC, y=-log10(adj.P.Val)),colour=expression$color) +
    ggtitle("ER associated") +
    xlab("log2 fold change") + xlim(c(-6,6)) +
    ylab("-log10 adjusted p-value") +
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) +
    geom_text_repel(data=subset(expression,logFC<(-0.273)),
      nudge_x = -6 - subset(expression,logFC<(-0.273))$logFC, 
#      segment.size  = 0.6,
      segment.colour = "gray40",
      direction     = "y",
      xlim=c(-6,6),ylim = c(1,10),
      aes(x = logFC, y = -log10(adj.P.Val), 
        label = genelabels)) +
    geom_text_repel(data=subset(expression,logFC>-0.273&logFC<0),
      nudge_x = -4 - subset(expression,logFC>-0.273&logFC<0)$logFC, 
#      segment.size  = 0.6,
      segment.colour = "gray40",
      direction     = "y",
      xlim=c(-4,6), ylim = c(1,10),
      aes(x = logFC, y = -log10(adj.P.Val), 
        label = genelabels)) +
    geom_text_repel(data=subset(expression,logFC<0.28&logFC>0),
      nudge_x = 4 - subset(expression,logFC<0.28&logFC>0)$logFC,       
#      segment.size  = 0.6,
      direction     = "y",
      segment.colour = "gray40",

      xlim=c(-6,4),ylim = c(1,10),
      aes(x = logFC, y = -log10(adj.P.Val), 
        label = genelabels)) +
    geom_text_repel(data=subset(expression,logFC>0.28),
      nudge_x = 6 - subset(expression,logFC>0.28)$logFC,       
#      segment.size  = 0.6,
      direction     = "y",
      segment.colour = "gray40",

      xlim=c(-6,6),ylim = c(1,10),
      aes(x = logFC, y = -log10(adj.P.Val), 
        label = genelabels)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) 
  print(p)
dev.off()

