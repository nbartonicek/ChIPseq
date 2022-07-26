
library(GenomicRanges)
library(ShortRead)
#library(R.utils)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
#library(RUVSeq)
library(org.Hs.eg.db)
library(DESeq)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPpeakAnno)
library(clusterProfiler)
library(pheatmap)

timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "human"
ensVer <- 84

library(fgsea)


homedir="../../../"



######## directory structure #######
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/API/")
projectname="JNCI"
robjectsDir = paste(resultsDir,"/",projectname,".Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
#names(chrs)=gsub("chr","",names(chrs))
#names(chrs)[25]="MT"
gr<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))

#Load both type and class lists
#load(paste0(annotationDir,"class.Rdata"))
#grClass<-class
#load(paste0(annotationDir,"type.Rdata"))
#grType<-type



#########################################
########## 1. Load in the files
#########################################

df<-read.table(paste0(resultsDir,"/JNCI/brainMetPairs.salmon.cts.txt"))

annotation<-read.table("jnci_annotation.txt",sep="\t",header=T)

#eliminate extra samples and reorder
df<-df[,!c(colnames(df) %in% c("X7P_RCS","X7M_RCS"))]
df<-df[,c(grep("P_",colnames(df)),grep("M_",colnames(df)))]


#trying out mds
group <- factor(paste0(c(annotation$ER,annotation$ER.1)))
y <- DGEList(df,group=group)
y <- calcNormFactors(y)

# without labels, indexes of samples are plotted.
col <- as.numeric(group)
pdf("../project_results/figures/samples_MDS_JNCI.pdf",width=12,height=8)
mds <- plotMDS(y, top=200, col=col,labels=c(paste0("P_",annotation$Case),paste0("M_",annotation$Case)))
dev.off()

#first group
ER="Pos"
Endo="Yes"
type="P"
group1<-paste0("ER_",ER,".","Endo_",".",Endo,"Type_",type)

#second group
ER.1="Pos"
Endo.1="Yes"
type.1="M"
group2<-paste0("ER_",ER.1,".","Endo_",".",Endo.1,"Type_",type.1)

index1<-1:42*c((annotation$ER==ER)&(annotation$Endo==Endo),rep(FALSE,21))
index1<-index1[!is.na(index1)]
index1<-index1[index1>0]
df1<-df[,index1]
#index2<-1:42*c(rep(FALSE,21),c((annotation$ER.1==ER.1)&(annotation$Endo.1==Endo.1)))
#index2<-index2[!is.na(index2)]
#index2<-index2[index2>0]
#df2<-df[,index2]

index2<-21+index1
df2<-df[,index2]

#for non ER:
index1<-c(1:21)[-index1]
index2<-21+index1
df1<-df[,index1]
df2<-df[,index2]


dfShort<-cbind(df1,df2)


nSamples1<-dim(df1)[2]
nSamples2<-dim(df2)[2]
nSamples<-nSamples1+nSamples2
#Annotate with entrez, aggregate
egENSEMBL <- toTable(org.Hs.egENSEMBL)
row.names(dfShort)=gsub("\\..*","",row.names(dfShort))
m <- match(row.names(dfShort), egENSEMBL$ensembl_id)
dfShort$EntrezGene<-egENSEMBL$gene_id[m]

egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(dfShort$EntrezGene, egSYMBOL$gene_id)
dfShort$symbol<-egSYMBOL$symbol[m]

#eliminate duplicated symbols
o <- order(rowSums(dfShort[,1:nSamples]), decreasing=TRUE)
dfShort=dfShort[o,]
d<-duplicated(dfShort$symbol)
dfShort<-dfShort[!d,]

#add the gene length
#data$gene_id=gsub("\\..*","",data$gene_id)
#df_merged<-merge(df,data[,c(1,4)],by.x=0,by.y="gene_id")
#df<-df_merged
#row.names(df)<-df[,1]
#df<-df[,-1]
#eliminate lowly expressed
include<-apply(dfShort[,1:nSamples],1,function(x){sum(x>=10)>=2})
dfShort<-dfShort[include,]
dfShort=dfShort[!is.na(dfShort$symbol),]
row.names(dfShort)<-dfShort$symbol

#change the id to entrez
#o <- order(rowSums(df[,1:6]), decreasing=TRUE)
#df=df[o,]
#d<-duplicated(df$EntrezGene)
#df<-df[!d,]
#row.names(df)<-df$EntrezGene


#provide the EC numbers
#egENZYME <- toTable(org.Mm.egENZYME)
#m <- match(row.names(dfShort), egENZYME$gene_id)
#
#dfShort[,"EC number"]<-egENZYME$ec_number[m]


treatment <- rep(c(type.1, type), times=c(nSamples2,nSamples1))
patient <- gsub("\\w_.*","",colnames(dfShort))[1:nSamples]

design <- model.matrix(~patient+treatment)

#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=dfShort[1:nSamples])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmQLFit(expr,design)
lrt <- glmQLFTest(fit)


# without labels, indexes of samples are plotted.
col <- as.numeric(factor(treatment))
pdf("../project_results/figures/samples_MDS_blocked_JNCI.pdf",width=12,height=8)
mds <- plotMDS(expr, top=200, col=col,labels=patient)
dev.off()

a=dfShort[,1:14]
a=a[row.names(a) %in% c("ELF5","ESR1"),]
a$genes<-row.names(a)




dataM<-melt(a)
dataM$loc<-rep(c("P","M"),each=14)
colnames(dataM)<-c("symbol","patient","counts","loc")
dataM$symbol<-factor(dataM$symbol,levels=c("ESR1","ELF5"))
dataM$patient<-gsub("(P|M)_","_",dataM$patient)
dataM$loc<-factor(dataM$loc,levels=c("P","M"))
samples<-colnames(fittedValues_genes)[1:7]
samples<-gsub("(P|M)_","_",samples)
dataM$patient<-factor(dataM$patient,levels=samples)
pdf("../project_results/figures/ELF5_ESR1_JNCI_raw.pdf",width=12,height=4)
p<-ggplot(dataM,aes(patient,counts,fill=loc))
p<-p+geom_bar(stat="identity",position="dodge")
p<-p+facet_wrap(~symbol,scales="free_y")
p<-p+ theme(axis.text.x = element_text(angle = 90))
p
dev.off()

samples<-dfShort[1:nSamples]

meta<-data.frame(row.names=colnames(samples),treatment=rev(treatment),patient=patient)

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData=samples, colData=meta, design=~patient+treatment)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


fittedValues<-normalized_counts
fittedValues_genes<-fittedValues[row.names(fittedValues) %in% c("ELF5","ESR1"),]
dataM<-melt(fittedValues_genes)
dataM$loc<-rep(c("P","M"),each=14)
colnames(dataM)<-c("symbol","patient","counts","loc")
dataM$patient<-gsub("(P|M)_","_",dataM$patient)
dataM$loc<-factor(dataM$loc,levels=c("P","M"))
samples<-colnames(fittedValues_genes)[1:7]
samples<-gsub("(P|M)_","_",samples)
dataM$patient<-factor(dataM$patient,levels=samples)
pdf("../project_results/figures/ELF5_ESR1_JNCI_DESeq2_normalised.pdf",width=12,height=4)
p<-ggplot(dataM,aes(patient,counts,fill=loc))
p<-p+geom_bar(stat="identity",position="dodge")
p<-p+facet_wrap(~symbol,scales="free_y")
p<-p+ theme(axis.text.x = element_text(angle = 90))
p
dev.off()




fittedValues<-fitted.values(lrt)
fittedValues_genes<-fittedValues[row.names(fittedValues) %in% c("ELF5","ESR1"),]
dataM<-melt(fittedValues_genes)
dataM$loc<-rep(c("P","M"),each=14)
colnames(dataM)<-c("symbol","patient","counts","loc")
dataM$patient<-gsub("(P|M)_","_",dataM$patient)
dataM$loc<-factor(dataM$loc,levels=c("P","M"))
samples<-colnames(fittedValues_genes)[1:7]
samples<-gsub("(P|M)_","_",samples)
dataM$patient<-factor(dataM$patient,levels=samples)
pdf("../project_results/figures/ELF5_ESR1_JNCI_normalised.pdf",width=12,height=4)
p<-ggplot(dataM,aes(patient,counts,fill=loc))
p<-p+geom_bar(stat="identity",position="dodge")
p<-p+facet_wrap(~symbol,scales="free_y")
p<-p+ theme(axis.text.x = element_text(angle = 90))
p
dev.off()



logFCs<-as.data.frame(topTags(lrt,n=dim(dfShort)[1]))
plus<-logFCs[logFCs$logFC>0,]
minus<-logFCs[logFCs$logFC<0,]
significantIDs<-c(row.names(plus[1:25,]),row.names(minus[1:25,]))
fittedValues<-fitted.values(lrt)
fittedValues50<-fittedValues[row.names(fittedValues) %in% significantIDs,]
colnames(fittedValues50)<-gsub("_R1","",colnames(fittedValues50))

#Organise



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

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
load(paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))

sampleName<-"ELF5Dox"
peak=cleanGRsPeaks[[sampleName]]
promoterOverlap<-sum(countOverlaps(peak,promoter)>0)

pathways<-list()
for(geneType in names(deGenesL)){
  deGenes<-deGenesL[[geneType]]
  deGenes<-deGenes[!is.na(deGenes)]
  gns<-bitr(deGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  promoterSubset=promoter[promoter$gene_id %in% gns$ENTREZID]
  mat<-findOverlaps(peak,promoterSubset)
  promoterSubset<-promoterSubset[unique(subjectHits(mat))]

  DEgenes_with_ChIP<-bitr(promoterSubset$gene_id, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")$SYMBOL
  pathways[[geneType]]<-deGenes
  pathways[[paste0(geneType,"_ChIP")]]<-DEgenes_with_ChIP
}



logFCs$rank<-ifelse(logFCs$logFC>=0,logFCs$F,-logFCs$F)

rank<-logFCs$rank
names(rank)<-row.names(logFCs)
rank=sort(rank)

fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=10, maxSize=500,
                  nperm=100000)

pdf(paste0(imageDir,"/JNCI_GSEA_byPatient.pdf"),width=12,height=8)
plotGseaTable(pathways, rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#do ranks 7 independent pairs
#do GSEA for gene and ChIP enrichment, as well as list B and Tamoxifen targets
#then plot genes




########## Normalize means by DESeq normalization factors, and dispersion with quantile normalization
library(DESeq)
countTable=dfShort[,1:14]
patient=as.factor(rep(c(1:7), 2))
treatment=rep(c("P","M"),each=7)
treatment[c(1,8)]<-c("4P","4M")
treatment<-as.factor(treatment)
short<-countTable

cds = newCountDataSet( short, treatment )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds, method="blind", sharingMode="fit-only" )
res = nbinomTest( cds, "4M", "4P" )
save(res,file="JINC.Rdata")
results[[i]]<-res

#load results for individual DESeq: 
load("../project_results/Robjects/DESeq_JNCI_nonER.Rdata")

#

#DE genes for T47D
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]




rank<-expression$t
names(rank)<-expression$hgnc_symbol
rank=sort(rank)
d<-duplicated(names(rank))
rank<-rank[!d]


pathways<-list()
for(i in 1:length(results)){
  sampleName<-colnames(dfShort)[i]
  sampleName<-gsub("P_.*","",sampleName)
  padj<-results[[i]]$padj
  padj[is.na(padj)]<-1
  geneList<-results[[1]]$id[padj<0.2]
  geneList<-geneList[geneList %in% names(rank)]
  pathways[[sampleName]]<-geneList
}

fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=10, maxSize=1000,
                  nperm=100000)

pdf(paste0(imageDir,"/JNCI_ERneg_GSEA.pdf"),width=12,height=4)
plotGseaTable(pathways, rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()



normalized.count <- counts(cds, normalized = TRUE)
save(normalized.count,file=paste0(robjectsDir,"normalized.count.Rdata"))
normalized.count<-normalizeQuantiles(normalized.count)



######## start with patient list




rank<-ifelse(res$log2FoldChange>0,-log1p(res$pval),log1p(res$pval))
rank<-res$log2FoldChange

names(rank)<-res$id
rank=sort(rank)
d<-duplicated(names(rank))
rank<-rank[!d]

rank[is.na(rank)]=0
rank[rank=="Inf"]=max(rank[rank!="Inf"],na.rm=T)
rank[rank=="-Inf"]=min(rank[rank!="-Inf"],na.rm=T)

tam<-read.table("../project_results/tables/FOXA1_doxEnriched_tamoxifenResistant_genes.xls",header=T,sep="\t")
genes<-c(unique(as.character(tam$SYMBOL)))
pathways<-list()
pathways[["TAM"]]<-genes
ER<-read.table("ER_associated.txt",header=T)[,1]
pathways[["ER"]]<-ER
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]
deGenesAll<-unique(as.character(expression$hgnc_symbol[(abs(expression$logFC)>=log2(1.5)) & (expression$adj.P.Val<=0.05)]))
deGenesPos<-unique(as.character(expression$hgnc_symbol[expression$logFC>=log2(1.5) & expression$adj.P.Val<=0.05]))
deGenesNeg<-unique(as.character(expression$hgnc_symbol[expression$logFC<=-log2(1.5) & expression$adj.P.Val<=0.05]))
pathways[["ELF5deGenesAll"]]<-deGenesAll
pathways[["ELF5deGenesPos"]]<-deGenesPos
pathways[["ELF5deGenesNeg"]]<-deGenesNeg
pathways[["CHIP_T47D_MCF7"]]<-genesChip
fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=10, maxSize=500,
                  nperm=100000)



pdf(paste0(imageDir,"/T47D_DEChip.pdf"),width=12,height=4)
plotGseaTable(pathways, rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()

#plot ER+ and ER- genes 
#load in listB and tamoxifen genes
#combine all the data
#pheatmap change


load("../project_results/Robjects/DESeq_JNCI.Rdata")
log2FCsER<-lapply(results,function(x){x$log2FoldChange})
log2FCsER_df<-do.call("cbind",log2FCsER)
row.names(log2FCsER_df)=results[[1]]$id



#first group
ER="Pos"
Endo="Yes"
type="P"
group1<-paste0("ER_",ER,".","Endo_",".",Endo,"Type_",type)
index1<-1:42*c((annotation$ER==ER)&(annotation$Endo==Endo),rep(FALSE,21))
index1<-index1[!is.na(index1)]
index1<-index1[index1>0]
df1<-df[,index1]

colnames(log2FCsER_df)<-names(df1)


load("../project_results/Robjects/DESeq_JNCI_nonER.Rdata")
log2FCsnonER<-lapply(results,function(x){x$log2FoldChange})
log2FCsnonER_df<-do.call("cbind",log2FCsnonER)
row.names(log2FCsnonER_df)=results[[1]]$id



#first group
ER="Pos"
Endo="Yes"
type="P"
group1<-paste0("ER_",ER,".","Endo_",".",Endo,"Type_",type)
index1<-1:42*c((annotation$ER==ER)&(annotation$Endo==Endo),rep(FALSE,21))
index1<-index1[!is.na(index1)]
index1<-index1[index1>0]
df1<-df[,c(1:21)[-index1]]

colnames(log2FCsnonER_df)<-names(df1)



#now select genes
tam<-read.table("../project_results/tables/FOXA1_doxEnriched_tamoxifenResistant_genes.xls",header=T,sep="\t")
genes<-c("ELF5","ESR1",unique(as.character(tam$SYMBOL)))

er<-log2FCsER_df[row.names(log2FCsER_df) %in% genes,]
noner<-log2FCsnonER_df[row.names(log2FCsnonER_df) %in% row.names(er),]

total<-cbind(er,noner)

write.table(total,file="../project_results/tables/JNCI_tamResistant_DE.xls",quote=F,sep="\t")

total[is.na(total)]=0
total[total=="Inf"]=max(total[total!="Inf"],na.rm=T)
total[total=="-Inf"]=min(total[total!="-Inf"],na.rm=T)

pdf(paste0(imageDir,"/JNCI_tamResitant_DE.heatmap.pdf"),width=10,height=12)
pheatmap(total)
dev.off()

ER<-read.table("ER_associated.txt",header=T)
genes<-c("ELF5","ESR1",unique(as.character(ER[,1])))



er<-log2FCsER_df[row.names(log2FCsER_df) %in% genes,]
noner<-log2FCsnonER_df[row.names(log2FCsnonER_df) %in% row.names(er),]

total<-cbind(er,noner)

write.table(total,file="../project_results/tables/JNCI_listA_DE.xls",quote=F,sep="\t")

total[is.na(total)]=0
total[total=="Inf"]=max(total[total!="Inf"],na.rm=T)
total[total=="-Inf"]=min(total[total!="-Inf"],na.rm=T)

pdf(paste0(imageDir,"/JNCI_listA_DE.heatmap.pdf"),width=10,height=12)
pheatmap(total)
dev.off()



#DE genes for T47D
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]



T47D<-read.table("../project_results/GSE30407/GSE30407_DE.csv",sep=",",header=T)
T47D=unique(T47D[,1])
rank<-expression$t
names(rank)<-expression$hgnc_symbol
rank=sort(rank)
d<-duplicated(names(rank))
rank<-rank[!d]

T47D=T47D[T47D %in% names(rank)]

pathways<-list()
pathways[["T47D"]]<-T47D
fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=10, maxSize=500,
                  nperm=100000)

pdf(paste0(imageDir,"/T47D_DEChip.pdf"),width=12,height=4)
plotGseaTable(pathways, rank, fgseaRes, 
              gseaParam = 0.5)
dev.off()




#make res with only samples 1,3,4,5,7

library(DESeq)
countTable=dfShort[,c(1,4,7,8,11,14)]
patient=as.factor(rep(c(1:7), 2))[c(1,4,7,8,11,14)]
treatment=rep(c("P","M"),each=7)[c(1,4,7,8,11,14)]
treatment[c(1,8)]<-c("4P","4M")
treatment<-as.factor(treatment)
short<-countTable

cds = newCountDataSet( short, treatment )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds, method="blind", sharingMode="fit-only" )
res = nbinomTest( cds, "4M", "4P" )




rank<-ifelse(res$log2FoldChange>0,-log1p(res$pval),log1p(res$pval))
rank<-res$log2FoldChange

names(rank)<-res$id
rank=sort(rank)
d<-duplicated(names(rank))
rank<-rank[!d]

rank[is.na(rank)]=0
rank[rank=="Inf"]=max(rank[rank!="Inf"],na.rm=T)
rank[rank=="-Inf"]=min(rank[rank!="-Inf"],na.rm=T)

tam<-read.table("../project_results/tables/FOXA1_doxEnriched_tamoxifenResistant_genes.xls",header=T,sep="\t")
genes<-c(unique(as.character(tam$SYMBOL)))
pathways<-list()
pathways[["TAM"]]<-genes
pathways[["TAM_HK"]]<-c("ANKRD32","ABHD10","INTS12","SIRT3","TATDN1","UBC","CAV2","ATP5E","HIST1H2BM","RAB27B","RPLP1","SLC12A9","REEP6","IFITM2","NDUFS6","TSSC4","TMSB15B","HIST1H3E","C5NK2A2","ATP6VOE2")
pathways[["TAM_XMen"]]<-c("ESR1","IGFBP5","ADAMTS9","ALDH3B2","RTN4RL1","MYH10","EDN1","HGFAC","DHRS2","PDGFA","MAOA","FLT4","INS-IGF2","BAMBI","E2F8","ALDH5A1","ALDH6A1","ALDH4A1","CDKN2A","LHX9","CCNA1","COL8A1","ALDH3A1","LAMA1","CDKN1C","KLF12","CDK14","PNMT","VEGFC","EDNRA","COL9A3","ALDH1L2","CDK6","TIMP2","ROBO1","KLF7","KLF13","COL4A5","SKP2","SLIT2","CDKN2C","NPY1R","C2CD4D","SOX8","GABRP","GFRA3","PTPRN2","ARNT2","ATP6V1C2","GRB14","NELL2")
ER<-read.table("ER_associated.txt",header=T)[,1]
pathways[["ER"]]<-ER
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]
deGenesAll<-unique(as.character(expression$hgnc_symbol[(abs(expression$logFC)>=log2(1.5)) & (expression$adj.P.Val<=0.05)]))
deGenesPos<-unique(as.character(expression$hgnc_symbol[expression$logFC>=log2(1.5) & expression$adj.P.Val<=0.05]))
deGenesNeg<-unique(as.character(expression$hgnc_symbol[expression$logFC<=-log2(1.5) & expression$adj.P.Val<=0.05]))
pathways[["ELF5deGenesAll"]]<-deGenesAll
pathways[["ELF5deGenesPos"]]<-deGenesPos
pathways[["ELF5deGenesNeg"]]<-deGenesNeg
pathways[["CHIP_T47D_MCF7"]]<-deGenesPos

fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=10, maxSize=500,
                  nperm=100000)




fittedValues<-normalized_counts
fittedValues_genes<-fittedValues[row.names(fittedValues) %in% ER,]

library(pheatmap)
fittedValues_genesL=log1p(fittedValues_genes)
pdf(paste0(imageDir,"/JNCI_ER.heatmap.pdf"),width=10,height=12)
pheatmap(fittedValues_genesL,cluster_cols=F)
dev.off()

tam<-unique(tam$SYMBOL)

fittedValues_genes<-fittedValues[row.names(fittedValues) %in% tam,]

library(pheatmap)
fittedValues_genesL=log1p(fittedValues_genes)
pdf(paste0(imageDir,"/JNCI_TAM.heatmap.pdf"),width=10,height=12)
pheatmap(fittedValues_genesL,cluster_cols=F)
dev.off()

library(pheatmap)
fittedValues_genesL=log1p(fittedValues_genes)
fittedValues_genesLSorted<-fittedValues_genesL[,c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
pdf(paste0(imageDir,"/JNCI_TAM_by_patient.heatmap.pdf"),width=10,height=12)
pheatmap(fittedValues_genesLSorted,cluster_cols=F)
dev.off()



DEgenes4<-res[res$padj<0.1,]
symbols<-unique(DEgenes4$id)
DEgenes4[!is.na(DEgenes4$id),c(1,3,4,6,7)]
symbols<-symbols[!is.na(symbols)]




