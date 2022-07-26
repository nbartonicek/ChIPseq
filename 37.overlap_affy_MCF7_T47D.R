
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
system(paste("mkdir -p",robjectsDir))
scriptsPath=paste(projectDir,"/scripts")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
#names(chrs)=gsub("chr","",names(chrs))
#names(chrs)[25]="MT"
gr<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))

#########################################
########## 1. Load in the files
#########################################


jnci<-paste0(robjectsDir,"JNCI.Rdata")
if(!file.exists(jnci)){
  df<-read.table(paste0(resultsDir,"/JNCI/brainMetPairs.salmon.cts.txt"))

  annotation<-read.table("jnci_annotation.txt",sep="\t",header=T)

  #eliminate extra samples and reorder
  df<-df[,!c(colnames(df) %in% c("X7P_RCS","X7M_RCS"))]
  df<-df[,c(grep("P_",colnames(df)),grep("M_",colnames(df)))]


  #trying out mds
  group <- factor(paste0(c(annotation$ER,annotation$ER.1)))
  y <- DGEList(df,group=group)
  y <- calcNormFactors(y)

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
  save(dfShort,file=jnci)
} else {load(jnci)}




############# load DE genes #############
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]
deGenesAll<-unique(as.character(expression$hgnc_symbol[(abs(expression$logFC)>=log2(1.5)) & (expression$adj.P.Val<=0.05)]))
deGenesPos<-unique(as.character(expression$hgnc_symbol[expression$logFC>=log2(1.5) & expression$adj.P.Val<=0.05]))
deGenesNeg<-unique(as.character(expression$hgnc_symbol[expression$logFC<=-log2(1.5) & expression$adj.P.Val<=0.05]))

deGenesL<-list()
deGenesL[["DE"]]<-deGenesAll
deGenesL[["pos"]]<-deGenesPos
deGenesL[["neg"]]<-deGenesNeg
deGenesNeg<-deGenesNeg[!is.na(deGenesNeg)]

############# RNAseq MCF7 absolute values #############

mcf7<-paste0(robjectsDir,"MCF7_RNAseq.Rdata")
if(!file.exists(mcf7)){
  inDir<-paste0(resultsDir,"/ELF5_RNAseq.rsem/")
  inFiles<-list.files(inDir,pattern="genes",full.names=T)

  results<-list()
  for(inFile in inFiles){
    data<-read.table(inFile,header=T)
    sampleName<-basename(inFile)
    sampleName<-gsub(".*_V5_","",sampleName)
    sampleName<-gsub("_.*","",sampleName)
    cat(sampleName)
    cat("\n")
    results[[sampleName]]<-data$expected_count
  }

  df<-do.call("cbind",results)
  df<-as.data.frame(df)
  row.names(df)<-data$gene_id

  egENSEMBL <- toTable(org.Hs.egENSEMBL)
  ensemblIDs=gsub("\\..*","",row.names(df))
  m <- match(ensemblIDs, egENSEMBL$ensembl_id)
  df$EntrezGene<-egENSEMBL$gene_id[m]

  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  m <- match(df$EntrezGene, egSYMBOL$gene_id)
  df$symbol<-egSYMBOL$symbol[m]

  #eliminate duplicated symbols
  o <- order(rowSums(df[,1:6]), decreasing=TRUE)
  df=df[o,]
  d<-duplicated(df$symbol)
  df<-df[!d,]

  #add the gene length
  #data$gene_id=gsub("\\..*","",data$gene_id)
  #df_merged<-merge(df,data[,c(1,4)],by.x=0,by.y="gene_id")
  #df<-df_merged
  #row.names(df)<-df[,1]
  #df<-df[,-1]
  #eliminate lowly expressed
  include<-apply(df[,1:6],1,function(x){sum(x>=10)>=2})
  df<-df[include,]
  df=df[!is.na(df$symbol),]
  row.names(df)<-df$symbol
  save(df,file=mcf7)
} else {load(mcf7)}


library(DESeq2)

treatment <- rep(c("wt","dox"), times=c(3,3))

design <- model.matrix(~treatment)

meta<-data.frame(row.names=colnames(df)[1:6],treatment=treatment)
dfInt<-sapply(df[,1:6],function(x){round(x,digits=0)})

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData=dfInt, colData=meta, design=~treatment)
dds <- estimateSizeFactors(dds)
normalized_counts_mcf7_RNAseq <- counts(dds, normalized=TRUE)
row.names(normalized_counts_mcf7_RNAseq)<-row.names(df)


############# plot heatmap ################

treatment <- rep(c("P", "M"), times=c(7,7))
nSamples<-14
patient <- gsub("\\w_.*","",colnames(dfShort))[1:nSamples]

design <- model.matrix(~patient+treatment)
samples<-dfShort[1:nSamples]

meta<-data.frame(row.names=colnames(samples),treatment=rev(treatment),patient=patient)

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData=samples, colData=meta, design=~patient+treatment)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)



normalized_counts_mcf7_RNAseq_short<-as.data.frame(normalized_counts_mcf7_RNAseq[row.names(normalized_counts_mcf7_RNAseq) %in% deGenesNeg,])
normalized_counts_short<-as.data.frame(normalized_counts[row.names(normalized_counts) %in% deGenesNeg,])
normalized_counts_short<-normalized_counts_short[,c(1,8,7,14)]

normalized_counts_mcf7_RNAseq_short$id<-row.names(normalized_counts_mcf7_RNAseq_short)
normalized_counts_short$id<-row.names(normalized_counts_short)

merged<-merge(normalized_counts_mcf7_RNAseq_short,normalized_counts_short)

row.names(merged)<-merged[,1]
merged<-merged[,-1]




pdf("../project_results/figures/JNCI_DESeq2_normalised_DEneg.pdf",width=8,height=32)
pheatmap(log1p(merged),cluster_cols=F)
dev.off()

affyMCF7<-read.table("../project_results/GSE30407/MCF7_DE_all.tsv",header=T)
affyT47D<-read.table("../project_results/GSE30407/T47D_DE_all.tsv",header=T)

expressionL<-list()
expression
for(type in c("logFC","t")){

  affyMCF7_short<-affyMCF7[,c("Gene.symbol",type)]
  affyT47D_short<-affyT47D[,c("Gene.symbol",type)]
  expression_short<-expression[,c("hgnc_symbol",type)]
  merged<-merge(expression_short,affyMCF7_short,by.x="hgnc_symbol",by.y="Gene.symbol")
  merged<-merge(merged,affyT47D_short,by.x="hgnc_symbol",by.y="Gene.symbol")
  colnames(merged)<-c("id","MCF7_RNAseq","MCF7_microarray","T47D_microarray")
  dataM<-melt(merged)
  colnames(dataM)<-c("gene","sample",type)
  pdf(paste0(imageDir,"/MCF7_RNAseq_vs_T47D_microarray_t.pdf"),width=12,height=8)
  p<-ggplot(merged,aes(MCF7_RNAseq,T47D_microarray))
  p<-p+geom_point()
  p<-p+geom_density_2d()
  p
  dev.off()


}






#try GSEA

#DE genes for T47D
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
expression=expression[expression$hgnc_symbol != "",]




rank<-expression$t
names(rank)<-expression$hgnc_symbol
rank=sort(rank)
d<-duplicated(names(rank))
rank<-rank[!d]

arrays<-list()

arrays[["MCF7_array"]]<-affyMCF7$Gene.symbol[affyMCF7$adj.P.Val<0.05]
arrays[["T47D_array"]]<-affyT47D$Gene.symbol[affyT47D$adj.P.Val<0.05]
fgseaRes <- fgsea(pathways = arrays, 
                  stats = rank,
                  minSize=10, maxSize=1100,
                  nperm=100000)




rank<-affyMCF7$t
names(rank)<-affyMCF7$Gene.symbol
rank=sort(rank)
d<-duplicated(names(rank))
rank<-rank[!d]

arrays<-list()

arrays[["T47D_array"]]<-affyT47D$Gene.symbol[affyT47D$adj.P.Val<0.05]
fgseaRes <- fgsea(pathways = arrays, 
                  stats = rank,
                  minSize=10, maxSize=1100,
                  nperm=100000)



