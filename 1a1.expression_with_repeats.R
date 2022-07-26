
library(GenomicRanges)
library(ShortRead)
#library(R.utils)
library("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
#library(RUVSeq)
library(org.Mm.eg.db)
library(DESeq)

timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84



homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../../../"
projectnames=c( "lx9", "NOD", "B6_CVB4", "CVB4_D4" )



######## directory structure #######
projectDir=paste0(homedir,"/projects/claudia")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/API/")
projectname="lx9"
robjectsDir = paste(resultsDir,"/",projectname,".Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
sizeDir=paste(resultsDir,"/",projectname,".libSize/",sep="")

outPath=paste0(resultsDir,"/",projectname,".repeatOverlap/")
system(paste("mkdir",outPath))
system(paste("mkdir",sizeDir))

chrs=seqlengths(Mmusculus)[!grepl("_",names(seqlengths(Mmusculus)))]
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
projectnames=c( "lx9", "B6_CVB4", "CVB4_D4" )
geneResults<-list()

for(projectname in projectnames){
  inPath=paste0(homedir,"/projects/claudia/project_results/",projectname,".rsemSequin/")

  geneFiles<-list.files(inPath,pattern="genes",full.names=T)


  for(file in geneFiles){
    sampleName<-basename(file)
    sampleName<-gsub("_M001.*","_M001",sampleName)
    cat(sampleName)
    cat("\n")
    data<-read.table(file,header=T)  
    geneResults[[sampleName]]<-as.integer(data$expected_count)
  }
  temp<-geneResults
}
df<-as.data.frame(geneResults)
row.names(df)<-data$gene_id

#load the repetitive elements
repeatResults<-list()

inPath="/share/ScratchGeneral/nenbar/projects/claudia/project_results/repeatOverlap"
repeatFiles<-list.files(inPath,pattern="type",full.names=T)
for(file in repeatFiles){
  sampleName<-basename(file)
  sampleName<-gsub("_type.Rdata","",sampleName)
  cat(sampleName)
  cat("\n")
  load(file)
  repeatResults[[sampleName]]<-countsType
}

dfRepeat<-do.call("cbind",repeatResults)
dfRepeat<-dfRepeat[,colnames(dfRepeat) %in% colnames(df)]

df<-rbind(df,dfRepeat)


sampleNames<-read.table("samples.txt",header=F)

for(i in 1:length(geneResults)){
  id=colnames(df)[i]
  id=gsub(".*_FD","FD",id)
  id=gsub("_.*","",id)
  newName<-as.character(sampleNames$V2[sampleNames$V1 %in% id])
  colnames(df)[i]=newName
}

#########################################
########## 1. Load in the files
#########################################



cols<-brewer.pal(6,"Paired")
pdf("../../project_results/figures/samples_pca_repeats.pdf",width=12,height=8)
pca<-princomp(df)
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, names(geneResults),pos = 1)
dev.off()

pdf("../../project_results/figures/samples_pca_components_repeats.pdf",width=12,height=8)
plot(pca)
dev.off()

#trying out mds
group=gl(4, c(6,3,3,3), labels = c("LX9", "B6", "B6_CVB4", "B6_CVB4_D4"))
group <- factor( rep(c("LX9", "B6", "B6_CVB4", "B6_CVB4_D4"), times=c(6,3,3,3)) )
y <- DGEList(df,group=group)
y <- calcNormFactors(y)

# without labels, indexes of samples are plotted.
col <- as.numeric(group)
pdf("../../project_results/figures/samples_MDS_repeats.pdf",width=12,height=8)
mds <- plotMDS(y, top=200, col=col,labels=group)
dev.off()

# or labels can be provided, here group indicators:
pdf("../../project_results/figures/samples_MDS_labels_repeats.pdf",width=12,height=8)

plotMDS(mds, col=col, labels=group)
dev.off()

#plot(pca)
#dev.off()

#Annotate with symbols, aggregate

#Annotate with entrez, aggregate
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df)=gsub("\\..*","",row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
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
include<-apply(df[,1:15],1,function(x){sum(x>=50)>=2})
df<-df[include,]
df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol

#change the id to entrez
#o <- order(rowSums(df[,1:6]), decreasing=TRUE)
#df=df[o,]
#d<-duplicated(df$EntrezGene)
#df<-df[!d,]
#row.names(df)<-df$EntrezGene


#provide the EC numbers
egENZYME <- toTable(org.Mm.egENZYME)
m <- match(row.names(df), egENZYME$gene_id)

df[,"EC number"]<-egENZYME$ec_number[m]




cols<-brewer.pal(6,"Paired")
pdf("../../project_results/figures/pca_clean_all_repeats.pdf",width=12,height=8)
pca<-princomp(df[1:6])
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, names(geneResults),pos = 1)
dev.off()

pdf("../../project_results/figures/pca_components_clean_all_repeats.pdf",width=12,height=8)
plot(pca)
dev.off()

#Make a table and submit to degust


treatment <- rep(c("inflammed", "control"), 3)
type <- c(rep(c("B6", "Nod", "SJL"), each=2))

design <- model.matrix(~type + treatment)

#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=df[1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt <- glmLRT(fit)


logFCs<-as.data.frame(topTags(lrt,n=1000))
plus<-logFCs[logFCs$logFC>0,]
minus<-logFCs[logFCs$logFC<0,]
significantIDs<-c(row.names(plus[1:25,]),row.names(minus[1:25,]))
fittedValues<-fitted.values(lrt)
fittedValues50<-fittedValues[row.names(fittedValues) %in% significantIDs,]
colnames(fittedValues50)<-gsub("_R1","",colnames(fittedValues50))

pdf("../../project_results/figures/logFC_genes_25plus25minus_repeats.pdf",width=12,height=8)
plot(pca)
#dev.off()

#expr_norm <- rpkm(expr, log=TRUE, gene.length=df$effective_length)
#save(expr_norm,file="normalized_expression.Rdata")
#pdf("../../project_results/figures/pca_fpkm.pdf",width=12,height=8)
#pca<-princomp(expr_norm)
#plot(pca$loading,pch=19, cex=2,col=cols)
#text(pca$loading, names(geneResults),pos = 1)
#dev.off()

#pdf("../../project_results/figures/pca_components_fpkm.pdf",width=12,height=8)
#plot(pca)
#dev.off()

go <- goana(lrt, species="Mm")
topGO(go)
write.table(topGO(go,number=Inf),"GO_categories_repeats.txt",row.names=T,quote=F,sep="\t")

keg <- kegga(lrt, species="Mm")
topKEGG(keg)
write.table(topKEGG(keg,number=Inf),"KEGG_categories_repeats.txt",row.names=T,quote=F,sep="\t")
break()

row.names(df)<-df$symbol
expr <- DGEList(counts=df[1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt <- glmLRT(fit)
topTags(lrt)

d=topTags(lrt,n=10000,p.value=0.01)
write.table(d,"all_inf.txt",quote=F,sep="\t",row.names=T)




################## DESEQ of the individual samples

grouping <- colnames(df)[1:6]
#grouping <- factor(type)

countTable=df[,1:6]
condition=grouping
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds, method='blind',sharingMode="fit-only" )

results<-list()
for(i in c(1,3,5)){
    res = nbinomTest( cds, colnames(df)[i+1],colnames(df)[i])
    sampleName=gsub("_.*","",colnames(df)[i])
    res[res$baseMeanA==0,"log2FoldChange"]<-log2(res[res$baseMeanA==0,"baseMeanB"])
    res[(res$baseMeanA==0&res$baseMeanB==0),"log2FoldChange"]<-0
    res[res$baseMeanB==0,"log2FoldChange"]<-log2(res[res$baseMeanB==0,"baseMeanB"])
    res[(res$baseMeanA==0&res$baseMeanB==0),"log2FoldChange"]<-0
    results[[sampleName]]<-res

}

diffExprs<-do.call(cbind,results)
colnames(diffExprs)<-paste0(rep(names(results),each=8),"_",colnames(res))
colnames(diffExprs)[1]<-"id"
diffExprs<-diffExprs[,-grep("_id",colnames(diffExprs))]

write.table(diffExprs,"pairwise_diff_expression_repeats.txt",row.names=F,quote=F,sep="\t")



normalized.count <- counts(cds, normalized = TRUE)
normalized.count<-normalizeQuantiles(normalized.count)
save(normalized.count,file="../project_results//endo_1.tables/normalized.counts_repeats.Rdata")

############ heatmap




pdf("../project_results/endo_1.figures/normalized_pca2_repeats.pdf",width=12,height=8)
pca<-princomp(normalized.count)
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, names(geneResults),pos = 1)
dev.off()

annotation<-read.table("../project_results//endo_1.tables/annotation.txt",header=T)
#combine the names and put into degust
normalized.count<-as.data.frame(normalized.count)
normalized.count$id<-row.names(normalized.count)
merged<-merge(normalized.count,annotation,by.x="id",by.y="gene_id")
write.table(merged,"../project_results/endo_1.tables/normalized.count_repeats.txt",row.names=T,quote=F,sep="\t")











