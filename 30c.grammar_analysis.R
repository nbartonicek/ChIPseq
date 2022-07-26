#data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66733

library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(corrplot)
library(ChIPpeakAnno)
library(reshape2)
library(pheatmap)
library(UpSetR) 

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
chrGR<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))
homedir="/share/ScratchGeneral/nenbar"
#homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs/")
inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate/")

inBams=paste0(homedir,"/projects/Chris/project_results/ELF5.picard/")
outPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate_consensus/")
system(paste0("mkdir -p ",outPath))

######## directory structure #######
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/grammar/")
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
grammarDir=paste0(resultsDir,"/grammar/")
permRobjectsDir=paste(resultsDir,"/Robjects/grammarPermutations/",sep="")

#load(paste0(grammarDir,"all_consensus_clean.Rdata"))
load(paste0(grammarDir,"targets.Rdata"))

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

#now for each in enhancer, superenhancer, DE, HOT, resistant, sensitive
#get p-values with ELF5

inDir=permRobjectsDir
files<-list.files(inDir,full.names=T)
files=files[!grepl("HOT",files)]
pvals<-list()
overlaps<-list()
res<-list()
for(file in files){
	sampleName<-gsub(".Rdata","",basename(file))
	sampleName<-gsub("peakOverlap_","",sampleName)
	cat(".")
	load(file)
	pvals[[sampleName]]<-pt$cntOverlaps$pval 
	overlaps[[sampleName]]<-pt$cntOverlaps$observed 
	res[[sampleName]]<-data.frame(pval=pt$cntOverlaps$pval,overlap=pt$cntOverlaps$observed)
}

DE<-pvals[grepl("DE_",names(pvals))]
superenhancers<-pvals[grepl("superenhancer",names(pvals))]
enhancers<-pvals[grepl("^enhancer",names(pvals))]
sensitive<-pvals[grepl("sensitive",names(pvals))]
resistant<-pvals[grepl("resistant",names(pvals))]

DE_ol<-overlaps[grepl("DE_",names(overlaps))]
superenhancers_ol<-overlaps[grepl("superenhancer",names(overlaps))]
enhancers_ol<-overlaps[grepl("^enhancer",names(overlaps))]
sensitive_ol<-overlaps[grepl("sensitive",names(overlaps))]
resistant_ol<-overlaps[grepl("resistant",names(overlaps))]


#take all the significant
overlapsShort<-overlaps[unlist(pvals)<0.0001]
DE_ol<-unlist(overlapsShort[grepl("DE_",names(overlapsShort))])
superenhancers_ol<-unlist(overlapsShort[grepl("superenhancer",names(overlapsShort))])
enhancers_ol<-unlist(overlapsShort[grepl("^enhancer",names(overlapsShort))])
sensitive_ol<-unlist(overlapsShort[grepl("sensitive",names(overlapsShort))])
resistant_ol<-unlist(overlapsShort[grepl("resistant",names(overlapsShort))])

#take all the top 10 
genes<-list() 
genes[["DE"]]=sort(DE_ol,decreasing=T)[1:30]
genes[["superenhancer"]]=sort(superenhancers_ol,decreasing=T)[1:30]
genes[["enhancer"]]=sort(enhancers_ol,decreasing=T)[1:20]
genes[["sensitive"]]=sort(sensitive_ol,decreasing=T)[1:20]
genes[["resistant"]]=sort(resistant_ol,decreasing=T)[1:20]


#now for each first make upset plots, and then do the heatmaps
#consensusReps contains TFS
#targets are in "targets"

for(regionName in names(genes)){

	#union
	TFs<-gsub(paste0(regionName,"_"),"",names(genes[[regionName]]))
	allPeaks<-reduce(targets[[regionName]])
	#countOverlaps
	for(sampleName in TFs){
	  cat(sampleName)
	  cat("\n")
	  values(allPeaks)[[sampleName]]<-countOverlaps(allPeaks,consensusReps[sampleName])
	}
	#matrix
	df<-as.data.frame(values(allPeaks))

	#umap

	#pdf(paste0(imageDir,"grammar_",regionName,".pdf"),width=10,height=10)
	#upset(df, order.by = "freq",nsets=length(TFs))
	#dev.off()
}






#the problem is that everything is statistically significant.
#which is shit

df<-do.call("rbind",res)

df$type=gsub("_.*","",row.names(df))

df$logPval<-(log1p(1/df$pval))
df$TF<-gsub("\\.p.*","",row.names(df))
df$TF<-gsub("^.*_(.*?)","\\1",df$TF)
df=df[,c(3:5)]
wide<-reshape(df,idvar="TF",timevar="type",direction="wide")
wide=wide[,-1]
colnames(wide)=gsub("logPval.","",colnames(wide))


row.names(wide)=gsub("DE_","",row.names(wide))


pdf(paste0(imageDir,"pheatmap_pvals_.pdf"),width=8,height=16)
#pheatmap(dfLength,annotation_row=data.frame(class=ELF5$class,tamoxifen=ELF5$resistance))
pheatmap(wide,show_rownames = T,color = colorRampPalette(c("white", "blue", "darkblue"))(11))
dev.off()














































