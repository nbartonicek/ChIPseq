#data is downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66733

#library(GenomicRanges)
library(ggplot2)
#library(rtracklayer)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(data.table)
#library(corrplot)
#library(ChIPpeakAnno)
library(reshape2)
#library(pheatmap)
#library(UpSetR) 

library(migest)
library(circlize)
library(dplyr)
library(tidyr)

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
imageDir=paste0(resultsDir,"/figures/circos/")
system(paste0("mkdir -p ",imageDir))
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
grammarDir=paste0(resultsDir,"/grammar/")
permRobjectsDir=paste(resultsDir,"/Robjects/grammarPermutations/",sep="")

load(paste0(grammarDir,"all_consensus_clean.Rdata"))
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
#genes<-list() 
#genes[["DE"]]=sort(DE_ol,decreasing=T)[1:30]
#genes[["superenhancer"]]=sort(superenhancers_ol,decreasing=T)[1:30]
#genes[["enhancer"]]=sort(enhancers_ol,decreasing=T)[1:20]
#genes[["sensitive"]]=sort(sensitive_ol,decreasing=T)[1:20]
#genes[["resistant"]]=sort(resistant_ol,decreasing=T)[1:20]

genes<-list() 
genes[["DE"]]=sort(DE_ol,decreasing=T)
genes[["superenhancer"]]=sort(superenhancers_ol,decreasing=T)
genes[["enhancer"]]=sort(enhancers_ol,decreasing=T)
genes[["sensitive"]]=sort(sensitive_ol,decreasing=T)
genes[["resistant"]]=sort(resistant_ol,decreasing=T)

#now for each first make upset plots, and then do the heatmaps
#consensusReps contains TFS
#targets are in "targets"

results<-list()
width=2000
for(regionName in names(genes)){

	#union
	TFs<-gsub(paste0(regionName,"_"),"",names(genes[[regionName]]))
	TFs<-sort(TFs)
	allPeaks<-reduce(targets[[regionName]])
	#countOverlaps
	tf<-list()
	for(sampleName in TFs){
	  cat(".")
	  #find overlaps of TF (within 1kb)
	  mat<-findOverlaps(consensusReps[[sampleName]],allPeaks)
	  subsample<-consensusReps[[sampleName]][queryHits(mat)]
	  subsample<-resize(subsample,width=width,fix="center")
	  tf[[sampleName]]<-countOverlaps(consensusReps,subsample)
	}
	#matrix
	df<-do.call("rbind",tf)
	results[[regionName]]<-df
	#umap

	#pdf(paste0(imageDir,"grammar_",regionName,".pdf"),width=10,height=10)
	#upset(df, order.by = "freq",nsets=length(TFs))
	#dev.off()
}
save(results,file=paste0(robjectsDir,"circos_results_",width,".Rdata"))


subselection<-50

SE<-results[["superenhancer"]]
TE<-results[["enhancer"]]
DE<-results[["DE"]]

df1<-SE[,colnames(SE) %in% row.names(SE)]
df1<-df1[order(row.names(df1)),order(colnames(df1))]


#df1 <- df0 %>% select(1:3) %>% rename(order = V1, rgb = V2, region = V3) %>% mutate(region = gsub("_", " ", region))

m <- as.matrix(df1)

dimnames(m) <- list(orig = row.names(df1), dest = colnames(df1))
#m[m<=quantile(m,0.5)]<-0
cn<-colnames(m)
#cnShort<-cn[order(colSums(df1),decreasing=T)][1:subselection]
cnShort<-cn[order(colSums(df1),decreasing=T)]

#m<-m[rowSums(m)>2000,rowSums(m)>2000]
m<-m[row.names(m) %in% cnShort,colnames(m) %in% cnShort]
#diag(m)<-0
#m<-apply(m,1,function(x){x/max(x)})
#m[m<=quantile(m,0.65)]<-0

m1=m
pdf(paste0(imageDir,"filteredSE.pdf"),width=18,height=16)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(x = m, directional = 1, order = row.names(m), 
             annotationTrack = "grid", 
             transparency = 0,  annotationTrackHeight = c(0.1, 0.1),
             diffHeight  = -0.04)
dev.off()


pdf(paste0(imageDir,"simpleSE.pdf"),width=18,height=16)
chordDiagram(m)
dev.off()

mSE=m

df1<-TE[,colnames(TE) %in% row.names(TE)]
df1<-df1[order(row.names(df1)),order(colnames(df1))]


#df1 <- df0 %>% select(1:3) %>% rename(order = V1, rgb = V2, region = V3) %>% mutate(region = gsub("_", " ", region))

m <- as.matrix(df1)

dimnames(m) <- list(orig = row.names(df1), dest = colnames(df1))
#m[m<=quantile(m,0.5)]<-0
cn<-colnames(m)
cnShort<-cn[order(colSums(df1),decreasing=T)][1:subselection]

#m<-m[rowSums(m)>2000,rowSums(m)>2000]
#m<-m[row.names(m) %in% cnShort,colnames(m) %in% cnShort]
#diag(m)<-0

m<-apply(m,1,function(x){x/max(x)})
#m[m<=quantile(m,0.65)]<-0

pdf(paste0(imageDir,"filteredTE.pdf"),width=18,height=16)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(x = m, directional = 0, order = row.names(m), 
             annotationTrack = "grid", 
             transparency = 0,  annotationTrackHeight = c(0.1, 0.1),
             diffHeight  = -0.04)
dev.off()

mTE=m
pdf(paste0(imageDir,"simpleTE.pdf"),width=18,height=16)
chordDiagram(m)
dev.off()




circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  sector.index = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), 2.5, sector.index, facing = "bending")
  circos.axis("top", major.at = seq(0, max(xlim)), minor.ticks=1, labels.away.percentage = 0.2, labels.niceFacing = FALSE )
  }, bg.border = NA)







df1<-SE[,colnames(SE) %in% row.names(SE)]
df1<-df1[order(row.names(df1)),order(colnames(df1))]


#df1 <- df0 %>% select(1:3) %>% rename(order = V1, rgb = V2, region = V3) %>% mutate(region = gsub("_", " ", region))

m <- as.matrix(df1)

dimnames(m) <- list(orig = row.names(df1), dest = colnames(df1))
#m[m<=quantile(m,0.5)]<-0
cn<-colnames(m)
cnShort<-cn[order(colSums(df1),decreasing=T)][1:subselection]

#m<-m[rowSums(m)>2000,rowSums(m)>2000]
#m<-m[row.names(m) %in% cnShort,colnames(m) %in% cnShort]
diag(m)<-0
m<-apply(m,1,function(x){x/max(x)})
#m[m<=quantile(m,0.65)]<-0

m1=m
pdf(paste0(imageDir,"filteredDE.pdf"),width=18,height=16)
circos.clear()
par(mar = rep(0, 4), cex=0.9)
circos.par(start.degree = 90, gap.degree = 4)

chordDiagram(x = m, directional = 0, order = row.names(m), 
             annotationTrack = "grid", 
             transparency = 0,  annotationTrackHeight = c(0.1, 0.1),
             diffHeight  = -0.04)
dev.off()


pdf(paste0(imageDir,"pheatmap_pvals_.pdf"),width=8,height=16)
#pheatmap(dfLength,annotation_row=data.frame(class=ELF5$class,tamoxifen=ELF5$resistance))
pheatmap(wide,show_rownames = T,color = colorRampPalette(c("white", "blue", "darkblue"))(11))
dev.off()










elf5=cleanGRs[[1]]
elf5<-resize(elf5,width=1000,fix="center")

mat<-findOverlaps(elf5,targets[["enhancer"]])
te5<-elf5[queryHits(mat)]

mat<-findOverlaps(elf5,targets[["superenhancer"]])
se5<-elf5[queryHits(mat)]

mat<-findOverlaps(elf5,targets[["DE"]])
de5<-elf5[queryHits(mat)]

#now check how many different do they overlap

consensusReps=consensusReps[grepl("vehicle",names(consensusReps))]
names(consensusReps)<-toupper(names(consensusReps))

sampleNames<-gsub(".VEHICLE.*","",names(consensusReps))

cR<-split(consensusReps,sampleNames)
cR<-GRangesList(unlist(cR))


enhancerOverlaps<-countOverlaps(te5,cR)
superEnhancerOverlaps<-countOverlaps(se5,cR)
deOverlaps<-countOverlaps(de5,cR)





enh<-table(enhancerOverlaps)
sEnh<-table(superEnhancerOverlaps)
de<-table(deOverlaps)

dfE<-as.data.frame(enh,as.factor=F)
dfSE<-as.data.frame(sEnh,as.factor=F)
dfDE<-as.data.frame(de,as.factor=F)

merged<-merge(dfE,dfSE,by.x="enhancerOverlaps",by.y="superEnhancerOverlaps",all=T)
merged<-merge(merged,dfDE,by.x="enhancerOverlaps",by.y="deOverlaps",all=T)

merged[is.na(merged)]=0
colnames(merged)<-c("ELF5_partners_per_peak","enhancers","superenhancers","DE")
merged$enhancers=merged$enhancers/sum(merged$enhancers)
merged$superenhancers=merged$superenhancers/sum(merged$superenhancers)
merged$DE=merged$DE/sum(merged$DE)

merged$ELF5_partners_per_peak=as.character(merged$ELF5_partners_per_peak)

merged$ELF5_partners_per_peak=factor(merged$ELF5_partners_per_peak,levels=as.character(sort(as.numeric(merged$ELF5_partners_per_peak))))
#merged=merged[order(merged$ELF5_partners_per_peak),]
#merged$ELF5_partners_per_peak=as.numeric(merged$ELF5_partners_per_peak)
dataM<-melt(merged)
#dataM$ELF5_partners_per_peak=factor(dataM$ELF5_partners_per_peak,levels=sort(dataM$ELF5_partners_per_peak))
dataM$ELF5_partners_per_peak=as.numeric(dataM$ELF5_partners_per_peak)

colnames(dataM)<-c("ELF5_partners_per_peak","type","count")


pdf(paste0(imageDir,"connectivity_comparison.pdf"),width=12,height=8)
p<-ggplot(dataM,aes(ELF5_partners_per_peak,count,group=type))
p<-p+geom_bar(stat="identity",aes(fill = type), position='dodge')
#p<-p+geom_density()
p
dev.off()



pdf(paste0(imageDir,"connectivity_density.pdf"),width=12,height=8)
p<-ggplot(dataM,aes(x=ELF5_partners_per_peak,y=count,colour=type))
#p<-p+geom_bar(stat="identity",aes(fill = type), position='dodge')
p<-p+geom_smooth(aes(colour = type),span = 0.5)
p
dev.off()











countOverlaps(consensusReps,de5S)

#is cohesin really so unique?

#find Cohesin in one and other, then overlap with all, find overlaps

coh<-consensusReps[["COHESIN.estradiol.p7"]]

mat<-findOverlaps(coh,targets[["enhancer"]])
cohE=coh[unique(queryHits(mat))]


mat<-findOverlaps(coh,targets[["superenhancer"]])
cohSE=coh[unique(queryHits(mat))]


cohEOverlaps<-countOverlaps(cohE[1:188],consensusReps)
cohSEOverlaps<-countOverlaps(cohSE,consensusReps)























