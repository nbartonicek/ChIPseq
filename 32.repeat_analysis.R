

library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(org.Hs.eg.db)
library(circlize)

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
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")

system(paste0("mkdir -p ",cleanRobjectsDir))

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
load(paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))

for(file in list.files(peakRobjectsDir,pattern="Rdata",full.names=T)){
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


#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tss <- promoters(genes(txdb,columns=c("gene_id")), upstream=0, downstream=1)
#promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
geneid <- mapIds(org.Hs.eg.db, names(tss), "SYMBOL","ENTREZID")
tss=tss[!duplicated(geneid)]
names(tss)<-geneid[!duplicated(geneid)]
#eliminate duplicated


#DE genes
expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)

rank<-expression$t
names(rank)<-expression$hgnc_symbol
rank=sort(rank)
DEgenes<-rev(rank[c(1:50,c(length(rank)-49):length(rank))])

tssDE<-tss[names(tss) %in% names(DEgenes)]
tssDE<-tssDE[names(DEgenes)]

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 5000, mean_mode = "w0", w = 50)

genesDf<-data.frame(genes=row.names(mat1),direction="down",stringsAsFactors=F)
genesDf$direction[c(1:dim(genesDf)[1])[genesDf$genes %in% names(rank)[rank>0]]]="up"

pdf(paste0(imageDir,"dist_",sampleName,"_DE.pdf"),width=6,height=6)
col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1, split=genesDf$direction,name =sampleName)
dev.off()


######################

load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["ELF5Dox"]]
sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 5000, mean_mode = "w0", w = 50)

#genesDf<-data.frame(genes=row.names(mat1),direction="down",stringsAsFactors=F)
#genesDf$direction[c(1:dim(genesDf)[1])[genesDf$genes %in% names(rank)[rank>0]]]="up"
sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]
mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 5000, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 5000, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 5000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_FOXA1_ER_ELF5_all.pdf"),width=6,height=6)
col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1, col=col_fun, name ="FOXA1Dox")+
  EnrichedHeatmap(mat2, col=col_fun, name ="FOXA1NoDox")+
  EnrichedHeatmap(mat3, col=col_fun, name ="ERDox")+
  EnrichedHeatmap(mat4, col=col_fun, name ="ERNoDox")
dev.off()

########### plot diff stuff

load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["ELF5Dox"]]
sampleName="FOXA1_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 500, mean_mode = "w0", w = 50)

#genesDf<-data.frame(genes=row.names(mat1),direction="down",stringsAsFactors=F)
#genesDf$direction[c(1:dim(genesDf)[1])[genesDf$genes %in% names(rank)[rank>0]]]="up"
sampleName="FOXA1_diff_Dox_depleted"
gr=cleanGRsPeaks[[sampleName]]
mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 500, mean_mode = "w0", w = 50)

sampleName="ER_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 500, mean_mode = "w0", w = 50)

sampleName="ER_diff_Dox_depleted"
gr=cleanGRsPeaks[[sampleName]]
mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 500, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_FOXA1_ER_ELF5_diff.pdf"),width=6,height=6)
#col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1,name ="FOXA1Dox")+
  EnrichedHeatmap(mat2,name ="FOXA1NoDox")+
  EnrichedHeatmap(mat3,name ="ERDox")+
  EnrichedHeatmap(mat4,name ="ERNoDox")
dev.off()

########### plot diff stuff, but overlapped with ELF5+FOXA1 only


load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["ELF5Dox"]]
sampleName="FOXA1_diff_Dox_enriched"

gr=cleanGRsPeaks[[sampleName]]
mat<-findOverlaps(tssDE,gr)
tssDE<-tssDE[unique(queryHits(mat))]

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]

mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]
mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_FOXA1_ER_ELF5overlap_diff2.pdf"),width=12,height=8)
#col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1, col = c("white", "darkred"),name ="FOXA1_diff_Dox_enriched",column_title = "FOXA1NoDox")+
  EnrichedHeatmap(mat2, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
  EnrichedHeatmap(mat3, col = c("white", "darkred"),name ="ER_diff_Dox_enriched",column_title = "ERDox")+
  EnrichedHeatmap(mat4, col = c("white", "darkblue"),name ="ERNoDox",column_title = "ERNoDox")
dev.off()

########### plot diff stuff, but overlapped with ELF5+FOXA1 only and sorted by FOXA1NoDox


load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["ELF5Dox"]]
sampleName="FOXA1_diff_Dox_enriched"

gr=cleanGRsPeaks[[sampleName]]
mat<-findOverlaps(tssDE,gr)
tssDE<-tssDE[unique(queryHits(mat))]


sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]

mat1 = normalizeToMatrix(gr, tssDE, value_column = "pval", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]

mat2 = normalizeToMatrix(gr, tssDE, value_column = "pval", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="FOXA1_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
#gr$score=10*(gr$pval)
mat3 = normalizeToMatrix(gr, tssDE, value_column = "pval", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
#gr$score=10*(gr$pval)

mat4 = normalizeToMatrix(gr, tssDE, value_column = "pval", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
mat5 = normalizeToMatrix(gr, tssDE, value_column = "pval", 
    extend = 1500, mean_mode = "w0", w = 50)


sampleName="ER_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
mat6 = normalizeToMatrix(gr, tssDE, value_column = "pval", 
    extend = 1500, mean_mode = "w0", w = 50)


pdf(paste0(imageDir,"dist_FOXA1_ER_ELF5overlap_diff_FoxA1NoDoxSorted.pdf"),width=14,height=8)
col_fun1 = colorRamp2(quantile(mat1, c(0, 0.945)), c("white", "darkblue"))
col_fun2 = colorRamp2(quantile(mat2, c(0, 0.82)), c("white", "darkred"))
col_fun4 = colorRamp2(quantile(mat4, c(0, 0.99522)), c("white", "darkblue"))
col_fun5 = colorRamp2(quantile(mat5, c(0, 0.985)), c("white", "darkred"))

EnrichedHeatmap(mat1, col = col_fun1,name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
EnrichedHeatmap(mat2, col = col_fun2,name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat3, col = c("white", "darkgreen"),name ="FOXA1_diff_Dox_enriched",column_title = "FOXA1_diff_Dox_enriched")+
EnrichedHeatmap(mat4, col = col_fun4, name ="ERNoDox",column_title = "ERNoDox") +
EnrichedHeatmap(mat5, col = col_fun5, name ="ERDox",column_title = "ERDox") +
EnrichedHeatmap(mat6, col = c("white", "darkgreen"),name ="ER_diff_Dox_enriched",column_title = "ER_diff_Dox_enriched")
dev.off()

########### version without diff

########### plot diff stuff, but overlapped with ELF5+FOXA1 only and sorted by FOXA1NoDox


load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["ELF5Dox"]]
sampleName="FOXA1_diff_Dox_enriched"

gr=cleanGRsPeaks[[sampleName]]
mat<-findOverlaps(tssDE,gr)
tssDE<-tssDE[unique(queryHits(mat))]


sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)


pdf(paste0(imageDir,"dist_FOXA1_ER_ELF5overlap_Nodiff_FoxA1NoDoxSorted.pdf"),width=10,height=8)
col_fun1 = colorRamp2(quantile(mat1, c(0, 0.945)), c("white", "darkblue"))
col_fun2 = colorRamp2(quantile(mat2, c(0, 0.82)), c("white", "darkred"))
col_fun4 = colorRamp2(quantile(mat4, c(0, 0.99522)), c("white", "darkblue"))
col_fun5 = colorRamp2(quantile(mat5, c(0, 0.985)), c("white", "darkred"))

EnrichedHeatmap(mat1, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
EnrichedHeatmap(mat2, col = c("white", "darkred"),name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat4, col = c("white", "darkblue"), name ="ERNoDox",column_title = "ERNoDox") +
EnrichedHeatmap(mat5, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox")
dev.off()



########### plot diff stuff, but overlapped with ELF5+FOXA1 depleted only


load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["ELF5Dox"]]
sampleName="FOXA1_diff_Dox_depleted"

gr=cleanGRsPeaks[[sampleName]]
mat<-findOverlaps(tssDE,gr)
tssDE<-tssDE[unique(queryHits(mat))]

sampleName="FOXA1_diff_Dox_depleted"
gr=cleanGRsPeaks[[sampleName]]

mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]
mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_FOXA1depleted_ER_ELF5overlap_diff2.pdf"),width=12,height=6)
#col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1, col = c("darkred", "white"),name ="FOXA1_diff_Dox_enriched",column_title = "FOXA1_diff_Dox_depleted")+
  EnrichedHeatmap(mat2, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
  EnrichedHeatmap(mat3, col = c("white", "darkred"),name ="ERDox",column_title = "ERDox")+
  EnrichedHeatmap(mat4, col = c("white", "darkblue"),name ="ERNoDox",column_title = "ERNoDox")
dev.off()

########## ER centred

load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

tssDE<-cleanGRs[["ELF5Dox"]]
sampleName="ER_diff_Dox_enriched"

gr=cleanGRsPeaks[[sampleName]]
mat<-findOverlaps(tssDE,gr)
tssDE<-tssDE[unique(queryHits(mat))]

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]

mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]
mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1500, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ER_FOXA1_ELF5overlap_diff.pdf"),width=12,height=8)
#col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat1, col = c("white", "darkred"),name ="FOXA1_diff_Dox_enriched",column_title = "FOXA1NoDox")+
  EnrichedHeatmap(mat2, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
  EnrichedHeatmap(mat3, col = c("white", "darkred"),name ="ER_diff_Dox_enriched",column_title = "ERDox")+
  EnrichedHeatmap(mat4, col = c("white", "darkblue"),name ="ERNoDox",column_title = "ERNoDox")
dev.off()














