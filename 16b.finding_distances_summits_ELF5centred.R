

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

system(paste0("mkdir -p ",cleanRobjectsDir))

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


#first find overlaps between ELF5, FOXA1enriched and ER 
mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],cleanGRsPeaks[["FOXA1_diff_Dox_enriched"]])
both<-cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat))]
mat<-findOverlaps(both,cleanGRsPeaks[["ER_diff_Dox_enriched"]])
threeDiff<-both[unique(queryHits(mat))]

mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],cleanGRsPeaks[["FOXA1Dox"]])
both<-cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat))]
mat<-findOverlaps(both,cleanGRsPeaks[["ERDox"]])
threeAll<-both[unique(queryHits(mat))]

#then find the distance of FOXA1Dox and ERDox summits to ELF5 summits
matDiff<-findOverlaps(cleanGRs[["ELF5Dox"]],threeDiff)
matAll<-findOverlaps(cleanGRs[["ELF5Dox"]],threeAll)

ELF5SummitsDiff<-cleanGRs[["ELF5Dox"]][unique(queryHits(matDiff))]
ELF5SummitsAll<-cleanGRs[["ELF5Dox"]][unique(queryHits(matAll))]

w=2
extend=100

ELF5Dox<-cleanGRs[["ELF5Dox"]]
FOXA1Dox<-cleanGRs[["FOXA1Dox"]]
ERDox<-cleanGRs[["ERDox"]]
ELF5Dox$score=log10(ELF5Dox$score)
FOXA1Dox$score=log10(FOXA1Dox$score)
ERDox$score=log10(ERDox$score)

mat0Diff = normalizeToMatrix(ELF5Dox, ELF5SummitsDiff, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1Diff = normalizeToMatrix(FOXA1Dox, ELF5SummitsDiff, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2Diff = normalizeToMatrix(ERDox, ELF5SummitsDiff, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat0All = normalizeToMatrix(ELF5Dox, ELF5SummitsAll, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1All = normalizeToMatrix(FOXA1Dox, ELF5SummitsAll, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2All = normalizeToMatrix(ERDox, ELF5SummitsAll, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("dist_ELF5_summits_FOXA1_ER_diff_log_",w,"_",extend,".pdf")),width=6,height=6)
#col_fun = colorRamp2(quantile(mat1, c(0, 0.99)), c("white", "red"))
EnrichedHeatmap(mat0Diff,name ="ELF5Dox")+
EnrichedHeatmap(mat1Diff,name ="FOXA1Dox")+
EnrichedHeatmap(mat2Diff,name ="ERDox")
dev.off()



pdf(paste0(imageDir,paste0("dist_ELF5_summits_FOXA1_ER_all_log_",w,"_",extend,".pdf")),width=6,height=6)
col_fun1 = colorRamp2(quantile(mat0All, c(0, 0.9999)), c("white", "red"))
col_fun2 = colorRamp2(quantile(mat1All, c(0, 0.9999)), c("white", "red"))
col_fun3 = colorRamp2(quantile(mat2All, c(0, 0.9999)), c("white", "red"))
EnrichedHeatmap(mat0All,name ="ELF5Dox",col=col_fun1)+
EnrichedHeatmap(mat1All,name ="FOXA1Dox",col=col_fun2)+
EnrichedHeatmap(mat2All,name ="ERDox",col=col_fun3)
dev.off()

#hmm... need motifs to confirm the hypothesis












#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
tss <- promoters(genes(txdb,columns=c("gene_id")), upstream=0, downstream=1)
#promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
geneid <- mapIds(org.Hs.eg.db, names(tss), "SYMBOL","ENTREZID")
tss=tss[!duplicated(geneid)]
names(tss)<-geneid[!duplicated(geneid)]






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
EnrichedHeatmap(mat1,   col = c("white", "darkred"),name ="FOXA1_diff_Dox_enriched",column_title = "FOXA1NoDox")+
  EnrichedHeatmap(mat2, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
  EnrichedHeatmap(mat3, col = c("white", "darkred"),name ="ER_diff_Dox_enriched",column_title = "ERDox")+
  EnrichedHeatmap(mat4, col = c("white", "darkblue"),name ="ERNoDox",column_title = "ERNoDox")
dev.off()














