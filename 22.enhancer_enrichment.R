#enhancers were downloaded from this paper: https://www.nature.com/articles/s41598-017-02257-3#Sec26


library(GenomicRanges)
library(ChIPpeakAnno)
library(ggplot2)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(rGREAT)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(hier.part)
library(fgsea)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(plyr)

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]

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
imageDir=paste0(resultsDir,"/figures/enhancers/")
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")

system(paste0("mkdir -p ",imageDir))

#first import the enhancers
enh<-read.table(paste0(annotationDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)

#types=c("TE","SE")
#type="SE"
#enh=enh[enh$type %in% c(paste0(c("Proximal_","Distal_"),type)),]

enh$type=mapvalues(enh$type, from=unique(enh$type),to=c("enh_proximal","enh_distal","superEnh_distal","superEnh_proximal"))
#enh$type=mapvalues(enh$type, from=unique(enh$type),to=c("enh_prox","enh_dist"))
enh=enh[order(enh$type),]

enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end),type=enh$type)

liftover<-function(gr){
  gr.hg19=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg19ToHg38.over.chain"))
  seqlevelsStyle(gr.hg19) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg19, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}

enh.hg38<-liftover(enh.hg19)

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


enh.hg38L<-split(enh.hg38,enh.hg38$type)
res<-list()
for(enh in names(enh.hg38L)){
  cat(enh)
  cat("\n")
  cat(countOverlaps(cleanGRsPeaks,enh.hg38L[[enh]])/length(enh.hg38L[[enh]]))
  cat("\n")
  res[[enh]]<-countOverlaps(cleanGRsPeaks,enh.hg38L[[enh]])/length(enh.hg38L[[enh]])
}
df<-do.call("rbind",res)
write.table(df,"../project_results/tables/enhancer_overlap.xls",row.names=T,quote=F,sep="\t")

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene



load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
#tssDE<-cleanGRs[["ELF5Dox"]]
#mat<-findOverlaps(tssDE,enh.hg38.Norm)
#tssDE<-enh.hg38.Norm[unique(subjectHits(mat))]


tssDE1=enh.hg38[grepl("enh",enh.hg38$type)]
tssDE1$type=gsub(".*_","",tssDE1$type)

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat1 = normalizeToMatrix(gr, tssDE1, value_column = "score", 
    extend = 0, mean_mode = "w0", w = 10,smooth=F)

tssDE2=enh.hg38[grepl("super",enh.hg38$type)]
tssDE2$type=gsub(".*_","",tssDE2$type)

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat2 = normalizeToMatrix(gr, tssDE2, value_column = "score", 
    extend = 0, mean_mode = "w0", w = 10,smooth=F)



pdf(paste0(imageDir,"enhancers_wholewidth_EFL5_sorted_start.pdf"),width=4,height=6)
EnrichedHeatmap(mat1, use_raster = TRUE,split = tssDE1$type, col = c("white", "darkblue"),name ="Enhancers",column_title = "Enhancers")
dev.off()

pdf(paste0(imageDir,"superenhancers_wholewidth_EFL5_sorted_start.pdf"),width=4,height=6)
EnrichedHeatmap(mat2, use_raster = TRUE,split = tssDE2$type, col = c("white", "darkblue"),name ="Super Enhancers",column_title = "Super Enhancers")
dev.off()




###  plot overlap

tssDE=enh.hg38[grepl("enh",enh.hg38$type)]
tssDE$type=gsub(".*_","",tssDE1$type)
tssDE<-resize(tssDE,1,fix="center")


sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)


sampleName="H3K4me3NoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="H3K4me3Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat6 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"enhancers_rasterized.pdf"),width=14,height=8)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
EnrichedHeatmap(mat2,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"),name ="FOXA1Dox",column_title = "FOXA1Dox")+
#EnrichedHeatmap(mat2.5, top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="FOXA1_Dox_enriched",column_title = "FOXA1_Dox_enriched")+
EnrichedHeatmap(mat3,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"), name ="ERNoDox",column_title = "ERNoDox") +
EnrichedHeatmap(mat4,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox") +
#EnrichedHeatmap(mat4.5, top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="ER_Dox_enriched",column_title = "ER_Dox_enriched")
EnrichedHeatmap(mat5,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"), name ="H3K4me3NoDox",column_title = "H3K4me3NoDox")+
EnrichedHeatmap(mat6,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"), name ="H3K4me3Dox",column_title = "H3K4me3Dox") 
#EnrichedHeatmap(mat6.5, col = c("white", "darkgreen"),name ="H3K4me3_Dox_enriched",column_title = "H3K4me3_Dox_enriched")
dev.off()



tssDE=enh.hg38[grepl("super",enh.hg38$type)]
tssDE$type=gsub(".*_","",tssDE1$type)
tssDE<-resize(tssDE,1,fix="center")


sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="FOXA1NoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat1 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat2 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat3 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat4 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)


sampleName="H3K4me3NoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="H3K4me3Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat6 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"superenhancers_rasterized.pdf"),width=14,height=8)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
EnrichedHeatmap(mat2,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"),name ="FOXA1Dox",column_title = "FOXA1Dox")+
#EnrichedHeatmap(mat2.5, top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="FOXA1_Dox_enriched",column_title = "FOXA1_Dox_enriched")+
EnrichedHeatmap(mat3,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"), name ="ERNoDox",column_title = "ERNoDox") +
EnrichedHeatmap(mat4,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox") +
#EnrichedHeatmap(mat4.5, top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="ER_Dox_enriched",column_title = "ER_Dox_enriched")
EnrichedHeatmap(mat5,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"), name ="H3K4me3NoDox",column_title = "H3K4me3NoDox")+
EnrichedHeatmap(mat6,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"), name ="H3K4me3Dox",column_title = "H3K4me3Dox") 
#EnrichedHeatmap(mat6.5, col = c("white", "darkgreen"),name ="H3K4me3_Dox_enriched",column_title = "H3K4me3_Dox_enriched")
dev.off()






