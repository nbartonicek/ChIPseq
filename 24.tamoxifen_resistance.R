library(rtracklayer)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(R.utils)
library(ChIPpeakAnno)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)

homedir="/share/ScratchGeneral/nenbar"
inPath=paste0(homedir,"/projects/Chris/project_results/GSE75372/")
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/tamoxifen/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste0(projectDir,"/scripts/")
  args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["id"]])){id = args$id} 
if (!is.null(args[["type"]])){type = args$type} 
sampleID="GSE75372"
ChIPs<-read.table(paste0(scriptsPath,sampleID,"_ids.txt"))
ChIPs$type<-c(rep("sensitive",1),rep("resistant",1))

peaks<-GRangesList()
for(i in 1:length(ChIPs$V2)){
  file=ChIPs$V2[i]
  sampleName=as.character(ChIPs$V1[i])
  cat(sampleName)
  cat("\n")
  gr<-import(paste0(inPath,file))
  #gr<-GRanges(seqnames=df$chr,IRanges(df$start,df$end))
  peaks[[sampleName]]<-gr
}

allPeaks<-unlist(peaks)
allPeaks<-reduce(allPeaks)
peakCount<-list()
for(peak in names(peaks)){
  peakCount[[peak]]<-1*(countOverlaps(allPeaks,peaks[[peak]])>0)
}
df<-do.call("cbind",peakCount)
df=as.data.frame(df)
#find those are present in resistant or absent in resistant
df$resistance=df[,1]-df[,2]
#there are 144571 peaks
#the sensitive (present in >1 out of 3)

sensitive=allPeaks[df$resistance==1]
resistant=allPeaks[df$resistance==-1]
#27227

#liftover
liftover<-function(gr){
  gr.hg19=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg19ToHg38.over.chain"))
  seqlevelsStyle(gr.hg19) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg19, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}


sensitiveLift<-liftover(sensitive)
resistantLift<-liftover(resistant)

#load the peaks and differentially expressed genes
#load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))
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

ERpeaks<-cleanGRsPeaks[c(1,2,3,8,9,10,11,12,13)]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)))
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
removeOl <- function(.ele){
        ol <- findOverlaps(.ele, hotGR)
        if(length(ol)>0) .ele <- .ele[-unique(queryHits(ol))]
        .ele
}
sensitiveLiftClean<-removeOl(sensitiveLift)
resistantLiftClean<-removeOl(resistantLift)
values(sensitiveLiftClean)<-NULL
values(resistantLiftClean)<-NULL


poolL<-list()
poolL[["sensitive"]] <- new("permPool", grs=GRangesList(wgEncodeTfbsV3), N=length(sensitiveLiftClean))
poolL[["resistant"]] <- new("permPool", grs=GRangesList(wgEncodeTfbsV3), N=length(resistantLiftClean))


#create a cool figure for enrichment
#get the foxa1 peaks that are sensitive / resistant
#centre elf5 peaks on them
#then plot diff fox and er


### first plot overlap

tamFox<-c(sensitiveLiftClean,resistantLiftClean)
values(tamFox)$type=rep(c("sensitive","resistant"),times=c(length(sensitiveLiftClean),length(resistantLiftClean)))

tssDE=tamFox

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

sampleName="FOXA1_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat2.5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
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

sampleName="ER_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat4.5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)
mat4.5[mat4.5<0]=0

#sampleName="H3K4me3NoDox"
#gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
#mat5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
#    extend = 1000, mean_mode = "w0", w = 50)
#
#sampleName="H3K4me3Dox"
#gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
#mat6 = normalizeToMatrix(gr, tssDE, value_column = "score", 
#    extend = 1000, mean_mode = "w0", w = 50)

#sampleName="H3K4me3_diff_Dox_enriched"
#gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
#mat6.5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
#    extend = 1000, mean_mode = "w0", w = 50)


pdf(paste0(imageDir,"FOXA1_peaks_EFL5_sorted_rasterized.pdf"),width=14,height=8)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
EnrichedHeatmap(mat2,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"),name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat2.5, top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="FOXA1_Dox_enriched",column_title = "FOXA1_Dox_enriched")+
EnrichedHeatmap(mat3,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"), name ="ERNoDox",column_title = "ERNoDox") +
EnrichedHeatmap(mat4,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox") +
EnrichedHeatmap(mat4.5, top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="ER_Dox_enriched",column_title = "ER_Dox_enriched")
#EnrichedHeatmap(mat5, col = c("white", "darkblue"), name ="H3K4me3NoDox",column_title = "H3K4me3NoDox")+
#EnrichedHeatmap(mat6, col = c("white", "darkred"), name ="H3K4me3Dox",column_title = "H3K4me3Dox") 
#EnrichedHeatmap(mat6.5, col = c("white", "darkgreen"),name ="H3K4me3_Dox_enriched",column_title = "H3K4me3_Dox_enriched")
dev.off()





### first plot overlap

tamFox<-c(sensitiveLiftClean,resistantLiftClean)
values(tamFox)$type=rep(c("sensitive","resistant"),times=c(length(sensitiveLiftClean),length(resistantLiftClean)))

tssDE=tamFox
sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
mat<-findOverlaps(tssDE,gr)
tssDE=tssDE[unique(queryHits(mat))]



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

sampleName="FOXA1_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat2.5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

sampleName="FOXA1_diff_Dox_depleted"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(-gr$score)

mat2.6 = normalizeToMatrix(gr, tssDE, value_column = "score", 
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

sampleName="ER_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat4.5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 1000, mean_mode = "w0", w = 50)

mat4.5[mat4.5<0]=0
#sampleName="H3K4me3NoDox"
#gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
#mat5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
#    extend = 1000, mean_mode = "w0", w = 50)
#
#sampleName="H3K4me3Dox"
#gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
#mat6 = normalizeToMatrix(gr, tssDE, value_column = "score", 
#    extend = 1000, mean_mode = "w0", w = 50)

#sampleName="H3K4me3_diff_Dox_enriched"
#gr=cleanGRsPeaks[[sampleName]]
#gr$score=log10(gr$score)
#mat6.5 = normalizeToMatrix(gr, tssDE, value_column = "score", 
#    extend = 1000, mean_mode = "w0", w = 50)


pdf(paste0(imageDir,"FOXA1_peaks_EFL5_overlapOnly_rasterized2.pdf"),width=14,height=8)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
EnrichedHeatmap(mat2,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"),name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat2.5, top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="FOXA1_Dox_enriched",column_title = "FOXA1_Dox_enriched")+
EnrichedHeatmap(mat2.6, top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="FOXA1_Dox_depleted",column_title = "FOXA1_Dox_depleted")+
EnrichedHeatmap(mat3,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"), name ="ERNoDox",column_title = "ERNoDox") +
EnrichedHeatmap(mat4,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox") +
EnrichedHeatmap(mat4.5, top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="ER_Dox_enriched",column_title = "ER_Dox_enriched")
#EnrichedHeatmap(mat6, col = c("white", "darkred"), name ="H3K4me3Dox",column_title = "H3K4me3Dox") 
#EnrichedHeatmap(mat6.5, col = c("white", "darkgreen"),name ="H3K4me3_Dox_enriched",column_title = "H3K4me3_Dox_enriched")
dev.off()

pdf(paste0(imageDir,"FOXA1_peaks_EFL5_overlapOnly_rasterized2.pdf"),width=16,height=8)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"),name ="FOXA1NoDox",column_title = "FOXA1NoDox")+
EnrichedHeatmap(mat2,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"),name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat2.5, top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="FOXA1_Dox_enriched",column_title = "FOXA1_Dox_enriched")+
EnrichedHeatmap(mat2.6, top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="FOXA1_Dox_depleted",column_title = "FOXA1_Dox_depleted")+
EnrichedHeatmap(mat3,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"), name ="ERNoDox",column_title = "ERNoDox") +
EnrichedHeatmap(mat4,   top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox") +
EnrichedHeatmap(mat4.5, top_annotation = HeatmapAnnotation(line = anno_enriched(yaxis_facing="left",gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkgreen"),name ="ER_Dox_enriched",column_title = "ER_Dox_enriched")
#EnrichedHeatmap(mat5, col = c("white", "darkblue"), name ="H3K4me3NoDox",column_title = "H3K4me3NoDox")+
#EnrichedHeatmap(mat6, col = c("white", "darkred"), name ="H3K4me3Dox",column_title = "H3K4me3Dox") 
#EnrichedHeatmap(mat6.5, col = c("white", "darkgreen"),name ="H3K4me3_Dox_enriched",column_title = "H3K4me3_Dox_enriched")
dev.off()

break()

sampleName="FOXA1_diff_Dox_enriched"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)

mat<-findOverlaps(gr,tssDE)
gr=gr[unique(queryHits(mat))]


peakAnno <- annotatePeak(gr, tssRegion=c(-2500, 2500), TxDb=txdb, annoDb="org.Hs.eg.db")
df=as.data.frame(peakAnno@anno)

write.table(df,file="../project_results/tables/FOXA1_doxEnriched_tamoxifenResistant_genes.xls",quote=F,sep="\t")

#Finding statistical difference between resistant and sensitive

#first load the diff peaks

load("/share/ScratchGeneral/nenbar/projects/Chris/project_results/Robjects/diff/FOXA1_diff.Rdata")
foxa1Diff<-experiment.DB

load("/share/ScratchGeneral/nenbar/projects/Chris/project_results/Robjects/diff/ER_diff.Rdata")
ERDiff<-experiment.DB

#Then find the sensitive and resistant

mat<-findOverlaps(foxa1Diff,sensitiveLiftClean)
foxa1DiffSens<-foxa1Diff[unique(queryHits(mat))]
mat<-findOverlaps(foxa1Diff,resistantLiftClean)
foxa1DiffRes<-foxa1Diff[unique(queryHits(mat))]
mat<-findOverlaps(ERDiff,sensitiveLiftClean)
ERDiffSens<-ERDiff[unique(queryHits(mat))]
mat<-findOverlaps(ERDiff,resistantLiftClean)
ERDiffRes<-ERDiff[unique(queryHits(mat))]


inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate/")
macsFiles<-list.files(inPath_separate,pattern="Peak",full.names=T)
macsFilesF<-macsFiles[grepl("FOXA1",macsFiles)]

results_separate_FoxA1<-GRangesList()
for(macsFile in macsFilesF){

  mypeaks <- read.delim(macsFile,header=F,stringsAsFactors=F)
  sampleName=gsub("_peaks.*Peak","",basename(macsFile))
  cat(sampleName)
  cat("\n")
  colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
  mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
  gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*",score=mypeaks$score,pval=mypeaks$log10.qvalue)
  names(gr)<-1:length(gr)
  results_separate_FoxA1[[sampleName]]<-gr

}


inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate/")
macsFiles<-list.files(inPath_separate,pattern="Peak",full.names=T)
macsFilesF<-macsFiles[grepl("ER",macsFiles)]

results_separate_ER<-GRangesList()
for(macsFile in macsFilesF){

  mypeaks <- read.delim(macsFile,header=F,stringsAsFactors=F)
  sampleName=gsub("_peaks.*Peak","",basename(macsFile))
  cat(sampleName)
  cat("\n")
  colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
  mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
  gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*",score=mypeaks$score,pval=mypeaks$log10.qvalue)
  names(gr)<-1:length(gr)
  results_separate_ER[[sampleName]]<-gr

}



#Count results_separate 



































