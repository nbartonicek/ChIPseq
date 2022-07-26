

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
imageDir=paste0(resultsDir,"/figures/repeats/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")

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



targets<-GRangesList()
for(class in c("enhancer","superenhancer","HOT","DE","resistant","sensitive")){
      cat(".")
      if(class=="enhancer"){
        #first import the enhancers
        enh<-read.table(paste0(expressionDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
        enh=enh[enh$type %in% "Distal_TE",]
        enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end))
        grTarget<-liftover(enh.hg19,"hg19")
        mat<-findOverlaps(grTarget,regions[[1]])
        regions[[paste0("ELF5_",class)]]<-regions[[1]][unique(subjectHits(mat))]
      } else if (grepl("superenhancer",class)){
        enh<-read.table(paste0(expressionDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
        enh=enh[enh$type %in% "Distal_SE",]
        enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end))
        grTarget<-liftover(enh.hg19,"hg19")
        mat<-findOverlaps(grTarget,regions[[1]])
        regions[[paste0("ELF5_",class)]]<-regions[[1]][unique(subjectHits(mat))]


      } else if (grepl("HOT",class)){
        #TSS distribution
                
        mat<-findOverlaps(hotGR,regions[[1]])
        regions[[paste0("ELF5_",class)]]<-regions[[1]][unique(subjectHits(mat))]
        regions[[paste0("ELF5_non",class)]]<-regions[[1]][-unique(subjectHits(mat))]


      } else if (grepl("DE",class)){
        #TSS distribution
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        tss <- promoters(genes(txdb,columns=c("gene_id")), upstream=5000, downstream=5000)
        #promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
        geneid <- mapIds(org.Hs.eg.db, names(tss), "SYMBOL","ENTREZID")
        tss=tss[!duplicated(geneid)]
        names(tss)<-geneid[!duplicated(geneid)]

        #DE genes
        expression<-read.table(paste0(rnaseqDir,"MCF7.plusDox_vs_minusDox.limma.txt"),header=T,sep="\t",stringsAsFactors=F)
        DEgenes<-unique(expression$hgnc_symbol[expression$adj.P.Val<0.05&abs(expression$logFC)>1])
        tssDE<-tss[names(tss) %in% DEgenes]
        grTarget <- tssDE
        values(grTarget)<-NULL
        mat<-findOverlaps(grTarget,regions[[1]])
        regions[[paste0("ELF5_",class)]]<-regions[[1]][unique(subjectHits(mat))]

      } else if (class=="resistant") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
        values(grTarget)<-NULL
        mat<-findOverlaps(grTarget,regions[[1]])
        regions[[paste0("ELF5_",class)]]<-regions[[1]][unique(subjectHits(mat))]

      } else if (class=="sensitive") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
        values(grTarget)<-NULL
        mat<-findOverlaps(grTarget,regions[[1]])
        regions[[paste0("ELF5_",class)]]<-regions[[1]][unique(subjectHits(mat))]

      } 
}



########### plot or significantly changed repeats

load(paste0(annotationDir,"repeats_type.Rdata"))
typeGR<-grL
load(paste0(annotationDir,"repeats_class.Rdata"))
classGR<-grL
#load(paste0(annotationDir,"repeats.Rdata"))

SINEs<-classGR[grepl("SINE",names(classGR))]
SINEs<-unlist(SINEs)
SINEs<-GRanges(paste0("chr",as.character(seqnames(SINEs))),IRanges(start=start(SINEs),end=end(SINEs)))

LINEs<-classGR[grepl("L2",names(classGR))]
LINEs<-unlist(LINEs)
LINEs<-GRanges(paste0("chr",as.character(seqnames(LINEs))),IRanges(start=start(LINEs),end=end(LINEs)))

load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))


####### first plot all repeats

tssDE<-cleanGRs[["ELF5Dox"]]

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEs"
gr=SINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ELF5_SINE_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEs",column_title = "SINEs")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()

#######  plot all repeats but narrow peak

tssDE<-cleanGRs[["ELF5Dox"]]

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 200, mean_mode = "w0", w = 50)

sampleName="SINEs"
gr=SINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,400,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 200, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,400,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 200, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ELF5_SINE_LINE_400.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEs",column_title = "SINEs")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()


########### plot diff stuff, but overlapped with ELF5+FOXA1 only


tssDE<-cleanGRs[["ELF5Dox"]]

SINEsMIRb<-typeGR[grepl("MIRb",names(typeGR))]
SINEsMIRb<-unlist(SINEsMIRb)
SINEsMIRb<-GRanges(paste0("chr",as.character(seqnames(SINEsMIRb))),IRanges(start=start(SINEsMIRb),end=end(SINEsMIRb)))


sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEsMIRb"
gr=SINEsMIRb
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ELF5_MirB_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEsMIRb",column_title = "SINEsMIRb")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()


########### ELF5 DE


tssDE<-cleanGRs[["ELF5Dox"]]

SINEsMIRb<-typeGR[grepl("MIRb",names(typeGR))]
SINEsMIRb<-unlist(SINEsMIRb)
SINEsMIRb<-GRanges(paste0("chr",as.character(seqnames(SINEsMIRb))),IRanges(start=start(SINEsMIRb),end=end(SINEsMIRb)))

mat<-findOverlaps(tssDE,regions[["ELF5_DE"]])
tssDE=tssDE[unique(queryHits(mat))]

sampleName="ELF5Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEsMIRb"
gr=SINEsMIRb
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ELF5DE_MirB_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="ELF5Dox",column_title = "ELF5Dox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEsMIRb",column_title = "SINEsMIRb")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()


####### FOXA1 peaks

tssDE<-cleanGRs[["FOXA1Dox"]]

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEs"
gr=SINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_FOXA1Dox_SINE_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEs",column_title = "SINEs")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()

####### Foxa1 short

tssDE<-cleanGRs[["FOXA1Dox"]]

sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 200, mean_mode = "w0", w = 50)

sampleName="SINEs"
gr=SINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,400,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 200, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,400,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 200, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_FOXA1Dox_SINE_LINE_400.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, use_raster = TRUE, col = c("white", "darkred"), name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat1, use_raster = TRUE, col = c("white", "darkblue"), name ="SINEs",column_title = "SINEs")+
EnrichedHeatmap(mat2, use_raster = TRUE, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()

####### Foxa1 MIRb


tssDE<-cleanGRs[["FOXA1Dox"]]

SINEsMIRb<-typeGR[grepl("MIRb",names(typeGR))]
SINEsMIRb<-unlist(SINEsMIRb)
SINEsMIRb<-GRanges(paste0("chr",as.character(seqnames(SINEsMIRb))),IRanges(start=start(SINEsMIRb),end=end(SINEsMIRb)))


sampleName="FOXA1Dox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEsMIRb"
gr=SINEsMIRb
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_FOXA1Dox_MirB_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, use_raster = TRUE, col = c("white", "darkred"), name ="FOXA1Dox",column_title = "FOXA1Dox")+
EnrichedHeatmap(mat1, use_raster = TRUE, col = c("white", "darkblue"), name ="SINEsMIRb",column_title = "SINEsMIRb")+
EnrichedHeatmap(mat2, use_raster = TRUE, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()


####### ER peaks

tssDE<-cleanGRs[["ERDox"]]

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEs"
gr=SINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ERDox_SINE_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEs",column_title = "SINEs")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()

####### ER short

tssDE<-cleanGRs[["ERDox"]]

sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 200, mean_mode = "w0", w = 50)

sampleName="SINEs"
gr=SINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,400,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 200, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,400,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 200, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ERDox_SINE_LINE_400.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEs",column_title = "SINEs")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()


####### ER NoDox peaks

tssDE<-cleanGRs[["ERNoDox"]]

sampleName="ERNoDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEs"
gr=SINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ERNoDox_SINE_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, col = c("white", "darkred"), name ="ERNoDox",column_title = "ERNoDox")+
EnrichedHeatmap(mat1, col = c("white", "darkblue"), name ="SINEs",column_title = "SINEs")+
EnrichedHeatmap(mat2, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()



####### ER MIRb


tssDE<-cleanGRs[["ERDox"]]

SINEsMIRb<-typeGR[grepl("MIRb",names(typeGR))]
SINEsMIRb<-unlist(SINEsMIRb)
SINEsMIRb<-GRanges(paste0("chr",as.character(seqnames(SINEsMIRb))),IRanges(start=start(SINEsMIRb),end=end(SINEsMIRb)))


sampleName="ERDox"
gr=cleanGRsPeaks[[sampleName]]
gr$score=log10(gr$score)
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score", 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="SINEsMIRb"
gr=SINEsMIRb
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat1 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

sampleName="LINEs"
gr=LINEs
names(gr)<-NULL
mat<-findOverlaps(gr,resize(tssDE,4000,fix="center"))
gr=gr[unique(queryHits(mat))]
mat2 = normalizeToMatrix(gr, tssDE, 
    extend = 2000, mean_mode = "w0", w = 50)

pdf(paste0(imageDir,"dist_ERDox_MirB_LINE.pdf"),width=5,height=6)
EnrichedHeatmap(mat0, use_raster = TRUE, col = c("white", "darkred"), name ="ERDox",column_title = "ERDox")+
EnrichedHeatmap(mat1, use_raster = TRUE, col = c("white", "darkblue"), name ="SINEsMIRb",column_title = "SINEsMIRb")+
EnrichedHeatmap(mat2, use_raster = TRUE, col = c("white", "darkgreen"), name ="LINEs",column_title = "LINEs")
dev.off()




