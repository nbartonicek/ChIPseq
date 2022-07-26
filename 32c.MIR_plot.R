

library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(org.Hs.eg.db)
library(circlize)
library(ChIPpeakAnno)
library(reshape2)

library(TFBSTools)
library(JASPAR2018)
library(Biostrings)
library(msa)

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
expressionDir=paste0(projectDir,"/annotation/")

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

########### plot or significantly changed repeats

repeats<-GRangesList()

load(paste0(annotationDir,"repeats_type.Rdata"))
typeGR<-grL
load(paste0(annotationDir,"repeats_class.Rdata"))
classGR<-grL

SINEs<-classGR[grepl("SINE",names(classGR))]
SINEs<-unlist(SINEs)
repeats[["SINEs"]]<-GRanges(paste0("chr",as.character(seqnames(SINEs))),IRanges(start=start(SINEs),end=end(SINEs)))

LINEs<-classGR[grepl("L2",names(classGR))]
LINEs<-unlist(LINEs)
repeats[["LINEs"]]<-GRanges(paste0("chr",as.character(seqnames(LINEs))),IRanges(start=start(LINEs),end=end(LINEs)))

SINEsMIRb<-typeGR[grepl("MIRb",names(typeGR))]
SINEsMIRb<-unlist(SINEsMIRb)
repeats[["MIRb"]]<-GRanges(paste0("chr",as.character(seqnames(SINEsMIRb))),IRanges(start=start(SINEsMIRb),end=end(SINEsMIRb)))

SINEsAluSx<-typeGR[grepl("AluSx",names(typeGR))]
SINEsAluSx<-unlist(SINEsAluSx)
repeats[["AluSx"]]<-GRanges(paste0("chr",as.character(seqnames(SINEsAluSx))),IRanges(start=start(SINEsAluSx),end=end(SINEsAluSx)))

Charlie1a<-typeGR[grepl("Charlie1a",names(typeGR))]
Charlie1a<-unlist(Charlie1a)
repeats[["Charlie1a"]]<-GRanges(paste0("chr",as.character(seqnames(Charlie1a))),IRanges(start=start(Charlie1a),end=end(Charlie1a)))



liftover<-function(gr,genome){
  if(genome=="hg38"){return(gr)} else {
    gr.old=gr
    ch = import.chain(paste0(projectDir,"/annotation/",genome,"ToHg38.over.chain"))
    seqlevelsStyle(gr.old) = "UCSC"  # necessary
    gr.hg38 = liftOver(gr.old, ch)
    gr.hg38  = unlist(gr.hg38 )
    genome(gr.hg38) = "hg38"
    return(gr.hg38)
  }
}


#first, how many of them overlap hotspots
data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)),"hg19")
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3,"hg19")   

regionWidth=1000

regions<-GRangesList()
for(i in 1:5){
  gr<-cleanGRs[[i]]
  values(gr)<-NULL
  gr<-resize(gr,regionWidth,fix="center")
  regions[[names(cleanGRs)[i]]]<-gr
}
for(sampleName in c("FOXA1_diff_Dox_enriched","FOXA1_diff_Dox_depleted","ER_diff_Dox_enriched","ER_diff_Dox_depleted")){
  gr=cleanGRsPeaks[[sampleName]]
  values(gr)<-NULL
  gr<-resize(gr,regionWidth,fix="center")

  regions[[sampleName]]<-gr
}
hotGR<-resize(hotGR,regionWidth,fix="center")
wgEncodeTfbsV3<-resize(wgEncodeTfbsV3,regionWidth,fix="center")

regions[["HOT"]]<-hotGR
regions[["wgEncodeTfbsV3"]]<-wgEncodeTfbsV3



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







resultsPercentWithRepeats<-list()
averageRepeats<-list()

for(region in names(regions)){
  cat(region)
  cat("\n")
  #c1<- sum(countOverlaps(regions[[region]],repeats[[1]])>0)
  #c2<- sum(countOverlaps(regions[[region]],repeats[[2]])>0)
  #c3<- sum(countOverlaps(regions[[region]],repeats[[3]])>0)
  #out<-c(c1,c2,c3)
  #names(out)<-names(repeats)
  #results[[region]]<-out/length(regions[[region]])
  resultsPercentWithRepeats[[region]]<-countOverlaps(repeats,regions[[region]])/length(regions[[region]])
  averageRepeats[[region]]<-sapply(repeats,function(x){sum(countOverlaps(x,regions[[region]]))})/length(regions[[region]])
}

df<-do.call("rbind",resultsPercentWithRepeats)
df
dfA<-do.call("rbind",averageRepeats)
dfA
dfA/df

write.table(df,paste0(resultsDir,"/tables/repeatOverlaps.xls"),quote=F,sep="\t")

temp=regions

#as a background use 1kb around all the encode TFBS (619678)
chrGR<-reduce(c(regions[["wgEncodeTfbsV3"]],regions[["ELF5Dox"]],regions[["FOXA1Dox"]],regions[["FOXA1NoDox"]],regions[["ERDox"]],regions[["ERNoDox"]]))
chrGR<-reduce(unlist(regions))
#chrGR<-reduce(c(regions[["ELF5Dox"]]))

repeatsClean<-endoapply(repeats,function(x){mat=findOverlaps(chrGR,x);x=x[unique(subjectHits(mat))];x=intersect(x,chrGR);x})
#chrGR=reduce(c(chrGR,unlist(repeatsClean)))


mat<-findOverlaps(cleanGRs[[1]],repeats[[3]])
mirOL<-cleanGRs[[1]][unique(queryHits(mat))]
mirOL

ELF5DE<-regions[["ELF5_DE"]]
mat<-findOverlaps(ELF5DE,resize(repeats[[3]],1000,fix="center"))
mirDEOL<-ELF5DE[unique(queryHits(mat))]
mirDEOL

mat<-findOverlaps(cleanGRs[[1]],ELF5DE)
ELF5DESum<-cleanGRs[[1]][unique(queryHits(mat))]
mat<-findOverlaps(ELF5DESum,repeats[[3]])
mirDEOLSummit<-ELF5DESum[unique(queryHits(mat))]
mirDEOLSummit

#so out of 28335 peaks, 1834 overlap, but only 499 summits are in MIR
#out of 119 ELF5 peaks in promotors of DE genes, 19 are in MIR (27 if MIRs extended by 1kb), but only 2 summits are in MIRs

#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[3]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(resize(repeats[[3]],5000+width(repeats[[3]]),fix="center"),cleanGRsPeaks[["ELF5Dox"]])

MIRsummits<-repeats[[3]][unique(queryHits(matSummits))]
MIRpeaks<-repeats[[3]][unique(queryHits(matPeaks))]
#jqums350-v4mw1jpo9r9d


sum(countOverlaps(cleanGRsPeaks[["ELF5Dox"]],repeats[["MIRb"]])>0)
sum(countOverlaps(cleanGRs[["ELF5Dox"]],repeats[["MIRb"]])>0)
repeats[["SINEs"]]

sum(countOverlaps(cleanGRs[["ELF5Dox"]],repeats[["SINEs"]])>0)


mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],repeats[["MIRb"]])
grQuery=cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat))]  
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, grQuery)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "ELF5"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
pwm=PFMatrixList[[1]]
i=75


hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )
lengths<-sapply(hits,length)
#adjust the GRanges
grQuery$lengths=lengths
starts<-sapply(hits,start)
grQuery$starts=starts
#forst repeat the gr by the amount of lengths
gr1kbL<-endoapply(split(grQuery,1:length(grQuery)),function(x){GenomicRanges::shift(rep(x,x$lengths),unlist(x$starts)-(width(x)/2-5))})
gr1kbNewShifted<-unlist(gr1kbL)

gr11L<-resize(gr1kbNewShifted,11,fix="center")


#gr1kbNew<-gr1kb[lengths>0]
#starts<-sapply(hits,start)
#starts=unlist(starts[lengths>0])
sum(countOverlaps(gr1kbL,repeats[["MIRb"]])>0)
sum(countOverlaps(gr11L,repeats[["MIRb"]])>0)

allRepeats<-GRanges(paste0("chr",as.character(seqnames(unlist(classGR)))),IRanges(start=start(unlist(classGR)),end=end(unlist(classGR))))
class<-unlist(classGR)$class
type<-unlist(classGR)$type

allRepeatsL<-split(allRepeats,class)
sum(countOverlaps(cleanGRsPeaks[["ELF5Dox"]],allRepeatsL)>0)
sort(countOverlaps(allRepeatsL,cleanGRsPeaks[["ELF5Dox"]]))

allRepeatsTypeL<-split(allRepeats,type)
sum(countOverlaps(cleanGRsPeaks[["ELF5Dox"]],allRepeatsTypeL)>0)
sort(countOverlaps(allRepeatsTypeL,cleanGRsPeaks[["ELF5Dox"]]))





seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kbNewShifted)
#hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )

names(seq) = paste0("SEQUENCEALL_", seq_along(seq))



#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[3]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[3]],cleanGRsPeaks[["ELF5Dox"]])

MIRsummits<-repeats[[3]][unique(queryHits(matSummits))]
MIRpeaks<-repeats[[3]][unique(queryHits(matPeaks))]

w=50
extend=500

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("dist_ELF5_ER_MirB_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()


#ER
matSummits<-findOverlaps(repeats[[3]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[3]],cleanGRsPeaks[["ELF5Dox"]])

MIRsummits<-repeats[[3]][unique(queryHits(matSummits))]
MIRpeaks<-repeats[[3]][unique(queryHits(matPeaks))]

w=50
extend=500

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]
ERNoDoxSummits<-resize(cleanGRs[["ERNoDox"]],1,fix="center")
ERNoDoxPeaks<-cleanGRsPeaks[["ERNoDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)
ERNoDoxSummits$score=log10(ERNoDoxSummits$score)
ERNoDoxPeaks$score=log10(ERNoDoxPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat4 = normalizeToMatrix(ERNoDoxSummits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat5 = normalizeToMatrix(ERNoDoxPeaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)

pdf(paste0(imageDir,paste0("ER_centred_dist_ELF5_ER_MirB_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1) +
  EnrichedHeatmap(mat4,name ="ERNoDoxSummits",column_title = "ERNoDoxSummits",col=col_fun1)+
  EnrichedHeatmap(mat5,name ="ERNoDoxPeaks",column_title = "ERNoDoxPeaks",col=col_fun1) +
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)
dev.off()

#ER test
matSummits<-findOverlaps(repeats[[3]],cleanGRs[["ERDox"]])
matPeaks<-findOverlaps(repeats[[3]],cleanGRsPeaks[["ERDox"]])

MIRsummits<-repeats[[3]][unique(queryHits(matSummits))]
MIRpeaks<-repeats[[3]][unique(queryHits(matPeaks))]

w=50
extend=500

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("ER_centred_ERoverlap_dist_ELF5_ER_MirB_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1) +
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)
dev.off()



#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[3]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[3]],cleanGRsPeaks[["ELF5Dox"]])

MIRsummits<-repeats[[3]][unique(queryHits(matSummits))]
MIRpeaks<-repeats[[3]][unique(queryHits(matPeaks))]

w=50
extend=500

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(gr11L, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("dist_ELF5_ER_MirB_motifs_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()


#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[4]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[4]],cleanGRsPeaks[["ELF5Dox"]])

AluSxsummits<-repeats[[4]][unique(queryHits(matSummits))]
AluSxpeaks<-repeats[[4]][unique(queryHits(matPeaks))]

w=50
extend=4000

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(gr11L, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("dist_ELF5_ER_AluSx_motifs_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()



#then find the distance of ELF5 summits and peaks to MIR

w=50
extend=2000

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, repeats[[4]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, repeats[[4]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, repeats[[4]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, repeats[[4]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("dist_ELF5_ER_AluSx_Total_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()






#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[4]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[4]],cleanGRsPeaks[["ELF5Dox"]])

AluSxsummits<-repeats[[4]][unique(queryHits(matSummits))]
AluSxpeaks<-repeats[[4]][unique(queryHits(matPeaks))]

w=50
extend=500

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("dist_ELF5_ER_AluSx_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()




#MIR_DE
ELF5DE<-regions[["ELF5_DE"]]
mat<-findOverlaps(ELF5DE,resize(repeats[[3]],1000,fix="center"))
mirDEOL<-ELF5DE[unique(queryHits(mat))]
mirDEOL

mat<-findOverlaps(cleanGRs[[1]],ELF5DE)
ELF5DESum<-cleanGRs[[1]][unique(queryHits(mat))]
mat<-findOverlaps(ELF5DESum,repeats[[3]])
mirDEOLSummit<-ELF5DESum[unique(queryHits(mat))]
mirDEOLSummit

#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[3]],ELF5DE)
matPeaks<-findOverlaps(repeats[[3]],ELF5DE)

MIRsummits<-repeats[[3]][unique(queryHits(matSummits))]
MIRpeaks<-repeats[[3]][unique(queryHits(matPeaks))]


w=50
extend=500

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("dist_ELF5DE_ER_MIRb_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()



################ FINAL!
matSummits<-findOverlaps(repeats[[3]],cleanGRsPeaks[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[3]],cleanGRsPeaks[["ELF5Dox"]])

MIRsummits<-repeats[[3]][unique(queryHits(matSummits))]
MIRpeaks<-repeats[[3]][unique(queryHits(matPeaks))]


w=50
extend=2000

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, MIRpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("final_MIRb_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()




#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[4]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[4]],cleanGRsPeaks[["ELF5Dox"]])

AluSxsummits<-repeats[[4]][unique(queryHits(matSummits))]
AluSxpeaks<-repeats[[4]][unique(queryHits(matPeaks))]

w=50
extend=2000

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, AluSxpeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("final_AluSx_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()



#then find the distance of ELF5 summits and peaks to MIR
matSummits<-findOverlaps(repeats[[5]],cleanGRs[["ELF5Dox"]])
matPeaks<-findOverlaps(repeats[[5]],cleanGRsPeaks[["ELF5Dox"]])

Charlie1asummits<-repeats[[5]][unique(queryHits(matSummits))]
Charlie1apeaks<-repeats[[5]][unique(queryHits(matPeaks))]

w=50
extend=2000

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, Charlie1apeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, Charlie1apeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, Charlie1apeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, Charlie1apeaks, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("final_Charlie1_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()



w=50
extend=1000

ELF5Summits<-resize(cleanGRs[["ELF5Dox"]],1,fix="center")
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-resize(cleanGRs[["ERDox"]],1,fix="center")
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Summits$score=log10(ELF5Summits$score)
ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)

mat0= normalizeToMatrix(ELF5Summits, repeats[[5]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, repeats[[5]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(ERSummits, repeats[[5]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, repeats[[5]], value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("final_Charlie1Test_Summits_log_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="ELF5Summits",column_title = "ELF5Summits",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="ERSummits",column_title = "ERSummits",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()

########### Final reverse 



#MIR_DE
ELF5DE<-regions[["ELF5_DE"]]
ELF5Summits<-cleanGRs[["ELF5Dox"]]
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]



mat<-findOverlaps(ELF5DE,ELF5Summits)
ELF5DESummits<-ELF5Summits[unique(subjectHits(mat))]
ELF5DESummits


w=20
extend=1000

ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-cleanGRs[["ERDox"]]
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)
MIRpeaks<-repeats[[3]]
MIRpeaks$score=5
AluSxPeaks<-repeats[[4]]
AluSxPeaks$score=5
CharliePeaks<-repeats[[5]]
CharliePeaks$score=5

mat0= normalizeToMatrix(MIRpeaks, ELF5DESummits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, ELF5DESummits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(AluSxPeaks, ELF5DESummits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, ELF5DESummits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("MIRb_centered_ELF5DE_AluSx_Charlie_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="MIRb",column_title = "MIRb",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="AluSx",column_title = "AluSx",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()


#MIR_DE
ELF5DE<-regions[["ELF5_DE"]]
ELF5Summits<-cleanGRs[["ELF5Dox"]]
ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]



mat<-findOverlaps(ELF5DE,ELF5Summits)
ELF5DESummits<-ELF5Summits[unique(subjectHits(mat))]
ELF5DESummits


w=20
extend=1000

ELF5Peaks<-cleanGRsPeaks[["ELF5Dox"]]
ERSummits<-cleanGRs[["ERDox"]]
ERPeaks<-cleanGRsPeaks[["ERDox"]]

ELF5Peaks$score=log10(ELF5Peaks$score)
ERSummits$score=log10(ERSummits$score)
ERPeaks$score=log10(ERPeaks$score)
MIRpeaks<-repeats[[3]]
MIRpeaks$score=5

mat<-findOverlaps(MIRpeaks,resize(ELF5Summits,2000,fix="center"))
MIRpeaks=MIRpeaks[unique(queryHits(mat))]

AluSxPeaks<-repeats[[4]]
AluSxPeaks$score=5
CharliePeaks<-repeats[[5]]
CharliePeaks$score=5

mat0= normalizeToMatrix(MIRpeaks, ELF5Summits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat1 = normalizeToMatrix(ELF5Peaks, ELF5Summits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat2 = normalizeToMatrix(AluSxPeaks, ELF5Summits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)
mat3 = normalizeToMatrix(ERPeaks, ELF5Summits, value_column = "score", 
    extend = extend, mean_mode = "w0", w = w)


pdf(paste0(imageDir,paste0("MIRb_centered_ELF5All_AluSx_Charlie_",w,"_",extend,".pdf")),width=8,height=6)
col_fun1 = colorRamp2(quantile(mat0, c(0, 0.99)), c("white", "red"))
  EnrichedHeatmap(mat0,name ="MIRb",column_title = "MIRb",col=col_fun1)+
  EnrichedHeatmap(mat1,name ="ELF5Peaks",column_title = "ELF5Peaks",col=col_fun1)+
  EnrichedHeatmap(mat2,name ="AluSx",column_title = "AluSx",col=col_fun1)+
  EnrichedHeatmap(mat3,name ="ERPeaks",column_title = "ERPeaks",col=col_fun1)
dev.off()

######################## sequence consensus #####################

mirB<-repeats[["MIRb"]]
mat<-findOverlaps(mirB,cleanGRsPeaks[["ELF5Dox"]])

mirBELF5<-mirB[unique(queryHits(mat))]

 
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, mirBELF5)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/MIRbELF5.fasta"))




mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],repeats[["MIRb"]])
grQuery=repeats[["MIRb"]][unique(subjectHits(mat))]  
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, grQuery)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "ELF5"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
pwm=PFMatrixList[[1]]
i=75


hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )
lengths<-sapply(hits,length)
#adjust the GRanges
grQuery$lengths=lengths
starts<-sapply(hits,start)
grQuery$starts=starts
#forst repeat the gr by the amount of lengths
gr1kbL<-endoapply(split(grQuery,1:length(grQuery)),function(x){GenomicRanges::shift(rep(x,x$lengths),unlist(x$starts)-(width(x)/2-5))})
gr1kbNewShifted<-unlist(gr1kbL)

gr11L<-resize(gr1kbNewShifted,11,fix="center")





seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr11L)
motifs=lapply(seq,as.character)
motifsUl<-unlist(motifs)
pdf(paste0(imageDir,"MIRbELF5_logo.pdf"),width=6,height=3)
p<-ggseqlogo( motifsUl,method="bits")
print(p)
dev.off()



