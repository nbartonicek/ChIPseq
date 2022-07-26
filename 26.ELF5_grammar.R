#enhancers were downloaded from this paper: https://www.nature.com/articles/s41598-017-02257-3#Sec26


library(GenomicRanges)
library(ggplot2)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(JASPAR2018)
library(Biostrings)
library(ChIPpeakAnno)
library(pheatmap)
library(corrplot)

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
imageDir=paste0(resultsDir,"/figures/grammar/")
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")

system(paste0("mkdir -p ",imageDir))

liftover<-function(gr){
  gr.hg19=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg19ToHg38.over.chain"))
  seqlevelsStyle(gr.hg19) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg19, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}

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

#so the idea is to take all ELF5 peaks and assign TF motifs to them
#and then see if any of the grouping or overlaps work
#have to watch out to eliminate the massive overlaps so clean that out first
#but that can be eliminated as an option

peakWidth=100


ELF5=cleanGRs[["ELF5Dox"]]
gr1kb<-resize(ELF5,fix="center",width=peakWidth)

cleaning="clean"
if(cleaning=="clean"){
  data(HOT.spots)
  data(wgEncodeTfbsV3)
  hotGR <- liftover(reduce(unlist(HOT.spots)))
  wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
  removeOl <- function(.ele){
          ol <- findOverlaps(.ele, hotGR)
          if(length(ol)>0) .ele <- .ele[-unique(queryHits(ol))]
          .ele
  }
  gr1kb<-removeOl(gr1kb)
}

#get sequences
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
      
#now get all the transcription factors
TFs<-c("ELF5","FOXA1","ESR1","EOMES","ESRRA","GATA3","CDX2","TEAD4","POU5F1","GRHL2","FLI1","KLF4","RUNX1","SP1","PBX1","FOXO3","FOSB::JUNB","RARA","TFAP2C","XBP1")

i=75

file=paste0(resultsDir,"/Robjects/ELF5_grammar_",i,"_",peakWidth,".Rdata")
if(!file.exists(file)){
  for(TF in TFs){
    cat(TF)
    cat("\n")
    opts <- list()
    opts[["species"]] <- 9606
    opts[["name"]] <- TF
    #opts[["type"]] <- "SELEX"
    opts[["all_versions"]] <- TRUE
    PFMatrixList <- getMatrixSet(JASPAR2018, opts)
    pwm=PFMatrixList[[1]]
    hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )
    lengths<-sapply(hits,length)
    #adjust the GRanges
    values(gr1kb)[paste0("lengths_",TF)]=lengths
    starts<-sapply(hits,start)        
    values(gr1kb)[paste0("starts_",TF)]=sapply(starts,function(x){paste(x,collapse=",")})
  }

  save(gr1kb,file=paste0(resultsDir,"/Robjects/ELF5_grammar_",i,"_",peakWidth,".Rdata"))
}else{load(file)}

#annotate ELF5
gr1kb$class<-"nongenic"
gr1kb$resistance<-"other"


for(class in c("enhancer","genic","resistant","sensitive")){
      if(class=="enhancer"){
        #first import the enhancers
        enh<-read.table(paste0(annotationDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
        enh=enh[enh$type %in% "Distal_TE",]
        enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end),type=enh$type)
        grTarget<-liftover(enh.hg19)
        mat<-findOverlaps(gr1kb,grTarget)
        gr1kb$class[unique(queryHits(mat))]<-"enhancer"
      } else if (grepl("genic",class)){
        #TSS distribution
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        grTarget <- promoters(genes(txdb,columns=c("gene_id")), upstream=10000, downstream=10000)
        mat<-findOverlaps(gr1kb,grTarget)
        gr1kb$class[unique(queryHits(mat))]<-"genic"
      } else if (class=="resistant") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
        mat<-findOverlaps(gr1kb,grTarget)
        gr1kb$resistance[unique(queryHits(mat))]<-"resistant"
      } else if (class=="sensitive") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
        mat<-findOverlaps(gr1kb,grTarget)
        gr1kb$resistance[unique(queryHits(mat))]<-"sensitive"
      } 
}

#then calculate two metrices
#median distance
#count

dfLength<-values(gr1kb)[grepl("length",colnames(values(gr1kb)))]
dfLength


#d <- dist(dfLength) # euclidean distances between the rows
#fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
#fit # view results
#
## plot solution 
#x <- fit$points[,1]
#y <- fit$points[,2]
#pdf(paste0(imageDir,"MDS.pdf"),width=12,height=8)
#plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
#  main="Metric  MDS", type="n")
#dev.off()
dfLength=as.data.frame(dfLength)
colnames(dfLength)<-gsub("lengths_","",colnames(dfLength))
dfLength[dfLength>10]=10

row.names(dfLength)=paste0("r_",1:length(gr1kb))


annotation=data.frame(class=gr1kb$class,tamoxifen=gr1kb$resistance)
row.names(annotation)=row.names(dfLength)
pdf(paste0(imageDir,"pheatmap_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)
#pheatmap(dfLength,annotation_row=data.frame(class=ELF5$class,tamoxifen=ELF5$resistance))
pheatmap(dfLength,show_rownames = F,annotation_row=annotation,color = colorRampPalette(c("white", "blue", "darkblue"))(11))
dev.off()

#add the genic as a category
dfLength$genic=1*(annotation$class=="genic")
dfLength$tmxf=1*(annotation$tamoxifen=="resistant")

pdf(paste0(imageDir,"pheatmap_genicInfo_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)
#pheatmap(dfLength,annotation_row=data.frame(class=ELF5$class,tamoxifen=ELF5$resistance))
pheatmap(dfLength,show_rownames = F,annotation_row=annotation,color = colorRampPalette(c("white", "blue", "darkblue"))(11))
dev.off()


#add overlaps with FOXA1 and ER as a category
foxa1Peaks<-reduce(c(cleanGRsPeaks[["FOXA1Dox"]],cleanGRsPeaks[["FOXA1NoDox"]]))
dfLength$FOXA1_peaks<-countOverlaps(gr1kb,foxa1Peaks)

erPeaks<-reduce(c(cleanGRsPeaks[["ERDox"]],cleanGRsPeaks[["ERNoDox"]]))
dfLength$ER_peaks<-countOverlaps(gr1kb,erPeaks)

#plot correlations

pdf(paste0(imageDir,"pheatmap_genic_peaks_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)
#pheatmap(dfLength,annotation_row=data.frame(class=ELF5$class,tamoxifen=ELF5$resistance))
pheatmap(dfLength,show_rownames = F,annotation_row=annotation,color = colorRampPalette(c("white", "blue", "darkblue"))(11))
dev.off()


#add overlaps with FOXA1 and ER as a category
foxa1Peaks<-reduce(c(cleanGRsPeaks[["FOXA1Dox"]]))
dfLength$FOXA1_peaks<-countOverlaps(gr1kb,foxa1Peaks)

erPeaks<-reduce(c(cleanGRsPeaks[["ERDox"]]))
dfLength$ER_peaks<-countOverlaps(gr1kb,erPeaks)


pdf(paste0(imageDir,"pheatmap_genic_Doxpeaks_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)
#pheatmap(dfLength,annotation_row=data.frame(class=ELF5$class,tamoxifen=ELF5$resistance))
pheatmap(dfLength,show_rownames = F,annotation_row=annotation,color = colorRampPalette(c("white", "blue", "darkblue"))(11))
dev.off()

#plot correlations

pdf(paste0(imageDir,"correlations_genic_Doxpeaks_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
M=cor(dfLength)
res1 <- cor.mtest(dfLength, conf.level = .95)

corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = T)
dev.off()



pdf(paste0(imageDir,"correlations_Onlygenic_Doxpeaks_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
M=cor(dfLength[dfLength$genic==1,-which(colnames(dfLength)=="genic")])
res1 <- cor.mtest(dfLength[dfLength$genic==1,-which(colnames(dfLength)=="genic")], conf.level = .95)

corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = T)
dev.off()



pdf(paste0(imageDir,"correlations_Nongenic_Doxpeaks_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
M=cor(dfLength[dfLength$genic==0,-which(colnames(dfLength)=="genic")])
res1 <- cor.mtest(dfLength[dfLength$genic==1,-which(colnames(dfLength)=="genic")], conf.level = .95)

corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = T)
dev.off()





#plot only for genic
dfLength[dfLength>0]=1



pdf(paste0(imageDir,"correlations_genic_Doxpeaks_binary_",i,"_",cleaning,"_",peakWidth,".pdf"),width=12,height=12)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
M=cor(dfLength)
res1 <- cor.mtest(dfLength, conf.level = .95)

corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .7,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = T)
dev.off()
#plot only for genic







