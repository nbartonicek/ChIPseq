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

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
chrGR<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))
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
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
grammarDir=paste0(resultsDir,"/grammar/")
system(paste0("mkdir -p ",imageDir))


sampleFile<-paste0(grammarDir,"ChIP_experiments.csv")
samples<-read.csv(sampleFile,header=T)
samples<-samples[,1:8]

#for each sample in 
#bed.gz download,unzip,load
#bed download, load
#bedGraph.gz unzip, load?
#bw download
#narrowPeak download, import
#narrowPeak.gz download, unzip, import
#txt gz download, fix, load


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

#for(type in c("bed.gz","bed")){
#for(type in c("narrowPeak.gz","narrowPeak")){
grL<-GRangesList()
load(paste0(grammarDir,"all.Rdata"))

#eliminate the encode replicates
#eliminate=c(72,74,76,78,80,82,84,86,88,90,92)
#grL=grL[-eliminate]

#for(type in c("bed.gz","bed","bw","narrowPeak.gz","narrowPeak","txt.gz")){
#for(type in c("bedGraph","txt.gz","bed.gz","bed","bw","narrowPeak.gz","narrowPeak")){
for(type in c("GRCh38.bed.gz")){
#for(type in c("bed.gz")){
  sampleData<-samples[samples$format %in% type,]
  sampleNames<-paste0(sampleData$TF,".",sampleData$treatment,".p",sampleData$ProjectID,".r",sampleData$replicate)
  for(id in 1:length(sampleNames)){
    cat(sampleNames[id])
    cat("\n")
    url=as.character(sampleData[id,"download"])
    fileName=basename(url)
    outDir=paste0(grammarDir,sampleData[id,"ProjectID"],"/")
    outFile=paste0(outDir,sampleNames[id],".",type)
    system(paste0("mkdir -p ",outDir))
    #download function
    if(grepl("bedGraph",type)){
      inDir=paste0(grammarDir,"GSE60272/")
      files<-list.files(inDir,full.names=T)
      inFile=files[grepl(as.character(sampleData[id,"download"]),files)]
      mypeaks=fread(inFile)
      mypeaks=as.data.frame(mypeaks)
      colnames(mypeaks)[1:4] <- c("chrom", "chromStart", "chromEnd", "score")
      mypeaks=mypeaks[mypeaks$score>50,]
      gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*") 
      gr=reduce(gr)
    } else {
      system(paste0("wget ",url," -O ",outFile))
    }
    #unzip
    if(grepl("gz",outFile)){
      system(paste0("gunzip ",outFile))
      outFile=gsub(".gz","",outFile)
    }
    #import
    if(grepl("bed$",outFile)&grepl("nnn|41820",sampleData[id,"download"])){
      mypeaks <- read.table(outFile,header=F,sep="\t",encoding="latin1")
      colnames(mypeaks)[1:3] <- c("chrom", "chromStart", "chromEnd")
      mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
      gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*") 
    }
    if(grepl("bed$",outFile)&!(grepl("nnn|41820",sampleData[id,"download"]))){
      gr=import(outFile)
    }
    if(grepl("narrowPeak",outFile)){
      mypeaks <- read.delim(outFile,header=F,stringsAsFactors=F)
      colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
      mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
      gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*") 
    }
    if(grepl("txt",outFile)&grepl("GSE72848",sampleData[id,"download"])){
      mypeaks <- read.table(outFile,header=T,sep="\t")
      colnames(mypeaks)[1:3] <- c("chrom", "chromStart", "chromEnd")
      mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
      mypeaks<-mypeaks[mypeaks[,paste0("Is",as.character(sampleData[id,"TF"]))]=="Y",]
      gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*") 
    }
    if(grepl("txt",outFile)&(!grepl("GSE72848",sampleData[id,"download"]))){
      mypeaks <- read.table(outFile,header=F,sep="\t",encoding="latin1",skip=27)
      colnames(mypeaks)[1:10] <- c("chrom", "chromStart", "chromEnd", "length","abs_summit", "pileup", "fold.enrichment","log10.pvalue", "log10.qvalue","name")
      mypeaks<-mypeaks[!grepl("_",mypeaks$chrom),]
      gr<-GRanges(seqnames=mypeaks$chrom,IRanges(start=mypeaks$chromStart,end=mypeaks$chromEnd),strand="*") 
    }
    if(grepl("bw",outFile)){
      temp=import(outFile)
      gr=temp[temp$score>=(-log10(0.05))]
      gr=reduce(gr)
    }
    genome=as.character(sampleData[id,"genome"])
    gr<-liftover(gr,genome)
    values(gr)=NULL
    save(gr,file=paste0(grammarDir,sampleNames[id],".Rdata"))
    grL[[sampleNames[id]]]<-gr
  }
}
save(grL,file=paste0(grammarDir,"all.Rdata"))

#convert to granges, liftover, save as Rdata

outBedDir<-paste0(resultsDir,"/grammarBeds/")
system(paste0("mkdir -p ",outBedDir))

for(sampleName in names(grL)){
  cat(sampleName)
  cat("\n")
  export(grL[[sampleName]],paste0(outBedDir,sampleName,".bed"))
}

grLShort=grL[sapply(grL,length)<100000]

grLShort=grLShort[grepl("vehicle|shCTL_Veh",names(grLShort))]



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
  #values(grPos)<-NULL
  #values(grNeg)<-NULL
  cleanGRsPeaks[[paste0(sampleName,"_Dox_enriched")]]<-grPos
  cleanGRsPeaks[[paste0(sampleName,"_Dox_depleted")]]<-grNeg
}

cleanGRsPeaksShort<-GRangesList()
for(sampleName in names(cleanGRsPeaks)[c(1,8,9,12,13)]){
  gr=cleanGRsPeaks[[sampleName]]
  values(gr)<-NULL
  grLShort[[sampleName]]=gr
}

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene




allPeaks<-reduce(unlist(c(grLShort)))
for(sampleName in as.character(names(grLShort))){
  cat(sampleName)
  cat("\n")
  values(allPeaks)[sampleName]<-countOverlaps(allPeaks,grLShort[sampleName])
}

df<-as.data.frame(values(allPeaks))

#first, how many of them overlap hotspots
data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)),"hg19")
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3,"hg19")   

#what percent of peaks overlap HOT regions
countsL<-list()
for(i in 1:length(grLShort)){
  mat<-findOverlaps(grLShort[[i]],hotGR)
  countsL[[i]]<-length(unique(queryHits(mat)))
}
percentHotOverlap<-unlist(countsL)/sapply(grLShort,length)

countOverlaps(hotGR,grLShort)

pdf(paste0(imageDir,"correlations_all_first_vehicle.pdf"),width=20,height=20)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
M=cor(df)
res1 <- cor.mtest(df, conf.level = .95)

corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .5,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = T)
dev.off()


pdf(paste0(imageDir,"correlations_all_first_estradiol.pdf"),width=20,height=20)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
M=cor(df)
res1 <- cor.mtest(df, conf.level = .95)

corrplot(M, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .5,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         p.mat = res1$p, sig.level = 0.01, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag = T)
dev.off()



#now clean up the data

gr1kb=allPeaks
gr1kb$class<-"nongenic"
gr1kb$resistance<-"other"

for(class in c("enhancer","genic","DE","resistant","sensitive")){
      if(class=="enhancer"){
        #first import the enhancers
        enh<-read.table(paste0(annotationDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
        enh=enh[enh$type %in% "Distal_TE",]
        enh.hg19<-GRanges(seqnames=enh$chrom,IRanges(start=enh$star,end=enh$end),type=enh$type)
        grTarget<-liftover(enh.hg19,"hg19")
        mat<-findOverlaps(gr1kb,grTarget)
        gr1kb$class[unique(queryHits(mat))]<-"enhancer"
      } else if (grepl("genic",class)){
        #TSS distribution
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        grTarget <- promoters(genes(txdb,columns=c("gene_id")), upstream=5000, downstream=5000)
        mat<-findOverlaps(gr1kb,grTarget)
        gr1kb$class[unique(queryHits(mat))]<-"genic"
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
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        grTarget <- tssDE
        mat<-findOverlaps(gr1kb,grTarget)
        gr1kb$class[unique(queryHits(mat))]<-"DE"
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



#Now for each class plot correlations

#plot correlations
for(dataClass in unique(gr1kb$class)){
  cat(dataClass)
  temp=gr1kb[gr1kb$class==dataClass]
  df=as.data.frame(values(temp))
  df=df[,-c(98,99)]
  df=df[,colSums(df)>0]
  pdf(paste0(imageDir,"correlations_",dataClass,"_vehicle.pdf"),width=30,height=30)

  col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
  M=cor(df)
  res1 <- cor.mtest(df, conf.level = .95)

  corrplot(M, method = "color", col = col(200),
           type = "upper", order = "hclust", number.cex = .5,
           addCoef.col = "black", # Add coefficient of correlation
           tl.col = "black", tl.srt = 90, # Text label color and rotation
           # Combine with significance
           p.mat = res1$p, sig.level = 0.01, insig = "blank", 
           # hide correlation coefficient on the principal diagonal
           diag = T)
  dev.off()

}


#plot correlations
dataClass="DE"

temp=gr1kb[gr1kb$class==dataClass]
df=as.data.frame(values(temp))
df=df[,-c(98,99)]
df=df[,colSums(df)>0]
M=cor(df)
res1 <- cor.mtest(df, conf.level = .95)

dataClass="enhancer"
temp=gr1kb[gr1kb$class==dataClass]
df=as.data.frame(values(temp))
df=df[,-c(98,99)]
df=df[,colnames(df) %in% colnames(M)]
M2=cor(df)
res2 <- cor.mtest(df, conf.level = .95)

totalM=M-M2

pdf(paste0(imageDir,"correlations_diff_genic_enhancer_vehicle.pdf"),width=30,height=30)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(totalM, method = "color", col = col(200),
         type = "upper", order = "hclust", number.cex = .5,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col = "black", tl.srt = 90, # Text label color and rotation
         # Combine with significance
         # hide correlation coefficient on the principal diagonal
         diag = T)
dev.off()




































































