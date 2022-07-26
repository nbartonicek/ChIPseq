

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


oddsRatio<-function(region,regionDatabase){

  #for a given nucleotide in region, what are the odds it has an epigenetic mark
  totalNuc<-sum(width(reduce(region)))
  reducedsnpDatabase=reduce(regionDatabase)

  #mat<-findOverlaps(region,reducedsnpDatabase)
  #total hits is the length of intersect
  strand(region)="*"
  strand(reducedsnpDatabase)="*"

  totalHits<-sum(width(intersect(reduce(region),reducedsnpDatabase)))

  #totalHits<-length(unique(subjectHits(mat)))
  firstRatio<-totalHits/(totalNuc-totalHits)

  nonRegion<-sum(as.numeric(width(chrGR)))-totalNuc
  genomeHits<-sum(as.numeric(width(reducedsnpDatabase)))-totalHits
  secondRatio<-genomeHits/(nonRegion-genomeHits)

  odds<-firstRatio/secondRatio
  se<-sqrt(1/totalNuc+1/(totalNuc-totalHits)+1/(nonRegion-genomeHits)+1/genomeHits)
  #cat(odds)
  #cat("\n")
  #cat(se)
  #cat("\n")
  result<-list()
  result[[1]]=odds
  result[[2]]=se


        #p-value
        SNPenrichment <-
     matrix(c(as.numeric(totalHits), as.numeric(totalNuc), as.numeric(genomeHits), as.numeric(nonRegion)),
            nrow = 2,
            dimnames = list(Regions = c("WithSNP", "Total"),
                            Genome = c("WithSNP", "Total")))
        result[[3]]=chisq.test(SNPenrichment)$p.value

        return(result)

}

resultsOdds<-list()
resultsSE<-list()
resultsP<-list()

regionsTemp<-regions[c("wgEncodeTfbsV3","FOXA1Dox","FOXA1NoDox","ERDox","ERNoDox","ELF5Dox","ELF5_nonHOT","ELF5_HOT","ELF5_enhancer","ELF5_superenhancer","ELF5_DE")]
#regions<-regions[c("ELF5Dox","ELF5_nonHOT","ELF5_enhancer","ELF5_superenhancer","ELF5_HOT","ELF5_DE")]

for(i in 1:length(regionsTemp)){
        #cat(i)

        temp<-list()
        for(j in 1:length(repeatsClean)){
                db=repeatsClean[[j]]
                db<-reduce(db)
                strand(db)<-"*"

                stats<-unlist(oddsRatio(regionsTemp[[i]],db))
                #names(stats)<-c("oddsRatio","SE")
                temp[[names(repeatsClean)[j]]]<-stats     
                cat(".")
        }
        resultsOdds[[names(regionsTemp)[i]]]<-do.call("cbind",temp)[1,]
        resultsSE[[names(regionsTemp)[i]]]<-do.call("cbind",temp)[2,]
        resultsP[[names(regionsTemp)[i]]]<-do.call("cbind",temp)[3,]
}

dfOdds<-do.call("rbind",resultsOdds)
dfSE<-do.call("rbind",resultsSE)
dfP<-do.call("rbind",resultsP)

dfOdds=as.data.frame(dfOdds)
dfSE=as.data.frame(dfSE)
dfP=as.data.frame(dfP)


break()
dfOdds$type=row.names(dfOdds)
dfSE$type=row.names(dfSE)

dataMOdds<-melt(dfOdds)
dataMOddsOld=dataMOdds
dataMSE<-melt(dfSE)
#dataMOdds$value=log(dataMOdds$value)

colnames(dataMOdds)<-c("class","region","odds_ratio")
colnames(dataMSE)<-c("class","region","odds_ratio")

dataMOdds$confidenceInterval<-1.96*dataMSE$odds_ratio
dataMOdds$class<-factor(dataMOdds$class,levels=names(regionsTemp))

pdf(paste0(imageDir,"oddsratios_repeats.pdf"),width=8,height=6)

limits <- aes(ymax = odds_ratio + confidenceInterval, ymin=odds_ratio - confidenceInterval)
p<-ggplot(dataMOdds,aes(class,odds_ratio,colour=region))
p<-p+geom_point(size=2)
p<-p+geom_errorbar(limits,width=0.2)
#p<-p+ylim(-1,1)
p<-p+scale_color_manual(values=c("darkred","darkblue","darkorange","purple"))
p<-p+theme(text = element_text(size=6))
p
dev.off()




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











