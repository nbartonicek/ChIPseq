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
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(plyr)
library(TFBSTools)
library(JASPAR2018)
library(Biostrings)

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

#TSS distribution
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


#Ok, so we are looking at centered, noncentered
#dox enriched vs depleted
#intergenic vs genic / tamoxifen or not 
#overlap with ELF5 or not

for(class in c("intergenic","genic","resistant","sensitive")){
  for(enriched in c("enriched","depleted")){
    for(centered in c("centered","noncentered")){
      if(class=="intergenic"){


        #first import the enhancers
        enh<-read.table(paste0(annotationDir,"mcf7_enhancers.csv"),sep=",",header=T,stringsAsFactors=F)
        enh=enh[enh$type %in% "Distal_TE",]

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
        grTarget<-liftover(enh.hg19)
      } else if (class=="genic") {
        #TSS distribution
        txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
        grTarget <- promoters(genes(txdb,columns=c("gene_id")), upstream=2500, downstream=2500)
      } else if (class=="resistant") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
      } else if (class=="sensitive") {
        load(paste0(robjectsDir,"tamoxifen.Rdata"))
        grTarget <- tam[[class]]
      } 

      grQuery=cleanGRsPeaks[[paste0("FOXA1_diff_Dox_",enriched)]]
      mat<-findOverlaps(grQuery,grTarget)
      gr<-grQuery[unique(queryHits(mat))]


      gr1kb<-resize(gr,fix="center",width=1000)

      seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
      names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

      if(centered=="centered"){
        opts <- list()
        opts[["species"]] <- 9606
        opts[["name"]] <- "FOXA1"
        #opts[["type"]] <- "SELEX"
        opts[["all_versions"]] <- TRUE
        PFMatrixList <- getMatrixSet(JASPAR2018, opts)
        pwm=PFMatrixList[[1]]
        i=80

        hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )
        lengths<-sapply(hits,length)
        #adjust the GRanges
        gr1kbNew<-gr1kb[lengths==1]
        starts<-sapply(hits,start)
        starts=unlist(starts[lengths==1])
        gr1kbNewShifted<-shift(gr1kbNew,starts-495)

        seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kbNewShifted)
        names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
      }      
      cat(paste0("FoxA1_",centered,"_",enriched,"_",class))
      cat("\n")
      Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/FoxA1_",centered,"_",enriched,"_",class,"_1000.bed"))

    }    
  }
}


