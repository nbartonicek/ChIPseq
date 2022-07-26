

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
library(TFBSTools)
library(JASPAR2018)
library(Biostrings)

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

#1. clean up from blacklist (186 in merged peaks, 50 in macs, 156 in Carroll data)
#2. plot bar plot for distribution and percentage overlap
#3. plot intensity/p-value/fdr vs presence in consensus/macs2/macs2 strict
#4. what is macs2 strict?
#liftover of blacklist sections
#data(blacklist_hg19)
#path = system.file(package="rtracklayer", "extdata", "../annotation/hg38ToHg19.over.chain")
#ch = import.chain( "../annotation/hg19ToHg38.over.chain")
#seqlevelsStyle(blacklist.hg19) = "UCSC"  # necessary
#blacklist = liftOver(blacklist.hg19, ch)
#blacklist = unlist(blacklist)
#genome(blacklist) = "hg38"

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


#first Diff
gr500<-resize(threeDiff,fix="center",width=500)
gr1kb<-resize(threeDiff,fix="center",width=1000)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCEDIFF_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5_FoxA1_ER_diff500bp.bed"))


seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
names(seq) = paste0("SEQUENCEDIFF_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5_FoxA1_ER_diff1kb.bed"))



#first Diff
gr500<-resize(threeAll,fix="center",width=500)
gr1kb<-resize(threeAll,fix="center",width=1000)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5_FoxA1_ER_all500bp.bed"))


seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5_FoxA1_ER_all1kb.bed"))

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "ELF5"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
ELF5pwm=PFMatrixList[[1]]

#find ELF5 hits
for(i in c(45,50,60,70,75,80,85,90,95,99)){


  hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0(i,"%")) )
  lengths<-sapply(hitsELF5,length)
  cat(i)
  cat("\t")
  cat(table(lengths))
  cat("\n")
}


seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("75%")) )
lengths<-sapply(hitsELF5,length)
cat(table(lengths))
cat("\n")
#adjust the GRanges

gr500New<-gr500[lengths==1]
starts<-sapply(hitsELF5,start)
starts=unlist(starts[lengths==1])
gr500NewShifted<-shift(gr500New,starts-245)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500NewShifted)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("75%")) )
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5motif_centred_500_75.bed"))

#do ELF5 and FOXA1 distance

mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],cleanGRsPeaks[["FOXA1Dox"]])
bothAll<-cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat))]
gr500<-resize(bothAll,fix="center",width=500)
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "ELF5"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
ELF5pwm=PFMatrixList[[1]]

hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("80%")) )
lengths<-sapply(hitsELF5,length)
cat(table(lengths))
cat("\n")
#adjust the GRanges

gr500New<-gr500[lengths==1]
starts<-sapply(hitsELF5,start)
starts=unlist(starts[lengths==1])
gr500NewShifted<-shift(gr500New,starts-245)



seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500NewShifted)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("80%")) )
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5motif_centred_FOXA1_noER_500_80.bed"))


######do ELF5 and ER distance

mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],cleanGRsPeaks[["ERDox"]])
bothAll<-cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat))]
gr500<-resize(bothAll,fix="center",width=500)
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "ELF5"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
ELF5pwm=PFMatrixList[[1]]


#find ELF5 hits
#for(i in c(60,70,75,80,85,90,95,99)){
#  hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0(i,"%")) )
#  lengths<-sapply(hitsELF5,length)
#  cat(i)
#  cat("\t")
#  cat(table(lengths))
#  cat("\n")
#}

hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("80%")) )
lengths<-sapply(hitsELF5,length)
cat(table(lengths))
cat("\n")
#adjust the GRanges

gr500New<-gr500[lengths==1]
starts<-sapply(hitsELF5,start)
starts=unlist(starts[lengths==1])
gr500NewShifted<-shift(gr500New,starts-245)



seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500NewShifted)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("80%")) )
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5motif_centered_ER_noFOXA1_500_80.bed"))



######do ELF5 and H3K4me3 distance

mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],cleanGRsPeaks[["H3K4me3Dox"]])
bothAll<-cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat))]
gr500<-resize(bothAll,fix="center",width=500)
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "ELF5"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
ELF5pwm=PFMatrixList[[1]]


#find ELF5 hits
#for(i in c(70,75,80,85,90,95,99)){
#  hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0(i,"%")) )
#  lengths<-sapply(hitsELF5,length)
#  cat(i)
#  cat("\t")
#  cat(table(lengths))
#  cat("\n")
#}

hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("85%")) )
lengths<-sapply(hitsELF5,length)
cat(table(lengths))
cat("\n")
#adjust the GRanges

gr500New<-gr500[lengths==1]
starts<-sapply(hitsELF5,start)
starts=unlist(starts[lengths==1])
gr500NewShifted<-shift(gr500New,starts-245)



seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500NewShifted)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0("85%")) )
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5motif_centered_H3K4me3_only_500_80.bed"))



######do FOXA1 and ER distance

mat<-findOverlaps(cleanGRsPeaks[["FOXA1Dox"]],cleanGRsPeaks[["ERDox"]])
bothAll<-cleanGRsPeaks[["FOXA1Dox"]][unique(queryHits(mat))]
gr500<-resize(bothAll,fix="center",width=500)
seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "FOXA1"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
pwm=PFMatrixList[[1]]


#find ELF5 hits
#for(i in c(70,75,80,85,90,95,99)){
#  hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )
#  lengths<-sapply(hits,length)
#  cat(i)
#  cat("\t")
#  cat(table(lengths))
#  cat("\n")
#}

i=85
hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )
lengths<-sapply(hits,length)
cat(table(lengths))
cat("\n")
#adjust the GRanges

gr500New<-gr500[lengths==1]
starts<-sapply(hits,start)
starts=unlist(starts[lengths==1])
gr500NewShifted<-shift(gr500New,starts-245)



seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr500NewShifted)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
hits <- lapply(seq, function(x) matchPWM(as.matrix(pwm), x, min.score=paste0(i,"%")) )
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/FOXA1motif_centered_ER_500_",i,".bed"))


##### find repeats?

mat<-findOverlaps(cleanGRsPeaks[["ELF5Dox"]],cleanGRsPeaks[["FOXA1Dox"]])
both<-cleanGRsPeaks[["ELF5Dox"]][unique(queryHits(mat))]
mat<-findOverlaps(both,cleanGRsPeaks[["ERDox"]])
threeAll<-both[unique(queryHits(mat))]

#first Diff
gr1kb<-resize(threeAll,fix="center",width=1000)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/ELF5_FoxA1_ER_all1kb.bed"))

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "ELF5"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
ELF5pwm=PFMatrixList[[1]]

#find ELF5 hits
#for(i in c(75,80,85,90,95,99)){
#
#
#  hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0(i,"%")) )
#  lengths<-sapply(hitsELF5,length)
#  cat(i)
#  cat("\t")
#  cat(table(lengths))
#  cat("\n")
#}

i=80

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kb)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))

hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0(i,"%")) )
lengths<-sapply(hitsELF5,length)
cat(table(lengths))
cat("\n")
#adjust the GRanges

gr1kbNew<-gr1kb[lengths==1]
starts<-sapply(hitsELF5,start)
starts=unlist(starts[lengths==1])
gr1kbNewShifted<-shift(gr1kbNew,starts-495)

seq = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, gr1kbNewShifted)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
hitsELF5 <- lapply(seq, function(x) matchPWM(as.matrix(ELF5pwm), x, min.score=paste0(i,"%")) )
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/repeatHunt_ELF5motif_centred_all_500_",i,".bed"))



