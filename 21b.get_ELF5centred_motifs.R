

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
require(ggseqlogo)

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
imageDir=paste0(resultsDir,"/figures/logos/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")

system(paste0("mkdir -p ",cleanRobjectsDir))
system(paste0("mkdir -p ",imageDir))

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




genicFile<-paste0(projectDir,"/project_results/fasta/ELF5_centered_genic_500.bed")
enhFile<-paste0(projectDir,"/project_results/fasta/ELF5_centered_enhancer_500.bed")
nongenicFile<-paste0(projectDir,"/project_results/fasta/ELF5_centered_nongenic_500.bed")

elf5gen<-as.character(readDNAStringSet(genicFile))
elf5enh<-as.character(readDNAStringSet(enhFile))
elf5nongen<-as.character(readDNAStringSet(nongenicFile))

#get substrings


pdf(paste0(imageDir,"ELF5_enhancer.pdf"),width=6,height=3)
ggseqlogo( bUl,method="bits")
dev.off()
