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
library(R.utils)

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
chrGR<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))
homedir="/share/ScratchGeneral/nenbar"
#homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs/")
inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate/")

inBams=paste0(homedir,"/projects/Chris/project_results/ELF5.picard/")
outPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate_consensus/")
system(paste0("mkdir -p ",outPath))


args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["class"]])){class = args$class} 
if (!is.null(args[["typeID"]])){typeID = args$typeID} 

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
permRobjectsDir=paste(resultsDir,"/Robjects/grammarPermutations/",sep="")

load(paste0(grammarDir,"all_consensus_clean.Rdata"))
load(paste0(grammarDir,"targets.Rdata"))

#load the peaks and differentially expressed genes
load(paste0(cleanRobjectsDir,"all_peaks.Rdata"))

#now for each in enhancer, superenhancer, DE, HOT, resistant, sensitive
#get p-values with ELF5
typeID=as.numeric(typeID)
type=names(consensusReps)[typeID]

DEregions=targets[[class]]
gr<-resize(cleanGRs[["ELF5Dox"]],2000,fix="center")

mat<-findOverlaps(gr,DEregions)
grDE<-gr[unique(queryHits(mat))]

mat<-findOverlaps(consensusReps[[type]],DEregions)
targetDE<-consensusReps[[type]][unique(queryHits(mat))]

pool <- new("permPool", grs=GRangesList(targets[[class]]), N=length(targetDE))
system.time(pt <- peakPermTest(grDE, targetDE, pool=pool, ntimes=1000000))


#pool <- new("permPool", grs=GRangesList(targets[[class]]), N=length(consensusReps[[type]]))
#pt <- peakPermTest(resize(cleanGRs[["ELF5Dox"]],2000,fix="center"), consensusReps[[type]], pool=pool, ntimes=20000)
save(pt,file=paste0(permRobjectsDir,"peakOverlap_",class,"_",type,".Rdata"))






















































