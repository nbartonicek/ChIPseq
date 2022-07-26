

library(ChIPQC)
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
library(ggplot2)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs/")
inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate_consensus/")

inBamsDir=paste0(homedir,"/projects/Chris/project_results/ELF5.picard/")
inBams<-list.files(inBamsDir,pattern="bam$",full.names=T)
controls<-inBams[grepl("Input",inBams)]
controlsNoDox<-controls[grepl("NoDox",controls)]
controlsDox<-controls[!grepl("NoDox",controls)]

inBams<-inBams[!grepl("Input",inBams)]

######## directory structure #######
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/macs_clean_consensus/")
system(paste0("mkdir -p ",imageDir))
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")

#create sample sheet
separateFiles<-list.files(inPath_separate,full.names=T,pattern="narrow|broad")

sampleNames<-gsub("_peaks.*","",basename(separateFiles))
treatment<-rep("Dox",length(sampleNames))
treatment[grepl("NoDox",sampleNames)]="NoDox"
uniqueIDs=gsub("\\d$","",sampleNames)

seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to"))
replicate=unlist(c(seq2(from = c(1), to = table(uniqueIDs), by = 1)))

bamControls=c(controlsDox[c(1,2,3,4,5,6,1,2)],controlsNoDox[1:4],controlsDox[c(3,4,1,2)],controlsNoDox[1:4],controlsDox[3:4],controlsNoDox[1:2])
samples<-data.frame(SampleID=sampleNames,
  Tissue=treatment,
  Condition=treatment,
  Factor=uniqueIDs,
  Replicate=replicate,
  bamReads=inBams,
  Peaks=separateFiles,
  ControlID=treatment,
  bamControl=bamControls,
  PeakCaller="narrow"
)

#liftover of blacklist sections
data(blacklist_hg19)
path = system.file(package="rtracklayer", "extdata", "../annotation/hg38ToHg19.over.chain")
ch = import.chain( "../annotation/hg19ToHg38.over.chain")
seqlevelsStyle(blacklist.hg19) = "UCSC"  # necessary
blacklist = liftOver(blacklist.hg19, ch)
blacklist = unlist(blacklist)
genome(blacklist) = "hg38"

samples<-samples[grepl("FOXA1",samples$SampleID),]

experiment = dba(sampleSheet=samples)
experiment <- dba.count(experiment)
experiment <- dba.contrast(experiment)
experiment <- dba.analyze(experiment)
experiment.DB <- dba.report(experiment)


pdf(paste0(imageDir,"diffBind_ER_heatmap.pdf"))
plot(experiment)
dev.off()


pdf(paste0(imageDir,"diffBind_ER_MA.pdf"))
dba.plotMA(experiment)
dev.off()


pdf(paste0(imageDir,"diffBind_ER_volcano.pdf"))
dba.plotVolcano(experiment)
dev.off()


pdf(paste0(imageDir,"diffBind_ER_boxplot.pdf"))
dba.plotBox(experiment)
dev.off()

#venn of diff
save(experiment.DB,file=paste0(robjectsDir,sampleName,"_diff.Rdata"))

#experiment <- dba.peakset(experiment, consensus=DBA_CONDITION, minOverlap=0.33)
#pdf(paste0(imageDir,"diffBind_ER_Venn.pdf"))
#dba.plotVenn(experiment,experiment$masks$consensus)
#dev.off()

#ChIPQCreport(experiment)
