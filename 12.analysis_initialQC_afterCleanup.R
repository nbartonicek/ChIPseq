

library(ChIPQC)
library(GenomicRanges)
library(VennDiagram)
library(ChIPpeakAnno)
library(ggplot2)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(reshape2)
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
  bamControl=bamControls
)

#liftover of blacklist sections
data(blacklist_hg19)
path = system.file(package="rtracklayer", "extdata", "../annotation/hg38ToHg19.over.chain")
ch = import.chain( "../annotation/hg19ToHg38.over.chain")
seqlevelsStyle(blacklist.hg19) = "UCSC"  # necessary
blacklist = liftOver(blacklist.hg19, ch)
blacklist = unlist(blacklist)
genome(blacklist) = "hg38"


experiment = ChIPQC(samples,annotation="hg38",chromosomes="chr14",blacklist=blacklist)
#ChIPQCreport(experiment)

pdf(paste0(imageDir,"CoverageHist.pdf"))
plotCoverageHist(experiment,facetBy=c("Factor","Condition"))
dev.off()
pdf(paste0(imageDir,"CC.pdf"))
plotCC(experiment,facetBy=c("Factor","Condition"))
dev.off()

pdf(paste0(imageDir,"Regi.pdf"),width=8,height=8)
plotRegi(experiment,facetBy=c("Factor"))+scale_fill_gradient2(low="white",high="darkblue")+ theme(text = element_text(size=10))
dev.off()


pdf(paste0(imageDir,"Regi.pdf"),width=8,height=8)
plotRegi(experiment)+scale_fill_gradient2(low="white",high="darkblue")+ theme(text = element_text(size=10)axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


df<-regi(experiment)
df=df[,-grep("^Dox|^NoDox",colnames(df))]
colnames(df)[is.na(colnames(df))]=paste0("input",1:8)
dataM<-melt(df)
colnames(dataM)<-c("region","sample","logOdds")
dataM$type=gsub("\\d$","",dataM$sample)
dataM$class=gsub("NoDox|Dox","",dataM$type)

pdf(paste0(imageDir,"Regi_custom.pdf"),width=8,height=8)
p<-ggplot(dataM,aes(region,sample))
p<-p+geom_tile(aes(fill=logOdds))
p<-p+scale_fill_gradient2(low="white",high="darkblue")
p<-p+facet_grid(type~.,drop = TRUE,scales = "free_y")
p<-p+theme(text = element_text(size=9),axis.text.x = element_text(angle = 90, hjust = 1))
p
dev.off()

pdf(paste0(imageDir,"peakProfile.pdf"))
plotPeakProfile(experiment,facetBy=c("Factor","Condition"))
dev.off()
pdf(paste0(imageDir,"FractionInPeaks"))
plotFrip(experiment,facetBy=c("Factor","Condition"))
dev.off()
pdf(paste0(imageDir,"heatmap.pdf"))
plotCorHeatmap(experiment,attributes=c("Factor","Condition","Replicate"))
dev.off()
pdf(paste0(imageDir,"PCA.pdf"))
plotPrincomp(experiment,attributes=c("Factor","Condition"))
dev.off()

