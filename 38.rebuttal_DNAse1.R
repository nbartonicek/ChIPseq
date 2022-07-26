library(rtracklayer)
library(GenomicRanges)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(LOLA)
library(data.table)

homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/GSE72252/")
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
lolaDir=paste0(homedir,"/projects/Chris/project_results/GSE72252/lola/regionDB/hg38")
#import csv files
data1<-read.table(paste0(inPath,"GSE72250_hs_RC_MCF7_DHS_unt_RC_GH1501_3xGH1571_3_FDR0.csv"),sep=",",header=T)
data2<-read.table(paste0(inPath,"GSE72250_hs_RC_MCF7_DHS_E2_RC_GH1503_3xGH1573_3_FDR0.csv"),sep=",",header=T)

gr1<-GRanges(seqnames=data1$chr,IRanges(start=data1$st,end=data1$ed))
gr2<-GRanges(seqnames=data2$chr,IRanges(start=data2$st,end=data2$ed))


liftover<-function(gr){
  gr.hg19=gr
  ch = import.chain(paste0(projectDir,"/annotation/hg19ToHg38.over.chain"))
  seqlevelsStyle(gr.hg19) = "UCSC"  # necessary
  gr.hg38 = liftOver(gr.hg19, ch)
  gr.hg38  = unlist(gr.hg38 )
  genome(gr.hg38) = "hg38"
  return(gr.hg38)
}

gr1hg38<-liftover(gr1)
gr2hg38<-liftover(gr2)

export(gr1hg38,paste0(inPath,"unt.bed"))
export(gr2hg38,paste0(inPath,"E2.bed"))

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

#get the FOXA1 enriched
gr<-cleanGRsPeaks[[16]]


sum(countOverlaps(gr,gr1hg38)>0)
sum(countOverlaps(gr,gr2hg38)>0)



gr<-cleanGRsPeaks[[4]]
sum(countOverlaps(gr,gr1hg38)>0)
sum(countOverlaps(gr,gr2hg38)>0)


gr<-cleanGRsPeaks[[5]]
sum(countOverlaps(gr,gr1hg38)>0)
sum(countOverlaps(gr,gr2hg38)>0)


gr<-cleanGRsPeaks[[17]]
sum(countOverlaps(gr,gr1hg38)>0)
sum(countOverlaps(gr,gr2hg38)>0)






#get the universe
files<-list.files(inPath,pattern="bedgraph",full.names=T)
results<-GRangesList()
for(file in files){
  cat(file)
  gr<-fread(file)
  gr<-as.data.frame(gr)
  gr<-GRanges(seqnames=gr[,1],IRanges(start=gr[,2],end=gr[,3]))
  gr<-reduce(gr)
  results[[basename(file)]]<-liftover(gr)
}

universe<-reduce(c(unlist(results),gr1hg38,gr2hg38))
universe<-cleanGRsPeaks[[1]]



dbPath = system.file(lolaDir, "hg38", package="LOLA")
regionDB = loadRegionDB(lolaDir)

load(paste0(cleanRobjectsDir,"all_peaks_regions.Rdata"))

userSets<-cleanGRsPeaks[[1]]

#only take those ELF5 regions that are fully in the universe
mat<-findOverlaps(userSets,universe,type="within")
userSets<-userSets[unique(queryHits(mat))]

locResults = runLOLA(userSets, universe, regionDB, cores=1)

foxDox<-cleanGRsPeaks[[4]]
foxNoDox<-cleanGRsPeaks[[5]]


#ELF5FoxDox
mat<-findOverlaps(userSets,foxDox)
elf5FoxDox<-userSets[unique(queryHits(mat))]
mat<-findOverlaps(userSets,foxNoDox)
elf5FoxNoDox<-userSets[unique(queryHits(mat))]

sum(countOverlaps(elf5FoxDox,gr1hg38)>0)
sum(countOverlaps(elf5FoxDox,gr2hg38)>0)


sum(countOverlaps(elf5FoxNoDox,gr1hg38)>0)
sum(countOverlaps(elf5FoxNoDox,gr2hg38)>0)

mat<-findOverlaps(userSets,foxDox)
elf5NoFoxDox<-userSets[-unique(queryHits(mat))]
mat<-findOverlaps(userSets,foxNoDox)
elf5NoFoxNoDox<-userSets[-unique(queryHits(mat))]

sum(countOverlaps(elf5NoFoxDox,gr1hg38)>0)
sum(countOverlaps(elf5NoFoxDox,gr2hg38)>0)

sum(countOverlaps(elf5NoFoxNoDox,gr1hg38)>0)
sum(countOverlaps(elf5NoFoxNoDox,gr2hg38)>0)


######## overlap of DNAse1 sites
mat<-findOverlaps(gr1hg38,foxDox)
dnaseUtFoxDox<-gr1hg38[unique(queryHits(mat))]
mat<-findOverlaps(gr2hg38,foxNoDox)
dnaseE2FoxDox<-gr2hg38[unique(queryHits(mat))]

sum(countOverlaps(userSets,dnaseUtFoxDox)>0)
sum(countOverlaps(userSets,dnaseE2FoxDox)>0)

######## overlap of DNAse1 sites
mat<-findOverlaps(userSets,gr1hg38)
elf5NoDNAseUt<-userSets[-unique(queryHits(mat))]
mat<-findOverlaps(userSets,gr2hg38)
elf5NoDNAseE2<-userSets[-unique(queryHits(mat))]

sum(countOverlaps(elf5NoDNAseUt,foxNoDox)>0)
sum(countOverlaps(elf5NoDNAseE2,foxNoDox)>0)




data(HOT.spots)
data(wgEncodeTfbsV3)
hotGR <- liftover(reduce(unlist(HOT.spots)))
wgEncodeTfbsV3<-liftover(wgEncodeTfbsV3)   
removeOl <- function(.ele){
        ol <- findOverlaps(.ele, hotGR)
        if(length(ol)>0) .ele <- .ele[-unique(queryHits(ol))]
        .ele
}
userSetsClean<-removeOl(userSets)
values(userSetsClean)<-NULL


poolL<-list()
poolL[["sensitive"]] <- new("permPool", grs=GRangesList(universe), N=length(userSets))


DNAsepeaks<-GRangesList()
DNAsepeaks[["ut"]]<-gr1hg38
DNAsepeaks[["E2"]]<-gr2hg38
pt <- peakPermTest(userSets, DNAsepeaks[[1]], pool=poolL[[1]], ntimes=1000)





ERpeaks<-cleanGRsPeaks[c(1,2,3,6,7,8,9,10,11,12,13)]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

sensitive<-enrichPeakOverlap(queryPeak     = sensitiveUniq,
                  targetPeak    = ERpeaks,
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 5000,
                  chainFile     = NULL,
                  verbose       = FALSE)


resistant<-enrichPeakOverlap(queryPeak     = resistantUniq,
                  targetPeak    = ERpeaks,
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 5000,
                  chainFile     = NULL,
                  verbose       = FALSE)










