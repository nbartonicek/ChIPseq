library(rGMAP)
library(data.table)
homedir="../../.."


######## directory structure #######
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/venn/")
annotationDir=paste0(projectDir,"/annotation/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
peakRobjectsDir = paste(resultsDir,"/Robjects/diff/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
rnaseqDir=paste0(resultsDir,"/RNAseq/")
inDir=paste0(resultsDir,"/GSE66733/Hi-C_MCF7_MCF10A_processed_HiCfiles/Heatmaps/chrxchr/40kb/")


cx=c(1:21,"X")
results=list()
for(chr in cx){
  inFiles=list.files(inDir,pattern=paste0("chr",chr,"__"),full.names=T)
  inFile=inFiles[grepl("MCF7-WT",inFiles)]

  df<-read.table(inFile,header=T,row.names=1)
  res = rGMAP(as.matrix(df), resl = 40000,dom_order=10)
  results[[chr]]=as.data.frame(res)
}
  
  pp = plotdom(as.matrix(df), res$hierTads, 0, 300, 40, 40000)
pp
