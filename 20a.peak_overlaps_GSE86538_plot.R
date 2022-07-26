
library(ggplot2)
library(ChIPpeakAnno)


homedir="/share/ScratchGeneral/nenbar"
inPath=paste0(homedir,"/projects/Chris/project_results/GSE86538/")
projectDir=paste0(homedir,"/projects/Chris")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
cleanRobjectsDir = paste(resultsDir,"/Robjects/cleanPeaks/",sep="")
robjectsDirPeakOverlap = paste(resultsDir,"/Robjects/",sep="")

files<-list.files(robjectsDirPeakOverlap,pattern="peakOverlap_",full.names=T)

results<-list()
for(file in files){
  sampleName=basename(file)
  sampleName=gsub(".Rdata","",sampleName)
  sampleName=gsub("peakOverlap_","",sampleName)
  cat(sampleName)
  cat("\n")
  load(file)
  results[[sampleName]]<-pt$cntOverlaps$pval
}

df<-do.call("rbind",results)
res<-grepl("res",row.names(df))
resDf<-data.frame(resistant=df[res,1],sensitive=df[!res,1])
row.names(resDf)=gsub("resistant_","",row.names(resDf))
resDf

resDfGSE75372=resDf[c(1:3,6:9),]

robjectsDirPeakOverlap = paste(resultsDir,"/Robjects/peakOverlap",sep="")
files<-list.files(robjectsDirPeakOverlap,pattern="peakOverlap_",full.names=T)
files=files[!grepl("GSE32222",files)]

results<-list()
for(file in files){
  sampleName=basename(file)
  sampleName=gsub(".Rdata","",sampleName)
  sampleName=gsub("peakOverlap_","",sampleName)
  cat(sampleName)
  cat("\n")
  load(file)
  results[[sampleName]]<-pt$cntOverlaps$pval
}

df<-do.call("rbind",results)
res<-grepl("res",row.names(df))
resDf<-data.frame(resistant=df[res,1],sensitive=df[!res,1])
row.names(resDf)=gsub("resistant_","",row.names(resDf))
resDf

resDfGSE86548=resDf[c(1:3,6:9),]



robjectsDirPeakOverlap = paste(resultsDir,"/Robjects/peakOverlap",sep="")
files<-list.files(robjectsDirPeakOverlap,pattern="peakOverlap_",full.names=T)
files=files[grepl("GSE32222",files)]

results<-list()
for(file in files){
  sampleName=basename(file)
  sampleName=gsub(".Rdata","",sampleName)
  sampleName=gsub("peakOverlap_","",sampleName)
  cat(sampleName)
  cat("\n")
  load(file)
  results[[sampleName]]<-pt$cntOverlaps$pval
}

df<-do.call("rbind",results)
res<-grepl("res",row.names(df))
resDf<-data.frame(resistant=df[res,1],sensitive=df[!res,1])
row.names(resDf)=gsub("GSE32222_","",row.names(resDf))
resDf
resDfGSE32222=resDf[c(1:3,6:9),]

df=cbind(resDfGSE32222,resDfGSE86548,resDfGSE75372)
colnames(df)=paste0(rep(c("GSE32222_ER_","GSE86548_ER_","GSE75372_FOX_"),each=2),colnames(df))
row.names(df)=gsub("resistant_","",row.names(df))

write.table(df,"../project_results/tables/GSE32222_86548_75372_1.xls",row.names=T,quote=F,sep="\t")


