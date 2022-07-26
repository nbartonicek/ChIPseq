


#homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/ELF5.macs/")
inPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate/")
outPath_separate=paste0(homedir,"/projects/Chris/project_results/ELF5.macs_separate_clean/")
system(paste0("mkdir -p ",outPath_separate))


macsFiles<-list.files(inPath_separate,pattern="narrow|broad",full.names=T)
results_separate<-GRangesList()
for(macsFile in macsFiles){

  mypeaks <- read.delim(macsFile,header=F,stringsAsFactors=F)
  colnames(mypeaks)[1:9] <- c("chrom", "chromStart", "chromEnd", "name","score", "strand", "fold.enrichment","log10.pvalue", "log10.qvalue")
  mypeaks<-mypeaks[mypeaks$score>=1000,]
  write.table(mypeaks,file=paste0(outPath_separate,basename(macsFile)),quote=F,row.names=F,col.names=F)
}




