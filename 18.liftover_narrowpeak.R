library(rtracklayer)
library(GenomicRanges)

homedir="../../.."
inPath=paste0(homedir,"/projects/Chris/project_results/beds/")
outPath=paste0(homedir,"/projects/Chris/project_results/beds/hg19/")
projectDir=paste0(homedir,"/projects/Chris")
system(paste0("mkdir ",outPath))

files<-list.files(inPath,full.names=T,pattern="Peak")

for(file in files){
	cat(basename(file))
	cat("\n")
	df<-read.table(file)
	gr<-GRanges(seqnames=df$V1,IRanges(start=df$V2,end=df$V3))
	values(gr)<-df[,4:dim(df)[2]]
	gr.hg38=gr
  	ch = import.chain(paste0(projectDir,"/annotation/hg38ToHg19.over.chain"))
  	seqlevelsStyle(gr.hg38) = "UCSC"  # necessary
  	gr.hg19 = liftOver(gr.hg38, ch)
  	gr.hg19  = unlist(gr.hg19 )
  	genome(gr.hg19) = "hg19"
  	df19<-as.data.frame(gr.hg19)
  	df19<-df19[,c(1,2,3,6:(dim(df)[2]+2))]
  	write.table(df19,paste0(outPath,basename(file)),quote=F,row.names=F,col.names=F,sep="\t")
}


#for(geneType in names(deGenesL)){
#    cat(geneType)
#    cat("\n")
#	deGenes<-deGenesL[[geneType]]
#	deGenes<-deGenes[!is.na(deGenes)]
#	gns<-bitr(deGenes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#	promoterSubset=promoter[promoter$gene_id %in% gns$ENTREZID]
#	merged<-merge(gns,promoterSubset$gene_id,by.x="ENTREZID",by.y=1)
#	promoterSubset$gene_id=merged$SYMBOL
#	names(promoterSubset)<-merged$SYMBOL
#	promoterSubset<-resize(promoterSubset,width=1000,fix="center")
#    export(promoterSubset,paste0(geneType,".bed"))  
#
#    gr.hg38=promoterSubset
#  	ch = import.chain(paste0(projectDir,"/annotation/hg38ToHg19.over.chain"))
#  	seqlevelsStyle(gr.hg38) = "UCSC"  # necessary
#  	gr.hg19 = liftOver(gr.hg38, ch)
#  	gr.hg19  = unlist(gr.hg19 )
#  	genome(gr.hg19) = "hg19"
#	export(gr.hg19,paste0(geneType,"_hg19.bed"))  
#
#}