#R --vanilla <bam2bw.R --inFile='LTP_D1_2h_Rat1_LDG' 

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("preprocessCore")
#biocLite("BSgenome.Hsapiens.NCBI.GRCh38")
#library("preprocessCore")
library(R.utils)
#library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicAlignments)

chrFile<-"/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg38_ercc/chrNameLength.txt"
chrDf<-read.table(chrFile,header=F)
chrDf=chrDf[!(grepl("_",chrDf[,1])),]
chrDf=chrDf[!(grepl("-",chrDf[,1])),]
chrDf=chrDf[order(chrDf[,1]),]

chrs<-chrDf[,2]
names(chrs)<-chrDf[,1]
#w <- as.integer(commandArgs(TRUE)[2])
#bamFile <- commandArgs(TRUE)[3]
#bigwigFile <- commandArgs(TRUE)[4]
####
args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["inFile"]])){inFile = args$inFile}


SampleName<-basename(inFile)
cat(SampleName,'\n')
Sample<-gsub('.bam','',SampleName)
outDir="/share/ScratchGeneral/nenbar/projects/Chris/project_results/ELF5.bigwig"
dir.create(outDir)


RNAbamTobw <- function(file, name, stranded=TRUE, firstStrand=TRUE, paired=TRUE) {
  require(rtracklayer)
   if (paired) { cat('Data is paired','\n')
       system.time(rs<-readGAlignmentPairs(file))
		#levels(seqnames(rs))[levels(seqnames(rs))=="MT"]<-'M'
	#	rs<-renameSeqlevels(rs, c(MT="M"))
	#	seqlevels(rs)<-paste0('chr',seqlevels(rs))
		#seqlengths(rs)<-chrs
        rs<-rs[seqnames(rs) %in% names(chrs)]
       rs.count <- length(rs)/1e6
       rs.cov <- unlist(grglist(rs))
       if (stranded) { cat('Data is stranded','\n')
           rs.cov <- mclapply(split(rs.cov, strand(rs.cov))[c("+", "-")], coverage,mc.cores=5)
           #if (!firstStrand) names(rs.cov) <- rev(names(rs.cov))
           export(rs.cov[["-"]]/rs.count, BigWigFile(paste(outDir,"/",name, "_plus.bw", sep="")))
           export(rs.cov[["+"]]/rs.count, BigWigFile(paste(outDir,"/",name, "_minus.bw", sep="")))
           export(rs.cov[["+"]]/-rs.count, BigWigFile(paste(outDir,"/",name, "_minusminus.bw", sep="")))
       } else export(coverage(rs.cov)/rs.count, BigWigFile(paste(outDir,"/",name, ".bw", sep="")))
	}
}
RNAbamTobw(file=inFile,Sample,stranded=FALSE,paired=TRUE,firstStrand=TRUE)
