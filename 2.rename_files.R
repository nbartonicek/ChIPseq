
homeDir="/share/ScratchGeneral/nenbar/projects/Chris"
sampleName="L5"
inDir=paste0(homeDir,"/raw_files/",sampleName,"/merged")
outDir=paste0(homeDir,"/raw_files/",sampleName,"/renamed")
system(paste0("mkdir ",outDir))

system(paste0("ls ",inDir,"/*.gz >filenames.txt"))
filenames<-read.table("filenames.txt",header=F)
conversions<-read.table("id_conversion.txt",header=F,colClasses="character")

for(i in 1:length(conversions[,1])){

        presentFiles<-as.character(filenames[grepl(as.character(paste0("-",conversions[i,1])),filenames[,1]),])
        cat(presentFiles)
        if(length(presentFiles)>0){
                for(presentFile in presentFiles){
                        if(grepl("_R1",presentFile)){
				system(paste0("mv ",presentFile," ",outDir,"/",conversions[i,2],"_R1.fastq.gz"))
			             #cat(paste0("dx mv ",presentFile," ",conversions[i,2],"_R1.fastq.gz"))
                        } else {
				system(paste0("mv ",presentFile," ",outDir,"/",conversions[i,2],"_R2.fastq.gz"))
			             #cat(paste0("dx mv ",presentFile," ",conversions[i,2],"_R2.fastq.gz"))
                        }
		
			#newName<-gsub(conversions[i,1],conversions[i,2],presentFile)
			#newName<-paste0(conversions[i,2],".",presentFile)
                        #system(paste("dx mv",presentFile,newName))
                        #cat(presentFile)
                        cat("\n")

                }
        }
}














