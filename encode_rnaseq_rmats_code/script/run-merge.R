args <- commandArgs( trailingOnly = TRUE )
fileDir <- args[1]
kdsample <- args[2]
outDir <- args[3]

#
patt = paste(kdsample,".*.count",sep="")
fs=list.files(fileDir,pattern=patt)

for (i in c("counts","fpkm","tpm")){
	mt <- do.call(cbind,lapply(fs, function(x){
		read.table(file.path(fileDir,x),header = T,sep = '\t',row.names = 1)[,i]
	}))

	mt<- data.frame(mt)
	row.names(mt) <- row.names(read.table(file.path(fileDir,fs[1]),header = T,sep = '\t',row.names = 1))
	names(mt) <- paste(gsub(".count","",fs),i,sep=".")
	write.table(mt,file=paste(outDir,kdsample,".", i,".merge.txt",sep=""),sep="\t")
}
