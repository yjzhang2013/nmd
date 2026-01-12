args <- commandArgs( trailingOnly = TRUE )
fileDir <- args[1]
cell <- args[2]

fileDir <- "/home/u1357/encode/rnaseq/k562/star_mapping_1/"

fs <- list.files(fileDir,pattern = ".SJ.out.tab",full.names = T)
# keep <- paste("chr",c(seq(1:22),"X","Y"),sep="")

for (i in 1:length(fs)){
  cat("\n",i)
  sj <- read.table(file = fs[i], sep="\t")
  sj <- sj[sj$V5>0 & sj$V6==0 & sj$V7>=3,]
  # sj <- sj[sj$V5>0 & sj$V6==0 & sj$V7>=3 & sj$V1%in%keep,]
  if (i==1){
    sjs <- sj
  }else{
    sjs <- rbind(sjs,sj)
  }
}

cat("\nsave")
fileName <- paste("/home/u1357/encode/sjmerge/merge_all_",cell,".SJ.out.tab",sep="")
write.table(sjs,file=fileName,sep="\t",row.names = F,col.names = F, quote = F)

