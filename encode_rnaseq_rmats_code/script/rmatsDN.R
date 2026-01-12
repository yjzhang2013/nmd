#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("run rmatsDN")

# Add command line arguments
#p <- add_argument(p, "--gtfExon", help="output prefix", type="character")

p <- add_argument(p, "--fileDir", help="input: rMATS ASEs file", type="character")
p <- add_argument(p, "--nRep", help="number of repeat")
p <- add_argument(p, "--outDir", help="output prefix", type="character")
p <- add_argument(p, "--cutoff", help="cutoff of ASE associated juction reads (SJC+IJC) ")
p <- add_argument(p, "--level", help="sample level for ramts",type="character")
p <- add_argument(p, "--novel", help="dose rmats results with novel events",type="logical")

# Parse the command line arguments
argv <- parse_args(p)

# 需要进行运算的整数型输入宜 as.numeric 转化
#gtfExon <- argv$gtfExon

fileDir <- argv$fileDir
n <- as.numeric(argv$nRep)
outPref <- argv$outDir
cutoff <- as.numeric(argv$cutoff)
level <- argv$level
novel <- argv$novel

# 
# fileDir="/home/u1357/encode/rnaseq/hepg2/rmats/ENCSR999MAE/"
# outPref="/home/u1357/encode/rnaseq/hepg2/rmatsDn/ENCSR999MAE/"
# dir.create(outPref,recursive = T)
# gtfExon="/home/u1357/encode/refv19/gencode.v19.processed.exon.txt"
# nRep=n=2
# cutoff=10
# level="nc,kd"
# novel=T


## def function
rmats2region <- function(infoAS,spliceType){
  infoAS$GeneID <- gsub('"','',infoAS$GeneID)
  infoAS$geneSymbol <- gsub('"','',infoAS$geneSymbol)
  infoAS$spliceType <- spliceType
  ## for exclusive event
  if(spliceType=="A3SS"){
    ## sequential locations of event exons
    infoAS$eventRegion <- ifelse(infoAS$strand=="+",
                                 paste(paste(infoAS$flankingES,infoAS$flankingEE,sep=":"),paste(infoAS$longExonStart_0base,infoAS$shortES,sep=":"),paste(infoAS$shortES,infoAS$shortEE,sep=":"),sep="-"),
                                 paste(paste(infoAS$shortES,infoAS$shortEE,sep=":"),paste(infoAS$shortEE,infoAS$longExonEnd,sep=":"),paste(infoAS$flankingES,infoAS$flankingEE,sep=":"),sep="-"))
    ## the longer alt exon in gtf annotation
    infoAS$altExon <- ifelse(infoAS$strand=="+",
                             paste(infoAS$longExonStart_0base,infoAS$shortES,sep=":"),
                             paste(infoAS$shortEE,infoAS$longExonEnd,sep=":"))
  }else if(spliceType=="A5SS"){
    infoAS$eventRegion <- ifelse(infoAS$strand=="+",
                                 paste(paste(infoAS$shortES,infoAS$shortEE,sep=":"),paste(infoAS$shortEE,infoAS$longExonEnd,sep=":"),paste(infoAS$flankingES,infoAS$flankingEE,sep=":"),sep="-"),
                                 paste(paste(infoAS$flankingES,infoAS$flankingEE,sep=":"),paste(infoAS$longExonStart_0base,infoAS$shortES,sep=":"),paste(infoAS$shortES,infoAS$shortEE,sep=":"),sep="-"))
    
    infoAS$altExon <- ifelse(infoAS$strand=="+",
                             paste(infoAS$shortEE,infoAS$longExonEnd,sep=":"),
                             paste(infoAS$longExonStart_0base,infoAS$shortES,sep=":"))
    
  }else if(spliceType=="RI"){
    infoAS$eventRegion <- paste(paste(infoAS$upstreamES,infoAS$upstreamEE,sep=":"),
                                paste(infoAS$upstreamEE,infoAS$downstreamES,sep=":"),
                                paste(infoAS$downstreamES,infoAS$downstreamEE,sep=":"),sep="-")
    infoAS$altExon <- paste(infoAS$upstreamEE,infoAS$downstreamES,sep=":")
    
  }else if(spliceType=="SE"){
    infoAS$eventRegion <- paste(paste(infoAS$upstreamES,infoAS$upstreamEE,sep=":"),
                                paste(infoAS$exonStart_0base,infoAS$exonEnd,sep=":"),
                                paste(infoAS$downstreamES,infoAS$downstreamEE,sep=":"),sep="-")
    infoAS$altExon <- paste(infoAS$exonStart_0base,infoAS$exonEnd,sep=":")
    
  }else if(spliceType=="MXE"){
    infoAS$eventRegion <- paste(paste(infoAS$upstreamES,infoAS$upstreamEE,sep=":"),
                                paste(infoAS$X1stExonStart_0base,infoAS$X1stExonEnd,sep=":"),
                                paste(infoAS$X2ndExonStart_0base,infoAS$X2ndExonEnd,sep=":"),
                                paste(infoAS$downstreamES,infoAS$downstreamEE,sep=":"),sep="-")
    
    infoAS$altExon <- paste(paste(infoAS$X1stExonStart_0base,infoAS$X1stExonEnd,sep=":"),
                            paste(infoAS$X2ndExonStart_0base,infoAS$X2ndExonEnd,sep=":"),sep="-")
  }
  return(infoAS)
}


rmats2bed <- function(infoAS,spliceType){
  ## for bed file
  if (spliceType=="MXE"){
    #
    infoAS$eventStart <- sapply(infoAS$eventRegion, function(x){y <- strsplit(x,':')[[1]][1]})
    infoAS$eventEnd <- sapply(infoAS$eventRegion, function(x){y <- strsplit(x,':')[[1]][5]})
    #
    infoAS$altStart <- sapply(infoAS$altExon, function(x){y <- strsplit(x,':')[[1]][1]})
    infoAS$altEnd <- sapply(infoAS$altExon, function(x){y <- strsplit(x,':')[[1]][3]})
    #
    infoAS$exonNum <- 4
    #
    infoAS$exonLength <- paste(as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][1],':')[[1]][2]}))-
                                 as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][1],':')[[1]][1]})),
                               as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][2],':')[[1]][2]}))-
                                 as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][2],':')[[1]][1]})),
                               as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][3],':')[[1]][2]}))-
                                 as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][3],':')[[1]][1]})),
                               as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][4],':')[[1]][2]}))-
                                 as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][4],':')[[1]][1]})),
                               "",sep=",")
    #
    infoAS$exonStart <- paste(0,
                              as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][2],':')[[1]][1]}))-
                                as.numeric(infoAS$eventStart),
                              as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][3],':')[[1]][1]}))-
                                as.numeric(infoAS$eventStart),
                              as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][4],':')[[1]][1]}))-
                                as.numeric(infoAS$eventStart),
                              "",sep=",")
  }else{
    #
    infoAS$eventStart <- sapply(infoAS$eventRegion, function(x){y <- strsplit(x,':')[[1]][1]})
    infoAS$eventEnd <- sapply(infoAS$eventRegion, function(x){y <- strsplit(x,':')[[1]][4]})
    #
    infoAS$altStart <- sapply(infoAS$altExon, function(x){y <- strsplit(x,':')[[1]][1]})
    infoAS$altEnd <- sapply(infoAS$altExon, function(x){y <- strsplit(x,':')[[1]][2]})
    #
    infoAS$exonNum <- 3
    #
    infoAS$exonLength <- paste(as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][1],':')[[1]][2]}))-
                                 as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][1],':')[[1]][1]})),
                               as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][2],':')[[1]][2]}))-
                                 as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][2],':')[[1]][1]})),
                               as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][3],':')[[1]][2]}))-
                                 as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][3],':')[[1]][1]})),
                               "",sep=",")
    #
    infoAS$exonStart <- paste(0,
                              as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][2],':')[[1]][1]}))-
                                as.numeric(infoAS$eventStart),
                              as.numeric(sapply(infoAS$eventRegion, function(x){y <- strsplit(strsplit(x,'-')[[1]][3],':')[[1]][1]}))-
                                as.numeric(infoAS$eventStart),
                              "",sep=",")
  }
  return(infoAS)
}


## novel AS
if (novel==T){
  a3ssnovel <- read.table(paste(fileDir,"fromGTF.novelSpliceSite.A3SS.txt",sep=""),sep="\t", header=T, quote="")
  a5ssnovel <- read.table(paste(fileDir,"fromGTF.novelSpliceSite.A5SS.txt",sep=""),sep="\t", header=T,quote="")
  mxenovel <- read.table(paste(fileDir,"fromGTF.novelSpliceSite.MXE.txt",sep=""),sep="\t", header=T,quote="")
  senovel <- read.table(paste(fileDir,"fromGTF.novelSpliceSite.SE.txt",sep=""),sep="\t", header=T,quote="")
  rinovel <- read.table(paste(fileDir,"fromGTF.novelSpliceSite.RI.txt",sep=""),sep="\t", header=T,quote="")
  
  a3ssnovel <- rmats2region(a3ssnovel,"A3SS")
  a5ssnovel <- rmats2region(a5ssnovel,"A5SS")
  mxenovel <- rmats2region(mxenovel,"MXE")
  senovel <- rmats2region(senovel,"SE")
  rinovel <- rmats2region(rinovel,"RI")
  
  colname = c("GeneID","ID","geneSymbol","chr","strand","spliceType","eventRegion","altExon")
  aseNovel <- rbind(a3ssnovel[,colname], a5ssnovel[,colname], mxenovel[,colname], senovel[,colname], rinovel[,colname])
  aseNovel$idType <- paste(aseNovel$geneSymbol, aseNovel$spliceType,aseNovel$ID,sep="_")
}


# input data
a3ss <- read.table(paste(fileDir,"A3SS.MATS.JC.txt",sep=""),sep="\t", header=T, quote="")
a5ss <- read.table(paste(fileDir,"A5SS.MATS.JC.txt",sep=""),sep="\t", header=T,quote="")
mxe <- read.table(paste(fileDir,"MXE.MATS.JC.txt",sep=""),sep="\t", header=T,quote="")
se <- read.table(paste(fileDir,"SE.MATS.JC.txt",sep=""),sep="\t", header=T,quote="")
ri <- read.table(paste(fileDir,"RI.MATS.JC.txt",sep=""),sep="\t", header=T,quote="")

a3ss <- rmats2region(a3ss,"A3SS")
a5ss <- rmats2region(a5ss,"A5SS")
mxe <- rmats2region(mxe,"MXE")
se <- rmats2region(se,"SE")
ri <- rmats2region(ri,"RI")

col = c("GeneID","ID","geneSymbol","chr","strand","spliceType","eventRegion","altExon","IJC_SAMPLE_1","SJC_SAMPLE_1","IJC_SAMPLE_2","SJC_SAMPLE_2","IncLevel1","IncLevel2","IncLevelDifference","PValue","FDR")
ase <- rbind(a3ss[,col], a5ss[,col], mxe[,col], se[,col], ri[,col])

ase$idType <- paste(ase$geneSymbol, ase$spliceType,ase$ID,sep="_")
if (novel==T){
  ase$novel <- ifelse(ase$idType%in%aseNovel$idType,"novel","known")
}
write.table(ase,file=paste(outPref,"all_ase_nofiltrated.txt",sep=""),sep="\t",row.names = FALSE,quote = FALSE)

cat("All ASEs Detected:\n")
table(ase$spliceType)
cat("Sum:")
cat(sum(table(ase$spliceType)),"\n")


## 拆解统计
level1 <- strsplit(level,split=",")[[1]][1]
level2 <- strsplit(level,split=",")[[1]][2]

# def function
stat4rmats <- function(ase,level1,level2){
  # psi
  ss1 <- data.frame(t(data.frame(lapply(as.character(ase[,"IncLevel1"]), function(x){
    y <- as.numeric(unlist(strsplit(x,',')[[1]][1:n]))
  }))))
  names(ss1) <- paste("PSI",level1,seq(1,2,1),sep="_")
  
  ss2 <- data.frame(t(data.frame(lapply(as.character(ase[,"IncLevel2"]), function(x){
    y <- as.numeric(unlist(strsplit(x,',')[[1]][1:n]))
  }))))
  names(ss2) <- paste("PSI",level2,seq(1,2,1),sep="_")
  
  ase <- cbind(ase,ss1,ss2)
  ase$IncLevel1 <- round(apply(ss1,1,mean),3)
  ase$IncLevel2 <- round(apply(ss2,1,mean),3)
  
  
  # IJC+SJC>10 per sample
  s1ijc <- data.frame(t(data.frame(lapply(as.character(ase[,"IJC_SAMPLE_1"]), function(x){
    y <- as.numeric(unlist(strsplit(x,',')[[1]][1:n]))
  }))))
  
  s1sjc <- data.frame(t(data.frame(lapply(as.character(ase[,"SJC_SAMPLE_1"]), function(x){
    y <- as.numeric(unlist(strsplit(x,',')[[1]][1:n]))
  }))))
  
  s2ijc <- data.frame(t(data.frame(lapply(as.character(ase[,"IJC_SAMPLE_2"]), function(x){
    y <- as.numeric(unlist(strsplit(x,',')[[1]][1:n]))
  }))))
  
  s2sjc <- data.frame(t(data.frame(lapply(as.character(ase[,"SJC_SAMPLE_2"]), function(x){
    y <- as.numeric(unlist(strsplit(x,',')[[1]][1:n]))
  }))))
  
  names(s1ijc) <- paste("IJC",level1,seq(1,2,1),sep="_")
  names(s2ijc) <- paste("IJC",level2,seq(1,2,1),sep="_")
  names(s1sjc) <- paste("SJC",level1,seq(1,2,1),sep="_")
  names(s2sjc) <- paste("SJC",level2,seq(1,2,1),sep="_")
  
  # stat reads across samples
  s1sijc <- s1ijc+s1sjc
  names(s1sijc) <- paste("SIJC",level1,seq(1,2,1),sep="_")
  s2sijc <- s2ijc+s2sjc
  names(s2sijc) <- paste("SIJC",level2,seq(1,2,1),sep="_")
  
  s12sijc <- data.frame(s1sijc,s2sijc)
  s12sijc$pctR10 <- round(apply(data.frame(s12sijc[,1:ncol(s12sijc)]>10),1,sum)/(2*n),2)
  s12sijc$pctR20 <- round(apply(data.frame(s12sijc[,1:ncol(s12sijc)]>20),1,sum)/(2*n),2)
  
  ase <- data.frame(ase,s1ijc,s1sjc,s2ijc,s2sjc,s12sijc)
  
  # IJC+SJC>10 per group
  ase$IJC_SAMPLE_1 <- round(as.numeric(lapply(
    as.character(ase[,"IJC_SAMPLE_1"]), function(x){
      y <- mean(as.numeric(unlist(strsplit(x,',')[[1]])))}
  )),1)

  ase$SJC_SAMPLE_1 <- round(as.numeric(lapply(
    as.character(ase[,"SJC_SAMPLE_1"]), function(x){
      y <- mean(as.numeric(unlist(strsplit(x,',')[[1]])))}
  )),1)

  ase$IJC_SAMPLE_2 <- round(as.numeric(lapply(
    as.character(ase[,"IJC_SAMPLE_2"]), function(x){
      y <- mean(as.numeric(unlist(strsplit(x,',')[[1]])))}
  )),1)

  ase$SJC_SAMPLE_2 <- round(as.numeric(lapply(
    as.character(ase[,"SJC_SAMPLE_2"]), function(x){
      y <- mean(as.numeric(unlist(strsplit(x,',')[[1]])))}
  )),1)
  
  return(ase)
}

ase <- stat4rmats(ase,level1,level2)

write.table(ase,file=paste(outPref,"all_ase_SJC_IJC_SIJC_PSI.txt",sep=""),sep="\t",row.names = FALSE,quote = FALSE)


## filtering by cutoff read counts
## total ase with sjc+ijc > cutoff
# aseTotal <- ase[ase$IJC_SAMPLE_1+ase$SJC_SAMPLE_1 > cutoff  & ase$IJC_SAMPLE_2+ase$SJC_SAMPLE_2 > cutoff,]
aseTotal <- ase[ase$pctR10==1,]
aseTotal <- aseTotal[!is.na(aseTotal$IncLevel1) & !is.na(aseTotal$IncLevel2),]

## filtering of novel ase
## nc or kd, max(min(ncSJC,ncIJC),min(kdSJC,kdIJC)) >10
if (novel==T){
  aseTotal$ncJCmin <- apply(aseTotal[,grepl("SJC_nc",names(aseTotal)) | grepl("IJC_nc",names(aseTotal))],1,min)
  aseTotal$kdJCmin <- apply(aseTotal[,grepl("SJC_kd",names(aseTotal)) | grepl("IJC_kd",names(aseTotal))],1,min)
  aseTotal$maxJCmin <- apply(data.frame(aseTotal[,c("ncJCmin","kdJCmin")]),1,max)
  aseTotal2 <- aseTotal[aseTotal$novel=="known" | aseTotal$maxJCmin>20,]
  aseTotal <- aseTotal[aseTotal$novel=="known" | aseTotal$maxJCmin>10,]
}

names <- names(ase)[grepl("SJC",names(ase)) | grepl("IJC",names(ase)) | grepl("SIJC",names(ase)) | grepl("PSI",names(ase))]
aseTotal <- aseTotal[,!names(aseTotal)%in%names]
aseTotal2 <- aseTotal2[,!names(aseTotal2)%in%names]
write.table(aseTotal,file=paste(outPref,"total_ase_minJC_10.txt",sep=""),sep="\t",row.names = FALSE,quote = FALSE)
write.table(aseTotal2,file=paste(outPref,"total_ase_minJC_20.txt",sep=""),sep="\t",row.names = FALSE,quote = FALSE)

cat(paste("Total ASEs after filtered with reads more than",cutoff,":",sep=""))
table(aseTotal$novel,aseTotal$spliceType)
cat("Sum:")
cat(sum(table(aseTotal$spliceType)),"\n")


# 筛选 sig ase
aseSig <- aseTotal[abs(aseTotal$IncLevelDifference)>0.05 & aseTotal$FDR<0.05,]
write.table(aseSig,file=paste(outPref,"total_ase_sig.txt",sep=""),sep="\t",row.names = FALSE,quote = FALSE)


cat("Significant DASEs(abs Diff > 0.05 & FDR < 0.05):")
table(aseSig$spliceType,aseSig$novel)
cat("Sum:")
cat(sum(table(aseSig$spliceType)),"\n")

## bed files
# bedName <- c("ID","geneSymbol","spliceType","chr","eventStart","eventEnd","strand","altStart","altEnd","exonNum","exonLength","exonStart")
# aseBed <- rbind(a3ss[,bedName], a5ss[,bedName], mxe[,bedName], se[,bedName], ri[,bedName])
# aseBed$idType <- paste(aseBed$geneSymbol,aseBed$spliceType,aseBed$ID,sep="_")
# aseBed <- aseBed[aseBed$idType%in%aseTotal$idType,]
# aseBed$score <- 0
# aseBed$altColor <- ifelse(aseBed$idType%in%aseSig$idType,
#                           ifelse(aseBed$idType%in%aseSigs$idType,"255,140,0","34,139,34"),
#                           "207,207,207")
# 
# if (novel==T){
#   aseBed$altColor <- ifelse(aseBed$idType%in%aseSig$idType,
#                             ifelse(aseBed$idType%in%aseTotal$novel,"255,140,0","34,139,34"),
#                             "207,207,207")
# }
# 
# bedName <- c("chr","eventStart","eventEnd","idType","score","strand","altStart","altEnd","altColor","exonNum","exonLength","exonStart")
# aseBed <- aseBed[aseBed$idType%in%aseTotal$idType,bedName]
# aseBedsig <- aseBed[aseBed$idType%in%aseSig$idType,]
# 
# write.table(aseBed,file=paste(outPref,"total_ase_ucsc.bed",sep=""),sep="\t",row.names = FALSE,quote = FALSE,col.names = FALSE)
# write.table(aseBedsig,file=paste(outPref,"sig_ase_ucsc.bed",sep=""),sep="\t",row.names = FALSE,quote = FALSE,col.names = FALSE)


## sig ase for sahami-plot
a3ssSig <- a3ss[paste(a3ss$geneSymbol, a3ss$spliceType,a3ss$ID,sep="_")%in%aseSig$idType,names(a3ss)!="spliceType"]
a5ssSig <- a5ss[paste(a5ss$geneSymbol, a5ss$spliceType,a5ss$ID,sep="_")%in%aseSig$idType,names(a5ss)!="spliceType"]
mxeSig <- mxe[paste(mxe$geneSymbol, mxe$spliceType,mxe$ID,sep="_")%in%aseSig$idType,names(mxe)!="spliceType"]
seSig <- se[paste(se$geneSymbol, se$spliceType,se$ID,sep="_")%in%aseSig$idType,names(se)!="spliceType"]
riSig <- ri[paste(ri$geneSymbol, ri$spliceType,ri$ID,sep="_")%in%aseSig$idType,names(ri)!="spliceType"]

write.table(a3ssSig, paste(outPref,"sig.A3SS.MATS.JC.txt",sep=""),sep="\t", row.names = FALSE,quote = FALSE)
write.table(a5ssSig, paste(outPref,"sig.A5SS.MATS.JC.txt",sep=""),sep="\t", row.names = FALSE,quote = FALSE)
write.table(mxeSig, paste(outPref,"sig.MXE.MATS.JC.txt",sep=""),sep="\t", row.names = FALSE,quote = FALSE)
write.table(seSig, paste(outPref,"sig.SE.MATS.JC.txt",sep=""),sep="\t", row.names = FALSE,quote = FALSE)
write.table(riSig, paste(outPref,"sig.RI.MATS.JC.txt",sep=""),sep="\t", row.names = FALSE,quote = FALSE)


## 以下合并所有数据后, 整体分析！
## splice score ----
# UP/DOWN stream, 根据根据基因组位置，统一按照正义链5‘->3',正反链没有方向差异
# splice site, 有5'与3'之分, 因此正反链有方向差异，按各自的5'->3'

# # A3
# a3ss$idType <- paste(a3ss$geneSymbol,a3ss$spliceType,a3ss$ID,sep="_")
# a3_1st3ss <- a3ss
# a3_1st3ss$start0 <- ifelse(a3_1st3ss$strand=="+",a3_1st3ss$longExonStart_0base-20,a3_1st3ss$longExonEnd-3)
# a3_1st3ss$end <- ifelse(a3_1st3ss$strand=="+",a3_1st3ss$longExonStart_0base+3,a3_1st3ss$longExonEnd+20)
# a3_1st3ss$ssType <- "1st3ss"
# 
# a3_2nd3ss <- a3ss
# a3_2nd3ss$start0 <- ifelse(a3_2nd3ss$strand=="+",a3_2nd3ss$shortES-20,a3_2nd3ss$shortEE-3)
# a3_2nd3ss$end <- ifelse(a3_2nd3ss$strand=="+",a3_2nd3ss$shortES+3,a3_2nd3ss$shortEE+20)
# a3_2nd3ss$ssType <- "2nd3ss"
# 
# # A5
# a5ss$idType <- paste(a5ss$geneSymbol,a5ss$spliceType,a5ss$ID,sep="_")
# a5_1st5ss <- a5ss
# a5_1st5ss$start0 <- ifelse(a5_1st5ss$strand=="+",a5_1st5ss$shortEE-3,a5_1st5ss$shortES-6)
# a5_1st5ss$end <- ifelse(a5_1st5ss$strand=="+",a5_1st5ss$shortEE+6,a5_1st5ss$shortES+3)
# a5_1st5ss$ssType <- "1st5ss"
# 
# a5_2nd5ss <- a5ss
# a5_2nd5ss$start0 <- ifelse(a5_2nd5ss$strand=="+",a5_2nd5ss$longExonEnd-3,a5_2nd5ss$longExonStart_0base-6)
# a5_2nd5ss$end <- ifelse(a5_2nd5ss$strand=="+",a5_2nd5ss$longExonEnd+6,a5_2nd5ss$longExonStart_0base+3)
# a5_2nd5ss$ssType <- "2nd5ss"
# 
# # SE
# # 3ss
# se$idType <- paste(se$geneSymbol,se$spliceType,se$ID,sep="_")
# 
# se_1st3ss <- se
# se_1st3ss$start0 <- ifelse(se_1st3ss$strand=="+",se_1st3ss$exonStart_0base-20,se_1st3ss$exonEnd-3)
# se_1st3ss$end <- ifelse(se_1st3ss$strand=="+",se_1st3ss$exonStart_0base+3,se_1st3ss$exonEnd+20)
# se_1st3ss$ssType <- "1st3ss"
# 
# se_2nd3ss <- se
# se_2nd3ss$start0 <- ifelse(se_2nd3ss$strand=="+",se_2nd3ss$downstreamES-20,se_2nd3ss$upstreamEE-3)
# se_2nd3ss$end <- ifelse(se_2nd3ss$strand=="+",se_2nd3ss$downstreamES+3,se_2nd3ss$upstreamEE+20)
# se_2nd3ss$ssType <- "2nd3ss"
# 
# # 5ss
# se_1st5ss <- se
# se_1st5ss$start0 <- ifelse(se_1st5ss$strand=="+",se_1st5ss$upstreamEE-3,se_1st5ss$downstreamES-6)
# se_1st5ss$end <- ifelse(se_1st5ss$strand=="+",se_1st5ss$upstreamEE+6,se_1st5ss$downstreamES+3)
# se_1st5ss$ssType <- "1st5ss"
# 
# se_2nd5ss <- se
# se_2nd5ss$start0 <- ifelse(se_2nd5ss$strand=="+",se_2nd5ss$exonEnd-3,se_2nd5ss$exonStart_0base-6)
# se_2nd5ss$end <- ifelse(se_2nd5ss$strand=="+",se_2nd5ss$exonEnd+6,se_2nd5ss$exonStart_0base+3)
# se_2nd5ss$ssType <- "2nd5ss"
# 
# 
# # MXE
# # 3ss
# mxe$idType <- paste(mxe$geneSymbol,mxe$spliceType,mxe$ID,sep="_")
# 
# mxe_1st3ss <- mxe
# mxe_1st3ss$start0 <- ifelse(mxe_1st3ss$strand=="+",mxe_1st3ss$X1stExonStart_0base -20,mxe_1st3ss$X2ndExonEnd-3)
# mxe_1st3ss$end <- ifelse(mxe_1st3ss$strand=="+",mxe_1st3ss$X1stExonStart_0base+3,mxe_1st3ss$X2ndExonEnd+20)
# mxe_1st3ss$ssType <- "1st3ss"
# 
# mxe_2nd3ss <- mxe
# mxe_2nd3ss$start0 <- ifelse(mxe_2nd3ss$strand=="+",mxe_2nd3ss$X2ndExonStart_0base-20,mxe_2nd3ss$X1stExonEnd-3)
# mxe_2nd3ss$end <- ifelse(mxe_2nd3ss$strand=="+",mxe_2nd3ss$X2ndExonStart_0base+3,mxe_2nd3ss$X1stExonEnd+20)
# mxe_2nd3ss$ssType <- "2nd3ss"
# 
# mxe_3rd3ss <- mxe
# mxe_3rd3ss$start0 <- ifelse(mxe_3rd3ss$strand=="+",mxe_3rd3ss$downstreamES-20,mxe_3rd3ss$upstreamEE-3)
# mxe_3rd3ss$end <- ifelse(mxe_3rd3ss$strand=="+",mxe_3rd3ss$downstreamES+3,mxe_3rd3ss$upstreamEE+20)
# mxe_3rd3ss$ssType <- "3rd3ss"
# 
# 
# # 5ss
# mxe_1st5ss <- mxe
# mxe_1st5ss$start0 <- ifelse(mxe_1st5ss$strand=="+",mxe_1st5ss$upstreamEE-3,mxe_1st5ss$downstreamES-6)
# mxe_1st5ss$end <- ifelse(mxe_1st5ss$strand=="+",mxe_1st5ss$upstreamEE+6,mxe_1st5ss$downstreamES+3)
# mxe_1st5ss$ssType <- "1st5ss"
# 
# mxe_2nd5ss <- mxe
# mxe_2nd5ss$start0 <- ifelse(mxe_2nd5ss$strand=="+",mxe_2nd5ss$X1stExonEnd-3,mxe_2nd5ss$X2ndExonStart_0base-6)
# mxe_2nd5ss$end <- ifelse(mxe_2nd5ss$strand=="+",mxe_2nd5ss$X1stExonEnd+6,mxe_2nd5ss$X2ndExonStart_0base+3)
# mxe_2nd5ss$ssType <- "2nd5ss"
# 
# mxe_3rd5ss <- mxe
# mxe_3rd5ss$start0 <- ifelse(mxe_3rd5ss$strand=="+",mxe_3rd5ss$X2ndExonEnd-3,mxe_3rd5ss$X1stExonStart_0base-6)
# mxe_3rd5ss$end <- ifelse(mxe_3rd5ss$strand=="+",mxe_3rd5ss$X2ndExonEnd+6,mxe_3rd5ss$X1stExonStart_0base+3)
# mxe_3rd5ss$ssType <- "3rd5ss"
# 
# 
# # RI
# # Note: riExon include upstream, intron and dowstream
# # 3ss
# ri$idType <- paste(ri$geneSymbol,ri$spliceType,ri$ID,sep="_")
# 
# ri_1st3ss <- ri
# ri_1st3ss$start0 <- ifelse(ri_1st3ss$strand=="+",ri_1st3ss$riExonStart_0base-20,ri_1st3ss$riExonEnd-3)
# ri_1st3ss$end <- ifelse(ri_1st3ss$strand=="+",ri_1st3ss$riExonStart_0base+3,ri_1st3ss$riExonEnd+20)
# ri_1st3ss$ssType <- "1st3ss"
# 
# ri_2nd3ss <- ri
# ri_2nd3ss$start0 <- ifelse(ri_2nd3ss$strand=="+",ri_2nd3ss$downstreamES-20,ri_2nd3ss$upstreamEE-3)
# ri_2nd3ss$end <- ifelse(ri_2nd3ss$strand=="+",ri_2nd3ss$downstreamES+3,ri_2nd3ss$upstreamEE+20)
# ri_2nd3ss$ssType <- "2nd3ss"
# 
# # 5ss
# ri_1st5ss <- ri
# ri_1st5ss$start0 <- ifelse(ri_1st5ss$strand=="+",ri_1st5ss$upstreamEE-3,ri_1st5ss$downstreamES-6)
# ri_1st5ss$end <- ifelse(ri_1st5ss$strand=="+",ri_1st5ss$upstreamEE+6,ri_1st5ss$downstreamES+3)
# ri_1st5ss$ssType <- "1st5ss"
# 
# ri_2nd5ss <- ri
# ri_2nd5ss$start0 <- ifelse(ri_2nd5ss$strand=="+",ri_2nd5ss$riExonEnd-3,ri_2nd5ss$riExonStart_0base-6)
# ri_2nd5ss$end <- ifelse(ri_2nd5ss$strand=="+",ri_2nd5ss$riExonEnd+6,ri_2nd5ss$riExonStart_0base+3)
# ri_2nd5ss$ssType <- "2nd5ss"
# 
# # 合并
# col = c("chr","start0","end","idType","ssType","strand")
# ss3 <- rbind(a3_1st3ss[,col],a3_2nd3ss[,col],
#              se_1st3ss[,col],se_2nd3ss[,col],
#              mxe_1st3ss[,col],mxe_2nd3ss[,col],mxe_3rd3ss[,col],
#              ri_1st3ss[,col],ri_2nd3ss[,col])
# 
# ss5 <- rbind(a5_1st5ss[,col],a5_2nd5ss[,col],
#              se_1st5ss[,col],se_2nd5ss[,col],
#              mxe_1st5ss[,col],mxe_2nd5ss[,col],mxe_3rd5ss[,col],
#              ri_1st5ss[,col],ri_2nd5ss[,col])
# 
# #ss3$chr <- gsub("chr","",ss3$chr)
# #ss5$chr <- gsub("chr","",ss5$chr)
# 
# write.table(ss3,file=paste(outPref,"_splice_site_3ss_location_for_score.bed",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote = FALSE)
# write.table(ss5,file=paste(outPref,"_splice_site_5ss_location_for_score.bed",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote = FALSE)
# 
# 
# 
# ## GC content ----
# # 应该包括上下游外显子情况
# 
# # A3
# a3_5intron <- a3ss
# a3_5intron$start0 <- ifelse(a3_5intron$strand=="+",a3_5intron$flankingEE,a3_5intron$longExonEnd)
# a3_5intron$end <- ifelse(a3_5intron$strand=="+",a3_5intron$longExonStart_0base,a3_5intron$flankingES)
# a3_5intron$exonType <- "a3_5intron"
# 
# a3_exonSE <- a3ss
# a3_exonSE$start0 <- ifelse(a3_exonSE$strand=="+",a3_exonSE$longExonStart_0base,a3_exonSE$shortEE)
# a3_exonSE$end <- ifelse(a3_exonSE$strand=="+",a3_exonSE$shortES,a3_exonSE$longExonEnd)
# a3_exonSE$exonType <- "a3_exonSE"
# 
# a3_3exon <- a3ss
# a3_3exon$start0 <- ifelse(a3_3exon$strand=="+",a3_3exon$shortES,a3_3exon$longExonStart_0base)
# a3_3exon$end <- ifelse(a3_3exon$strand=="+",a3_3exon$longExonEnd,a3_3exon$shortEE)
# a3_3exon$exonType <- "a3_3exon"
# 
# a3_gc <- rbind(a3_5intron,a3_exonSE,a3_3exon)
# 
# 
# # A5
# a5_5exon <- a5ss
# a5_5exon$start0 <- ifelse(a5_5exon$strand=="+",a5_5exon$longExonStart_0base,a5_5exon$shortES)
# a5_5exon$end <- ifelse(a5_5exon$strand=="+",a5_5exon$shortEE,a5_5exon$longExonEnd)
# a5_5exon$exonType <- "a5_5exon"
# 
# a5_exonSE <- a5ss
# a5_exonSE$start0 <- ifelse(a5_exonSE$strand=="+",a5_exonSE$shortEE,a5_exonSE$longExonStart_0base)
# a5_exonSE$end <- ifelse(a5_exonSE$strand=="+",a5_exonSE$longExonEnd,a5_exonSE$shortES)
# a5_exonSE$exonType <- "a5_exonSE"
# 
# a5_3intron <- a5ss
# a5_3intron$start0 <- ifelse(a5_3intron$strand=="+",a5_3intron$longExonEnd,a5_3intron$flankingEE)
# a5_3intron$end <- ifelse(a5_3intron$strand=="+",a5_3intron$flankingES,a5_3intron$longExonStart_0base)
# a5_3intron$exonType <- "a5_3intron"
# 
# a5_gc <- rbind(a5_5exon, a5_exonSE, a5_3intron)
# 
# 
# # SE
# se_5intron <- se
# se_5intron$start0 <- ifelse(se_5intron$strand=="+",se_5intron$upstreamEE,se_5intron$exonEnd)
# se_5intron$end <- ifelse(se_5intron$strand=="+",se_5intron$exonStart_0base,se_5intron$downstreamES)
# se_5intron$exonType <- "se_5intron"
# 
# se_exonSE <- se
# se_exonSE$start0 <- se_exonSE$exonStart_0base
# se_exonSE$end <- se_exonSE$exonEnd
# se_exonSE$exonType <- "se_exonSE"
# 
# se_3intron <- se
# se_3intron$start0 <- ifelse(se_3intron$strand=="+",se_3intron$exonEnd,se_3intron$upstreamEE)
# se_3intron$end <- ifelse(se_3intron$strand=="+",se_3intron$downstreamES,se_3intron$exonStart_0base)
# se_3intron$exonType <- "se_3intron"
# 
# se_gc <- rbind(se_5intron,se_exonSE,se_3intron)
# 
# 
# # MXE
# mxe_5intron <- mxe
# mxe_5intron$start0 <- ifelse(mxe_5intron$strand=="+",mxe_5intron$upstreamEE,mxe_5intron$X2ndExonEnd)
# mxe_5intron$end <- ifelse(mxe_5intron$strand=="+",mxe_5intron$X1stExonStart_0base,mxe_5intron$downstreamES)
# mxe_5intron$exonType <- "mxe_5intron"
# 
# mxe_5exonSE <- mxe
# mxe_5exonSE$start0 <- ifelse(mxe_5exonSE$strand=="+",mxe_5exonSE$X1stExonStart_0base,mxe_5exonSE$X2ndExonStart_0base)
# mxe_5exonSE$end <- ifelse(mxe_5exonSE$strand=="+",mxe_5exonSE$X1stExonEnd,mxe_5exonSE$X2ndExonEnd)
# mxe_5exonSE$exonType <- "mxe_5exonSE"
# 
# mxe_intron <- mxe
# mxe_intron$start0 <- mxe_intron$X1stExonEnd
# mxe_intron$end <- mxe_intron$X2ndExonStart_0base
# mxe_intron$exonType <- "mxe_intron"
# 
# mxe_3exonSE <- mxe
# mxe_3exonSE$start0 <- ifelse(mxe_3exonSE$strand=="+",mxe_3exonSE$X2ndExonStart_0base,mxe_3exonSE$X1stExonStart_0base)
# mxe_3exonSE$end <- ifelse(mxe_3exonSE$strand=="+",mxe_3exonSE$X2ndExonEnd,mxe_3exonSE$X1stExonEnd)
# mxe_3exonSE$exonType <- "mxe_3exonSE"
# 
# mxe_3intron <- mxe
# mxe_3intron$start0 <- ifelse(mxe_3intron$strand=="+",mxe_3intron$X2ndExonEnd,mxe_3intron$upstreamEE)
# mxe_3intron$end <- ifelse(mxe_3intron$strand=="+",mxe_3intron$downstreamES,mxe_3intron$X1stExonStart_0base)
# mxe_3intron$exonType <- "mxe_3intron"
# 
# mxe_gc <- rbind(mxe_5intron,mxe_5exonSE,mxe_intron,mxe_3exonSE,mxe_3intron)
# 
# 
# # RI
# ri_5exon <- ri
# ri_5exon$start0 <- ifelse(ri_5exon$strand=="+",ri_5exon$riExonStart_0base,ri_5exon$downstreamES)
# ri_5exon$end <- ifelse(ri_5exon$strand=="+",ri_5exon$upstreamEE,ri_5exon$riExonEnd)
# ri_5exon$exonType <- "ri_5exon"
# 
# ri_intronRI <- ri
# ri_intronRI$start0 <- ri_intronRI$upstreamEE
# ri_intronRI$end <- ri_intronRI$downstreamES
# ri_intronRI$exonType <- "ri_intronRI"
# 
# ri_3exon <- ri
# ri_3exon$start0 <- ifelse(ri_3exon$strand=="+",ri_3exon$downstreamES,ri_3exon$riExonStart_0base)
# ri_3exon$end <- ifelse(ri_3exon$strand=="+",ri_3exon$riExonEnd,ri_3exon$upstreamEE)
# ri_3exon$exonType <- "ri_3exon"
# 
# ri_gc <- rbind(ri_5exon,ri_intronRI,ri_3exon)
# 
# # 合并
# col = c("chr","start0","end","idType","exonType","strand")
# gcc <- rbind(a3_gc[,col],a5_gc[,col],se_gc[,col],mxe_gc[,col],ri_gc[,col])
# #gcc$chr <- gsub("chr","",gcc$chr)
# 
# write.table(gcc,file=paste(outPref,"_splice_exon_location_for_gc.bed",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote = FALSE)
# 
# 
# 
# ## Exon sequence----
# ## 包括上下游外显子情况
# 
# # A3
# a3_5exon <- a3ss
# a3_5exon$start0 <- ifelse(a3_5exon$strand=="+",a3_5exon$flankingES,a3_5exon$flankingES)
# a3_5exon$end <- ifelse(a3_5exon$strand=="+",a3_5exon$flankingEE,a3_5exon$flankingEE)
# a3_5exon$exonType <- "a3_5exon"
# 
# a3_exonSE <- a3ss
# a3_exonSE$start0 <- ifelse(a3_exonSE$strand=="+",a3_exonSE$longExonStart_0base,a3_exonSE$shortEE)
# a3_exonSE$end <- ifelse(a3_exonSE$strand=="+",a3_exonSE$shortES,a3_exonSE$longExonEnd)
# a3_exonSE$exonType <- "a3_exonSE"
# 
# a3_3exon <- a3ss
# a3_3exon$start0 <- ifelse(a3_3exon$strand=="+",a3_3exon$shortES,a3_3exon$shortES)
# a3_3exon$end <- ifelse(a3_3exon$strand=="+",a3_3exon$shortEE,a3_3exon$shortEE)
# a3_3exon$exonType <- "a3_3exon"
# 
# a3_exon <- rbind(a3_5exon,a3_exonSE,a3_3exon)
# 
# 
# # A5
# a5_5exon <- a5ss
# a5_5exon$start0 <- ifelse(a5_5exon$strand=="+",a5_5exon$shortES,a5_5exon$shortES)
# a5_5exon$end <- ifelse(a5_5exon$strand=="+",a5_5exon$shortEE,a5_5exon$longExonEnd)
# a5_5exon$exonType <- "a5_5exon"
# 
# a5_exonSE <- a5ss
# a5_exonSE$start0 <- ifelse(a5_exonSE$strand=="+",a5_exonSE$shortEE,a5_exonSE$longExonStart_0base)
# a5_exonSE$end <- ifelse(a5_exonSE$strand=="+",a5_exonSE$longExonEnd,a5_exonSE$shortES)
# a5_exonSE$exonType <- "a5_exonSE"
# 
# a5_3exon <- a5ss
# a5_3exon$start0 <- ifelse(a5_3exon$strand=="+",a5_3exon$flankingES,a5_3exon$flankingES)
# a5_3exon$end <- ifelse(a5_3exon$strand=="+",a5_3exon$flankingEE,a5_3exon$flankingEE)
# a5_3exon$exonType <- "a5_3exon"
# 
# a5_exon <- rbind(a5_5exon, a5_exonSE, a5_3exon)
# 
# 
# # SE
# se_5exon <- se
# se_5exon$start0 <- ifelse(se_5exon$strand=="+",se_5exon$upstreamES,se_5exon$downstreamES)
# se_5exon$end <- ifelse(se_5exon$strand=="+",se_5exon$upstreamEE,se_5exon$downstreamEE)
# se_5exon$exonType <- "se_5exon"
# 
# se_exonSE <- se
# se_exonSE$start0 <- se_exonSE$exonStart_0base
# se_exonSE$end <- se_exonSE$exonEnd
# se_exonSE$exonType <- "se_exonSE"
# 
# se_3exon <- se
# se_3exon$start0 <- ifelse(se_3exon$strand=="+",se_3exon$downstreamES,se_3exon$upstreamES)
# se_3exon$end <- ifelse(se_3exon$strand=="+",se_3exon$downstreamEE,se_3exon$upstreamEE)
# se_3exon$exonType <- "se_3exon"
# 
# se_exon <- rbind(se_5exon,se_exonSE,se_3exon)
# 
# 
# # MXE
# mxe_5exon <- mxe
# mxe_5exon$start0 <- ifelse(mxe_5exon$strand=="+",mxe_5exon$upstreamES,mxe_5exon$downstreamES)
# mxe_5exon$end <- ifelse(mxe_5exon$strand=="+",mxe_5exon$upstreamEE,mxe_5exon$downstreamEE)
# mxe_5exon$exonType <- "mxe_5exon"
# 
# mxe_5exonSE <- mxe
# mxe_5exonSE$start0 <- ifelse(mxe_5exonSE$strand=="+",mxe_5exonSE$X1stExonStart_0base,mxe_5exonSE$X2ndExonStart_0base)
# mxe_5exonSE$end <- ifelse(mxe_5exonSE$strand=="+",mxe_5exonSE$X1stExonEnd,mxe_5exonSE$X2ndExonEnd)
# mxe_5exonSE$exonType <- "mxe_5exonSE"
# 
# mxe_3exonSE <- mxe
# mxe_3exonSE$start0 <- ifelse(mxe_3exonSE$strand=="+",mxe_3exonSE$X2ndExonStart_0base,mxe_3exonSE$X1stExonStart_0base)
# mxe_3exonSE$end <- ifelse(mxe_3exonSE$strand=="+",mxe_3exonSE$X2ndExonEnd,mxe_3exonSE$X1stExonEnd)
# mxe_3exonSE$exonType <- "mxe_3exonSE"
# 
# mxe_3exon <- mxe
# mxe_3exon$start0 <- ifelse(mxe_3exon$strand=="+",mxe_3exon$downstreamES,mxe_3exon$upstreamES)
# mxe_3exon$end <- ifelse(mxe_3exon$strand=="+",mxe_3exon$downstreamEE,mxe_3exon$upstreamEE)
# mxe_3exon$exonType <- "mxe_3exon"
# 
# mxe_exon <- rbind(mxe_5exon,mxe_5exonSE,mxe_3exonSE,mxe_3exon)
# 
# 
# # RI
# ri_5exon <- ri
# ri_5exon$start0 <- ifelse(ri_5exon$strand=="+",ri_5exon$upstreamES,ri_5exon$downstreamES)
# ri_5exon$end <- ifelse(ri_5exon$strand=="+",ri_5exon$upstreamEE,ri_5exon$downstreamEE)
# ri_5exon$exonType <- "ri_5exon"
# 
# ri_exonSE <- ri
# ri_exonSE$start0 <- ri_exonSE$upstreamEE
# ri_exonSE$end <- ri_exonSE$downstreamES
# ri_exonSE$exonType <- "ri_exonSE"
# 
# ri_3exon <- ri
# ri_3exon$start0 <- ifelse(ri_3exon$strand=="+",ri_3exon$downstreamES,ri_3exon$upstreamES)
# ri_3exon$end <- ifelse(ri_3exon$strand=="+",ri_3exon$downstreamEE,ri_3exon$upstreamEE)
# ri_3exon$exonType <- "ri_3exon"
# 
# ri_exon <- rbind(ri_5exon,ri_exonSE,ri_3exon)
# 
# # 合并
# col = c("chr","start0","end","idType","exonType","strand")
# exonAll <- rbind(a3_exon[,col],a5_exon[,col],se_exon[,col],mxe_exon[,col],ri_exon[,col])
# #exonAll$chr <- gsub("chr","",exonAll$chr)
# exonAll$idType <- paste(exonAll$idType,exonAll$exonType,sep="::")
# write.table(exonAll,file=paste(outPref,"_splice_exon_location_for_sequence.bed",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote = FALSE)
# 
# 
# 
# ########################
# # isoforms switching----
# # gtf processing
# library(dplyr)
# library(stringr)
# # gtf1: exon info
# gtf1 <- read.table(gtfExon,sep="\t")
# gtf1$start0 <- gtf1$start-1
# gtf1$exonLoc <- paste(gtf1$start0,gtf1$end,sep=":")
# gtf1 <- gtf1[order(gtf1$transcriptId,gtf1$start0),]
# # gtf2: exons in transcript
# gtf2 <- gtf1 %>% group_by(transcriptId,geneId) %>% summarise(transcriptLoc = str_c(exonLoc, collapse = "-"))
# # gtf3: transcript info
# gtf3 <- gtf1[!duplicated(gtf1$transcriptId),]
# 
# # a3ss
# a3ss$exonInUp <- ifelse(a3ss$strand=="+",paste(a3ss$flankingES,a3ss$flankingEE,sep=":"),paste(a3ss$longExonStart_0base,a3ss$longExonEnd,sep=":"))
# a3ss$exonInDn <- ifelse(a3ss$strand=="+",paste(a3ss$longExonStart_0base,a3ss$longExonEnd,sep=":"),paste(a3ss$flankingES,a3ss$flankingEE,sep=":"))
# a3ss$exonIn <- paste(a3ss$exonInUp,a3ss$exonInDn,sep="-")
# a3ss$exonInPart <- ifelse(a3ss$strand=="+",paste(a3ss$flankingEE,a3ss$longExonStart_0base,sep="-"),paste(a3ss$longExonEnd,a3ss$flankingES,sep="-"))
# 
# a3ss$exonOutUp <- ifelse(a3ss$strand=="+",paste(a3ss$flankingES,a3ss$flankingEE,sep=":"),paste(a3ss$shortES,a3ss$shortEE,sep=":"))
# a3ss$exonOutDn <- ifelse(a3ss$strand=="+",paste(a3ss$shortES,a3ss$shortEE,sep=":"),paste(a3ss$flankingES,a3ss$flankingEE,sep=":"))
# a3ss$exonOut <- paste(a3ss$exonOutUp,a3ss$exonOutDn,sep="-")
# a3ss$exonOutPart <- ifelse(a3ss$strand=="+",paste(a3ss$flankingEE,a3ss$shortES,sep="-"),paste(a3ss$shortEE,a3ss$flankingES,sep="-"))
# 
# 
# # a5ss
# a5ss$exonInUp <- ifelse(a5ss$strand=="+",paste(a5ss$longExonStart_0base,a5ss$longExonEnd,sep=":"),paste(a5ss$flankingES,a5ss$flankingEE,sep=":"))
# a5ss$exonInDn <- ifelse(a5ss$strand=="+",paste(a5ss$flankingES,a5ss$flankingEE,sep=":"),paste(a5ss$longExonStart_0base,a5ss$longExonEnd,sep=":"))
# a5ss$exonIn <- paste(a5ss$exonInUp,a5ss$exonInDn,sep="-")
# a5ss$exonInPart <- ifelse(a5ss$strand=="+",paste(a5ss$longExonEnd,a5ss$flankingES,sep="-"),paste(a5ss$flankingEE,a5ss$longExonStart_0base,sep="-"))
# 
# a5ss$exonOutUp <- ifelse(a5ss$strand=="+",paste(a5ss$shortES,a5ss$shortEE,sep=":"),paste(a5ss$flankingES,a5ss$flankingEE,sep=":"))
# a5ss$exonOutDn <- ifelse(a5ss$strand=="+",paste(a5ss$flankingES,a5ss$flankingEE,sep=":"),paste(a5ss$shortES,a5ss$shortEE,sep=":"))
# a5ss$exonOut <- paste(a5ss$exonOutUp,a5ss$exonOutDn,sep="-")
# a5ss$exonOutPart <- ifelse(a5ss$strand=="+",paste(a5ss$shortEE,a5ss$flankingES,sep="-"),paste(a5ss$flankingEE,a5ss$shortES,sep="-"))
# 
# 
# # se 假定 SE 单个外显子
# se$exonInUp <- paste(se$upstreamES,se$upstreamEE,sep=":")
# se$exonInCe <- paste(se$exonStart_0base,se$exonEnd,sep=":")
# se$exonInDn <- paste(se$downstreamES,se$downstreamEE,sep=":")
# se$exonIn <- paste(se$exonInUp,se$exonInCe,se$exonInDn,sep="-")
# se$exonInPart <- paste(se$upstreamEE,se$exonInCe,se$downstreamES,sep="-")
# 
# # if SE多个外显子，则
# # se$exonInPartL <- paste(se$upstreamEE,se$exonStart_0base,sep="-")
# # se$exonInPartR <- paste(se$exonEnd,se$downstreamES,sep="-")
# 
# se$exonOutUp <- paste(se$upstreamES,se$upstreamEE,sep=":")
# se$exonOutDn <- paste(se$downstreamES,se$downstreamEE,sep=":")
# se$exonOut <- paste(se$exonOutUp,se$exonOutDn,sep="-")
# se$exonOutPart <- paste(se$upstreamEE,se$downstreamES,sep="-")
# 
# 
# # ri
# ri$exonIn <- paste(ri$riExonStart_0base,ri$riExonEnd,sep=":")
# ri$exonInPart <- paste(ri$riExonStart_0base,ri$riExonEnd,sep=":") #
# 
# ri$exonOutUp <- paste(ri$upstreamES,ri$upstreamEE,sep=":")
# ri$exonOutDn <- paste(ri$downstreamES,ri$downstreamEE,sep=":")
# ri$exonOut <- paste(ri$exonOutUp,ri$exonOutDn,sep="-")
# ri$exonOutPart <- paste(ri$upstreamEE,ri$downstreamES,sep="-")
# 
# 
# # mxe 
# # 5'exon in named splice in
# mxe$exonIn1st <- paste(mxe$upstreamES,mxe$upstreamEE,sep=":")
# mxe$exonIn2nd <- paste(mxe$X1stExonStart_0base,mxe$X1stExonEnd,sep=":")
# mxe$exonIn3rd <- paste(mxe$X2ndExonStart_0base,mxe$X2ndExonEnd,sep=":")
# mxe$exonIn4th <- paste(mxe$downstreamES,mxe$downstreamEE,sep=":")
# 
# mxe$exonIn <- paste(mxe$exonIn1st,mxe$exonIn2nd,mxe$exonIn4th,sep="-")
# mxe$exonInPart <- paste(mxe$upstreamEE,mxe$exonIn2nd,mxe$downstreamES,sep="-")
# 
# mxe$exonOut <- paste(mxe$exonIn1st,mxe$exonIn3rd,mxe$exonIn4th,sep="-")
# mxe$exonOutPart <- paste(mxe$upstreamEE,mxe$exonIn3rd,mxe$downstreamES,sep="-")
# 
# # exon location, will merge with gtf1, for exon number in transcript
# a3ss$exonLoc <- paste(a3ss$longExonStart_0base,a3ss$longExonEnd,sep=":")
# a5ss$exonLoc <- paste(a5ss$longExonStart_0base,a5ss$longExonEnd,sep=":")
# se$exonLoc <- paste(se$exonStart_0base,se$exonEnd,sep=":")
# ri$exonLoc <- paste(ri$riExonStart_0base,ri$riExonEnd,sep=":")
# mxe$exonLoc <- paste(mxe$X1stExonStart_0base,mxe$X1stExonEnd,sep=":")
# 
# col = c("GeneID","idType","exonIn","exonInPart","exonOut","exonOutPart","exonLoc")
# aseAll <- rbind(a3ss[,col],a5ss[,col],se[,col],mxe[,col],ri[,col])
# 
# # exon in transcript
# exonInfo <- merge(aseAll,gtf2,by.x="GeneID",by.y="geneId",all.x=T)
# 
# exonIn <- exonInfo[str_detect(exonInfo$transcriptLoc,exonInfo$exonIn) | str_detect(exonInfo$transcriptLoc,exonInfo$exonInPart),]
# exonIn$partORwhole <- ifelse(str_detect(exonIn$transcriptLoc,exonIn$exonIn),"wholeIn","partIn")
# exonIn$exonType <- "In"
# 
# exonOut <- exonInfo[str_detect(exonInfo$transcriptLoc,exonInfo$exonOut) | str_detect(exonInfo$transcriptLoc,exonInfo$exonOutPart),]
# exonOut$partORwhole <- ifelse(str_detect(exonOut$transcriptLoc,exonOut$exonOut),"wholeOut","partOut")
# exonOut$exonType <- "Out"
# 
# 
# # merge with transcript info
# colname = c("transcriptId","idType","partORwhole","exonType")
# gtfname = c("transcriptBiotype","transcriptId","transcriptName")
# 
# exonIn <- merge(exonIn[,colname],gtf3[,gtfname],by="transcriptId",all.x=T)
# exonOut <- merge(exonOut[,colname],gtf3[,gtfname],by="transcriptId",all.x=T)
# 
# # ase associated isoform switching
# # partIn/out 指至少juction匹配, 不管上下游exon的边界
# # wholeIn/out 指上下游exon都匹配
# # splice in, for mxe means 1st splice exon splice in
# exonPartInIso <- 
#   exonIn %>% 
#   group_by(idType) %>% 
#   summarise(inTranscriptId = str_c(transcriptId, collapse = ","),
#             inTranscriptName = str_c(transcriptName, collapse = ","),
#             inTranscriptBiotype = str_c(transcriptBiotype, collapse = ","))
# 
# exonWholeInIso <- 
#   exonIn[exonIn$partORwhole=="wholeIn",] %>% 
#   group_by(idType) %>% 
#   summarise(inTranscriptId = str_c(transcriptId, collapse = ","),
#             inTranscriptName = str_c(transcriptName, collapse = ","),
#             inTranscriptBiotype = str_c(transcriptBiotype, collapse = ","))
# 
# # splice out
# exonPartOutIso <- 
#   exonOut %>% 
#   group_by(idType) %>% 
#   summarise(outTranscriptId = str_c(transcriptId, collapse = ","),
#             outTranscriptName = str_c(transcriptName, collapse = ","),
#             outTranscriptBiotype = str_c(transcriptBiotype, collapse = ","))
# 
# exonWholeOutIso <- 
#   exonOut[exonOut$partORwhole=="wholeOut",] %>% 
#   group_by(idType) %>% 
#   summarise(outTranscriptId = str_c(transcriptId, collapse = ","),
#             outTranscriptName = str_c(transcriptName, collapse = ","),
#             outTranscriptBiotype = str_c(transcriptBiotype, collapse = ","))
# 
# # exon number of ase exon
# aseAll$exonSE <- paste(aseAll$GeneID,aseAll$exonLoc,sep=":")
# gtf1$exonSE <- paste(gtf1$geneId,gtf1$exonLoc,sep=":")
# aseExon <- merge(aseAll,gtf1,by="exonSE",all.x=T)
# aseExon$transcriptNameExon <- paste(aseExon$transcriptName,"( E",aseExon$exonNumber," )",sep="")
# aseExon <- aseExon %>% group_by(idType,geneBiotype) %>% summarise(transcriptNameExon = str_c(transcriptNameExon, collapse = ","))
# 
# # 
# aseInfoName <- c("idType","spliceType","GeneID","geneSymbol","IncLevel1","IncLevel2","IncLevelDifference","FDR")
# isoTransPart <- merge(exonPartInIso,exonPartOutIso,by="idType",all=T)
# isoTransPart <- merge(aseExon,isoTransPart,by="idType",all=T)
# isoTransPart <- merge(ase[,aseInfoName],isoTransPart,by="idType",all.x=T)
# isoTransPartSig <- isoTransPart[isoTransPart$idType%in%aseSig$idType,]
# 
# isoTransWhole <- merge(exonWholeInIso,exonWholeOutIso,by="idType",all=T)
# isoTransWhole <- merge(aseExon,isoTransWhole,by="idType",all=T)
# isoTransWhole <- merge(ase[,aseInfoName],isoTransWhole,by="idType",all.x=T)
# isoTransWholeSig <- isoTransWhole[isoTransWhole$idType%in%aseSig$idType,]
# 
# write.table(isoTransPart,file=paste(outPref,"_isoform_switch_juction_used_only.txt",sep=""),sep="\t",row.names=F, quote = F)
# write.table(isoTransWhole,file=paste(outPref,"_isoform_switch_by_whole_exon_boundary.txt",sep=""),sep="\t",row.names=F, quote = F)
# write.table(isoTransPartSig,file=paste(outPref,"_isoform_switch_juction_used_only_sigASE.txt",sep=""),sep="\t",row.names=F, quote = FALSE)
# write.table(isoTransWholeSig,file=paste(outPref,"_isoform_switch_by_whole_exon_boundary_sigASE.txt",sep=""),sep="\t",row.names=F, quote = FALSE)


