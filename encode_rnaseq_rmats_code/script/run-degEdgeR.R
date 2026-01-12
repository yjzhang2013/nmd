## edgeR 差异分析
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
library(edgeR)

# Create a parser
p <- arg_parser("run DEG by edgeR")

# Add command line arguments
p <- add_argument(p, "--counts", help="input counts", type="character")
p <- add_argument(p, "--gtfGene", help="gene info from gtf", type="character")
p <- add_argument(p, "--l1", help="sample group level1:interference", type="character")
p <- add_argument(p, "--l2", help="sample group level2:control", type="character")
p <- add_argument(p, "--outPref", help="output prefix", type="character",default="")

# Parse the command line arguments
argv <- parse_args(p)

counts <- argv$counts
gtfGene <- argv$gtfGene
outPref <- argv$outPref
level1 <- argv$l1
level2 <- argv$l2

# counts <- "/home/u1357/encode/rnaseq/hepg2/merge/ENCSR000SKS.counts.merge.txt"
# gtfGene <- "/home/u1357/encode/refv19/gencode.v19.processed.geneInfo.txt"
# outPref <- "/home/u1357/encode/rnaseq/hepg2/deg/ENCSR000SKS"
# level1 <- "nc"
# level2 <- "kd"


cat("Running DEG by edgeR:",level1,"vs",level2)


# input counts matrix
fcounts <- read.table(counts,sep="\t",header = T)
names(fcounts) <- gsub(".counts","",names(fcounts))

groupList <- unlist(strsplit(names(fcounts),"_"))[seq(2,11,3)]

# geneInfo
geneInfo <- read.table(gtfGene,sep="\t",header = T)


# DGEList
dge <- DGEList(fcounts, group=groupList)

# cpmCutoff/lcpmCutoff calculate for plot
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
cpmCutoff <- 10/M + 2/L
lcpmCutoff <- log2(10/M + 2/L)

# auto-delete low expressed genes
# keep <- filterByExpr(dge,group=groupList)
# keep <- dge[keep,, keep.lib.sizes=FALSE]
cat("\nTotal genes before filtering: \n",nrow(dge),"\n")
# cat("Total genes after filtering: \n",nrow(keep),"\n")

# normalization by TMM
# nkeep <- calcNormFactors(keep, method="TMM")
nkeep <- calcNormFactors(dge, method="TMM")


# differential analysis take consider of batch effects
design <- model.matrix(~0 + factor(groupList,level=c(level1,level2)))
rownames(design)<-colnames(nkeep)
colnames(design)<-c(level1,level2)

# evaluate common、trended、tagwise
nkeep <- estimateGLMCommonDisp(nkeep,design)
nkeep <- estimateGLMTrendedDisp(nkeep, design)
nkeep <- estimateGLMTagwiseDisp(nkeep, design)

# normalized cpm of keep genes
ncpmk <- data.frame(round(cpm(nkeep),2))
names(ncpmk) <- paste(names(ncpmk),".ncpm",sep="")
ncpmk$geneID <- row.names(ncpmk)

# glmFit and glmLRT
fit <- glmFit(nkeep, design)
contrast <- c(1,-1)
lrt <- glmLRT(fit, contrast=contrast)

nrDEG=data.frame(topTags(lrt, n=nrow(dge)))
nrDEG$absLogFC <- abs(nrDEG$logFC)

nrDEG$geneID <- row.names(nrDEG)
nrDEG = merge(nrDEG,ncpmk,by="geneID")

nrDEG <- merge(geneInfo,nrDEG,by="geneID",all.y=T)
nrDEG=nrDEG[order(nrDEG$FDR),]
write.table(nrDEG, file=paste(outPref,"_diff_gene.txt",sep=""),sep="\t",row.names = F)

# tpm <- read.table(file="/home/u1357/encode/rnaseq/hepg2/merge/ENCSR000SKS.tpm.merge.txt",sep="\t")
# tpm$geneID <- row.names(tpm)
# cts <- read.table(file="/home/u1357/encode/rnaseq/hepg2/merge/ENCSR000SKS.counts.merge.txt",sep="\t")
# cts$geneID <- row.names(cts)




########------####
####     plot     #### to compare
########------####
# # # compare before and after
# library(RColorBrewer)
# nsamples <- ncol(dge)
# col <- brewer.pal(nsamples, "Paired")
# 
# pdf(paste(outPref,"filteringByexpression.pdf",sep=""),width=6,height=8)
# par(mfrow=c(3,2))
# 
# # 01 before/after filtering
# # before
# lcpm <- cpm(dge, log=TRUE)
# plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm[,1])$y)), las=2, main="", xlab="")
# title(main=paste("Before filtering \n n=",nrow(dge),sep=""), xlab="Log-cpm")
# abline(v=lcpmCutoff, lty=3)
# for (i in 2:nsamples){
#   den <- density(lcpm[,i])
#   lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", colnames(dge), text.col=col, bty="n")
# # after
# lcpmk <- cpm(keep, log=TRUE)
# plot(density(lcpmk[,1]), col=col[1], lwd=2, ylim=c(0,max(density(lcpm[,1])$y)), las=2, main="", xlab="")
# title(main=paste("After filtering \n n=",nrow(keep),sep=""), xlab="Log-cpm")
# abline(v=lcpmCutoff, lty=3)
# for (i in 2:nsamples){
# den <- density(lcpmk[,i])
# lines(den$x, den$y, col=col[i], lwd=2)
# }
# legend("topright", colnames(keep), text.col=col, bty="n")
# 
# 
# # 02 before/after normalization
# #par(mfrow=c(1,2))
# # before
# lcpmk <- cpm(keep, log=TRUE)
# boxplot(lcpmk, las=2, col=col, main="")
# title(main="Before normalization",ylab="Log-cpm")
# 
# # after
# nlcpmk <- cpm(nkeep, log=TRUE)
# boxplot(nlcpmk, las=2, col=col, main="")
# title(main="After normalization",ylab="Log-cpm")
# 
# 
# # 03 MDS clustering before/after filtering and normalization
# # before
# # lcpm <- cpm(dge, log=TRUE)
# plotMDS(lcpm,labels=rownames(design),col=col)
# title(main="Before filtering+normalization")
# # after
# #nlcpmk <- cpm(nkeep, log=TRUE)
# plotMDS(nlcpmk,labels=rownames(design),col=col)
# title(main="After filtering+normalization")
# dev.off()
# 
# 
# # plot BCV before/after filtering and normalization
# # evaluate common、trended、tagwise of dge
# dge <- estimateGLMCommonDisp(dge,design)
# dge <- estimateGLMTrendedDisp(dge, design)
# dge <- estimateGLMTagwiseDisp(dge, design)
# 
# tiff(paste(outPref,"plotBCV.tif",sep=""),res=300,units='in',width=7,height=4)
# par(mfrow=c(1,2))
# # before
# plotBCV(dge)
# title(main="Before filtering+nomalization",ylab="Log-cpm")
# 
# # after
# plotBCV(nkeep)
# title(main="After filtering+nomalization",ylab="Log-cpm")
# dev.off()

