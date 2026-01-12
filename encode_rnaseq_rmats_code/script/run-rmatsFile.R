# rmats file

library(argparser, quietly=TRUE)
library(dplyr)
library(stringr)

# Create a parser
p <- arg_parser("run rmatsFile")

# Add command line arguments
p <- add_argument(p, "--fileDir", help="file directory", type="character")
p <- add_argument(p, "--bamDir", help="bam file directory", type="character")
p <- add_argument(p, "--l1", help="sample group level1:interference", type="character")
p <- add_argument(p, "--l2", help="sample group level2:control", type="character")


# Parse the command line arguments
argv <- parse_args(p)

fileDir <- argv$fileDir
bamPref <- argv$bamDir
level1 <- argv$l1
level2 <- argv$l2

sample <- read.table(paste(fileDir,"sampleGroup.txt",sep=""),sep="",header=T)
sample$level <- ifelse(sample$group==level1,"b1","b2")
sample$bamFile <- paste(bamPref,sample$sample,".Aligned.sortedByCoord.out.bam",sep="")
bx <- sample %>% group_by(group) %>% summarise(bamDir=str_c(bamFile,collapse=","))
write.table(bx[bx$group==level1,2],file=paste(fileDir,"b1.txt",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote = FALSE)
write.table(bx[bx$group==level2,2],file=paste(fileDir,"b2.txt",sep=""),sep="\t",row.names = FALSE,col.names=FALSE,quote = FALSE)
write.table(sample,file=paste(fileDir,"rmatsSampleGroup.txt",sep=""),sep="\t",row.names = FALSE,quote = FALSE)

cat("Set", level1,"as b1.txt")
cat("Set", level2,"as b2.txt")

