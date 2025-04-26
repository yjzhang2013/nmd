#########
######### exon length, splice site strength, and GC content
#########
tcga <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all_txAnno_repeatAnno.rds")
outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/siteScore/"
dir.create(outDir)

gcContent(tcga, 250, 50, outDir)
siteScoring(tcga, outDir)


siteScoring <- function(tcga, outDir){
  message("tcga region process...\n")
  tcga <- tcga %>% mutate(geneName = limma::strsplit2(idMap, "\\|")[,1],
                          eventRegion = limma::strsplit2(idMap, "\\|")[,3])
  
  tcga <- cbind(tcga, data.frame(limma::strsplit2(tcga$eventRegion, "[:, -]")))
  for (i in paste0("X",seq(8))){
    tcga[,i] <- as.numeric(tcga[,i])
  }
  
  elen <- 3
  ilen <- 6
  slen <- 20
  
  a3 <- tcga[tcga$spliceType=="A3",]
  a3 <- a3 %>% mutate(s5start1=ifelse(strand=="+", X2 - elen, X5 - ilen),
                      s5end1=ifelse(strand=="+", X2 + ilen, X5 + elen),
                      
                      s3start1=ifelse(strand=="+", X3 - slen, X4 - elen),
                      s3end1=ifelse(strand=="+", X3 + elen, X4 + slen),
                      
                      s3start2=ifelse(strand=="+", X5 - slen, X2 - elen),
                      s3end2=ifelse(strand=="+", X5 + elen, X2 + slen))
  
  ## get s1~s6 (A3, A5, ES, IR, ME), start & end
  a3 <- a3[,c("chr", "idMap", "strand",names(a3)[grepl("s[3,5].*[1,2]",names(a3))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start", regionType), "start", "end"),
                  regionType = gsub("start", ".", gsub("end", ".", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  a5 <- tcga[tcga$spliceType=="A5",]
  a5 <- a5 %>% mutate(s5start1=ifelse(strand=="-", X5 - ilen, X2 - elen),
                      s5end1=ifelse(strand=="-", X5 + elen, X2 + ilen),
                      
                      s5start2=ifelse(strand=="-", X3 - ilen, X4 - elen),
                      s5end2=ifelse(strand=="-", X3 + elen, X4 + ilen),
                      
                      s3start1=ifelse(strand=="-", X2 - elen, X5 - slen),
                      s3end1=ifelse(strand=="-", X2 + slen, X5 + elen))
  
  a5 <- a5[,c("chr", "idMap", "strand",names(a5)[grepl("s[3,5]{1,2}.*[1,2]",names(a5))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start", regionType), "start", "end"),
                  regionType = gsub("start", ".", gsub("end", ".", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  se <- tcga[tcga$spliceType=="ES",]
  se <- se %>% mutate(s5start1=ifelse(strand=="+", X2 - elen, X5 - ilen),
                      s5end1=ifelse(strand=="+", X2 + ilen, X5 + elen),
                      
                      s3start1=ifelse(strand=="+", X3 - slen, X4 - elen),
                      s3end1=ifelse(strand=="+", X3 + elen, X4 + slen),
                      
                      s5start2=ifelse(strand=="+", X4 - elen, X3 - ilen),
                      s5end2=ifelse(strand=="+", X4 + ilen, X3 + elen),
                      
                      s3start2=ifelse(strand=="+", X5 - slen, X2 - elen),
                      s3end2=ifelse(strand=="+", X5 + elen, X2 + slen))
  
  se <- se[,c("chr", "idMap", "strand",names(se)[grepl("s[3,5]{1,2}.*[1,2]",names(se))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start", regionType), "start", "end"),
                  regionType = gsub("start", ".", gsub("end", ".", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  
  ri <- tcga[tcga$spliceType=="IR",]
  ri <- ri %>% mutate(s5start1=ifelse(strand=="+", X2 - elen, X5 - ilen),
                      s5end1=ifelse(strand=="+", X2 + ilen, X5 + elen),
                      
                      s3start1=ifelse(strand=="+", X5 - slen, X2 - elen),
                      s3end1=ifelse(strand=="+", X5 + elen, X2 + slen))
  
  ri <- ri[,c("chr", "idMap", "strand",names(ri)[grepl("s[3,5]{1,2}.*[1,2]",names(ri))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start", regionType), "start", "end"),
                  regionType = gsub("start", ".", gsub("end", ".", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  me <- tcga[tcga$spliceType=="ME",]
  me <- me %>% mutate(s5start1=ifelse(strand=="+", X2 - elen, X7 - ilen),
                      s5end1=ifelse(strand=="+", X2 + ilen, X7 + elen),
                      
                      s3start1=ifelse(strand=="+", X3 - slen, X6 - elen),
                      s3end1=ifelse(strand=="+", X3 + elen, X6 + slen),
                      
                      s5start2=ifelse(strand=="+", X4 - elen, X5 - ilen),
                      s5end2=ifelse(strand=="+", X4 + ilen, X5 + elen),
                      
                      s3start2=ifelse(strand=="+", X5 - slen, X4 - elen),
                      s3end2=ifelse(strand=="+", X5 + elen, X4 + slen),
                      
                      s5start3=ifelse(strand=="+", X6 - elen, X3 - ilen),
                      s5end3=ifelse(strand=="+", X6 + ilen, X3 + elen),
                      
                      s3start3=ifelse(strand=="+", X7 - slen, X2 - elen),
                      s3end3=ifelse(strand=="+", X7 + elen, X2 + slen))
  
  me <- me[,c("chr", "idMap", "strand",names(me)[grepl("s[3,5]{1,2}.*[1-3]",names(me))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start", regionType), "start", "end"),
                  regionType = gsub("start", ".", gsub("end", ".", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  tcga <- rbind(a3, a5, se, ri, me) %>% 
    dplyr::mutate(exonLen = end-start) %>% 
    dplyr::filter((grepl("^s3", regionType) & exonLen==23) | grepl("^s5", regionType) & exonLen==9) %>% 
    select(chr, start, end, idMap, regionType, strand)
  
  
  s3 <- tcga %>% 
    dplyr::filter(grepl("^s3.", regionType))
  
  s5 <- tcga %>% 
    dplyr::filter(grepl("^s5.", regionType))
  
  options(scipen = 7)
  fileNames3 <- paste0(outDir, "tcga_ase_sitescoring_sites3.bed")
  write.table(s3, file=fileNames3, sep="\t", row.names = F, col.names = F, quote = F)
  
  fileNames5 <- paste0(outDir, "tcga_ase_sitescoring_sites5.bed")
  write.table(s5, file=fileNames5, sep="\t", row.names = F, col.names = F, quote = F)
  
  ## get fasta
  message("bedtools getfasta...\n")
  genomeFa <- "/home/u1357/RNAseq/refv19/GRCh37.p13.genome.fa"
  
  ## s3
  outputs3 <- paste0(outDir, "tcga_ase_sitescoring_sites3.fasta")
  path <- "/pub/anaconda3/bin/"
  cmd <- paste0(path, "bedtools getfasta -fi ", genomeFa," -bed ", fileNames3," -name -s >", outputs3)
  system(cmd)
  
  
  ## s5
  outputs5 <- paste0(outDir, "tcga_ase_sitescoring_sites5.fasta")
  path <- "/pub/anaconda3/bin/"
  cmd <- paste0(path, "bedtools getfasta -fi ", genomeFa," -bed ", fileNames5," -name -s >", outputs5)
  system(cmd)
  
  
  scores3 <- gsub(".fasta", "_score.txt", outputs3)
  scores5 <- gsub(".fasta", "_score.txt", outputs5)
  setwd("/home/u1357/software/maxentscore/")
  cmds3 <- paste0("perl /home/u1357/software/maxentscore/score3.pl ", outputs3, " > ", scores3)
  cmds5 <- paste0("perl /home/u1357/software/maxentscore/score5.pl ", outputs5, " >", scores5)
  system(cmds3)
  system(cmds5)

  message("all done\n")
}

gcContent <- function(tcga, ilen, elen, outDir){
  # ilen <- 250
  # elen <- 50
  message("tcga region process...\n")
  tcga <- tcga %>% mutate(geneName = limma::strsplit2(idMap, "\\|")[,1],
                          eventRegion = limma::strsplit2(idMap, "\\|")[,3],
                          spliceType = limma::strsplit2(idMap, "\\|")[,2])
  
  tcga <- cbind(tcga, data.frame(limma::strsplit2(tcga$eventRegion, "[:, -]")))
  for (i in paste0("X",seq(8))){
    tcga[,i] <- as.numeric(tcga[,i])
  }
  
  
  a3 <- tcga[tcga$spliceType=="A3",]
  a3 <- a3 %>% mutate(r1start=ifelse(strand=="+", pmax(X1, X2 - elen), X5),
                      r1end=ifelse(strand=="+", X2, pmin(X5 + elen, X6)),
                      
                      r2start=ifelse(strand=="+", X2, pmax(X4, X5 - ilen)),
                      r2end=ifelse(strand=="+", pmin(X2 + ilen, X3), X5),
                      
                      r3start=ifelse(strand=="+", pmax(X2, X3 - ilen), X4),
                      r3end=ifelse(strand=="+", X3, pmin(X4 + ilen, X5)),
                      
                      r4start=ifelse(strand=="+", X3, pmax(X1, X4 - elen)),
                      r4end=ifelse(strand=="+", pmin(X3 + elen, X6), X4),
                      
                      r5start=ifelse(strand=="+", pmax(X2, X5 - ilen), X2),
                      r5end=ifelse(strand=="+", X5, pmin(X2 + ilen, X5)),
                      
                      r6start=ifelse(strand=="+", X5, pmax(X1, X2 - elen)),
                      r6end=ifelse(strand=="+", pmin(X5 + elen, X6), X2))
  
  ## get r1~r12 (A3, A5, ES, IR, ME), start & end                    
  a3 <- a3[,c("chr", "idMap", "strand",names(a3)[grepl("r[0-9]{1,2}[a-z]{3,5}",names(a3))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  a5 <- tcga[tcga$spliceType=="A5",]
  a5 <- a5 %>% mutate(r1start=ifelse(strand=="-", X5, pmax(X1, X2 - elen)),
                      r1end=ifelse(strand=="-", pmin(X5 + elen, X6), X2),
                      
                      r2start=ifelse(strand=="-", pmax(X2, X5 - ilen), X2),
                      r2end=ifelse(strand=="-", X5, pmin(X2 + ilen, X5)),
                      
                      r3start=ifelse(strand=="-", X3, pmax(X1, X4 - elen)),
                      r3end=ifelse(strand=="-", pmin(X3 + elen, X6), X4),
                      
                      r4start=ifelse(strand=="-", pmax(X2, X3 - ilen), X4),
                      r4end=ifelse(strand=="-", X3, pmin(X4 + ilen, X5)),
                      
                      r5start=ifelse(strand=="-", X2, pmax(X4, X5 - ilen)),
                      r5end=ifelse(strand=="-", pmin(X2 + ilen, X3), X5),
                      
                      r6start=ifelse(strand=="-", pmax(X1, X2 - elen), X5),
                      r6end=ifelse(strand=="-", X2, pmin(X5 + elen, X6)))
  
  a5 <- a5[,c("chr", "idMap", "strand",names(a5)[grepl("r[0-9]{1,2}[a-z]{3,5}",names(a5))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  se <- tcga[tcga$spliceType=="ES",]
  se <- se %>% mutate(r1start=ifelse(strand=="+", pmax(X1, X2 - elen), X5),
                      r1end=ifelse(strand=="+", X2, pmin(X5 + elen, X6)),
                      
                      r2start=ifelse(strand=="+", X2, pmax(X4, X5 - ilen)),
                      r2end=ifelse(strand=="+", pmin(X2 + ilen, X3), X5),
                      
                      r3start=ifelse(strand=="+", pmax(X2, X3 - ilen), X4),
                      r3end=ifelse(strand=="+", X3, pmin(X4 + ilen, X5)),
                      
                      r4start=ifelse(strand=="+", X3, pmax(X3, X4 - elen)),
                      r4end=ifelse(strand=="+", pmin(X3 + elen, X4), X4),
                      
                      r5start=ifelse(strand=="+", pmax(X3, X4 - elen), X3),
                      r5end=ifelse(strand=="+", X4, pmin(X3 + elen, X4)),
                      
                      r6start=ifelse(strand=="+", X4, pmax(X2, X3 - ilen)),
                      r6end=ifelse(strand=="+", pmin(X4 + ilen, X5), X3),
                      
                      r7start=ifelse(strand=="+", pmax(X4, X5 - ilen), X2),
                      r7end=ifelse(strand=="+", X5, pmin(X2 + ilen, X3)),
                      
                      r8start=ifelse(strand=="+", X5, pmax(X1, X2 - elen)),
                      r8end=ifelse(strand=="+", pmin(X5 + elen, X6), X2))
  
  se <- se[,c("chr", "idMap", "strand",names(se)[grepl("r[0-9]{1,2}[a-z]{3,5}",names(se))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  
  ri <- tcga[tcga$spliceType=="IR",]
  ri <- ri %>% mutate(r1start=ifelse(strand=="+", pmax(X1, X2 - elen), X5),
                      r1end=ifelse(strand=="+", X2, pmin(X5 + elen, X6)),
                      
                      r2start=ifelse(strand=="+", X2, pmax(X2, X5 - ilen)),
                      r2end=ifelse(strand=="+", pmin(X2 + ilen, X5), X5),
                      
                      r3start=ifelse(strand=="+", pmax(X2, X5 - ilen), X2),
                      r3end=ifelse(strand=="+", X5, pmin(X2 + ilen, X5)),
                      
                      r4start=ifelse(strand=="+", X5, pmax(X1, X2 - elen)),
                      r4end=ifelse(strand=="+", pmin(X5 + elen, X6), X2))
  
  ri <- ri[,c("chr", "idMap", "strand",names(ri)[grepl("r[0-9]{1,2}[a-z]{3,5}",names(ri))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  me <- tcga[tcga$spliceType=="ME",]
  me <- me %>% mutate(r1start=ifelse(strand=="+", pmax(X1, X2 - elen), X7),
                      r1end=ifelse(strand=="+", X2, pmin(X7 + elen, X8)),
                      
                      r2start=ifelse(strand=="+", X2, pmax(X6, X7 - ilen)),
                      r2end=ifelse(strand=="+", pmin(X2 + ilen, X3), X7),
                      
                      r3start=ifelse(strand=="+", pmax(X2, X3 - ilen), X6),
                      r3end=ifelse(strand=="+", X3, pmin(X6 + ilen, X7)),
                      
                      r4start=ifelse(strand=="+", X3, pmax(X5, X6 - elen)),
                      r4end=ifelse(strand=="+", pmin(X3 + elen, X4), X6),
                      
                      r5start=ifelse(strand=="+", pmax(X3, X4 - elen), X5),
                      r5end=ifelse(strand=="+", X4, pmin(X5 + elen, X6)),
                      
                      r6start=ifelse(strand=="+", X4, pmax(X4, X5 - ilen)),
                      r6end=ifelse(strand=="+", pmin(X4 + ilen, X5), X5),
                      
                      r7start=ifelse(strand=="+", pmax(X4, X5 - ilen), X4),
                      r7end=ifelse(strand=="+", X5, pmin(X4 + ilen, X5)),
                      
                      r8start=ifelse(strand=="+", X5, pmax(X3, X4 - elen)),
                      r8end=ifelse(strand=="+", pmin(X5 + elen, X6), X4),
                      
                      r9start=ifelse(strand=="+", pmax(X5, X6 - elen), X3),
                      r9end=ifelse(strand=="+", X6, pmin(X3 + elen, X4)),
                      
                      r10start=ifelse(strand=="+", X6, pmax(X2, X3 - ilen)),
                      r10end=ifelse(strand=="+", pmin(X6 + ilen, X7), X3),
                      
                      r11start=ifelse(strand=="+", pmax(X6, X7 - ilen), X2),
                      r11end=ifelse(strand=="+", X7, pmin(X2 + ilen, X3)),
                      
                      r12start=ifelse(strand=="+", X7, pmax(X1, X2 - elen)),
                      r12end=ifelse(strand=="+", pmin(X7 + elen, X8), X2))
  
  me <- me[,c("chr", "idMap", "strand",names(me)[grepl("r[0-9]{1,2}[a-z]{3,5}",names(me))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  tcga <- rbind(a3, a5, se, ri, me) %>% select(chr, start, end, idMap, regionType, strand)
  tcga <- tcga[(tcga$end - tcga$start) > 20,]
  
  message("GC calculate")
  options(scipen = 7)
  genomeFa <- "/home/u1357/RNAseq/refv19/GRCh37.p13.genome.fa"
  fileName <- paste0(outDir, "tcga_ase_rmaps_i", ilen, "_e", elen, "_region_separated.bed")
  write.table(tcga, file=fileName, sep="\t", row.names = F, col.names = F, quote = F)
  
  output <- gsub(".bed","_gc.txt",fileName)
  path <- "/pub/anaconda3/bin/"
  cmd <- paste0(path, "bedtools nuc -fi ", genomeFa," -bed ", fileName," -s >", output)
  system(cmd)
  
  message("all done\n")
}




