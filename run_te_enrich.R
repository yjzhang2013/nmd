##### TE enrichment #####
## 
tcga <- read.table("/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all.bed", sep="\t")
colnames(tcga) <- c("chr", "start", "end", "idMap", "spliceType", "strand")
outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/"
dir.create(outDir)
rmapsRegionNew <- function(tcga, ilen, elen, outDir){
  # ilen <- 250
  # elen <- 50
  message("tcga region process...\n")
  tcga <- tcga %>% mutate(geneName = limma::strsplit2(idMap, "\\|")[,1],
                          eventRegion = limma::strsplit2(idMap, "\\|")[,3])
  
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
  
  options(scipen = 7)
  
  fileName <- paste0(outDir, "tcga_ase_rmaps_i", ilen, "_e", elen, "_region_separated.bed")
  write.table(tcga, file=fileName, sep="\t", row.names = F, col.names = F, quote = F)

  message("all done\n")
}

rmapsRegionNew(tcga, 250, 50, outDir)


## splice site region
rmapsSite <- function(tcga, ilen, elen, outDir){
  # ilen <- 7
  # elen <- 7
  message("tcga region process...\n")
  tcga <- tcga %>% mutate(geneName = limma::strsplit2(idMap, "\\|")[,1],
                          eventRegion = limma::strsplit2(idMap, "\\|")[,3])
  
  tcga <- cbind(tcga, data.frame(limma::strsplit2(tcga$eventRegion, "[:, -]")))
  for (i in paste0("X",seq(8))){
    tcga[,i] <- as.numeric(tcga[,i])
  }
  
  
  a3 <- tcga[tcga$spliceType=="A3",]
  a3 <- a3 %>% mutate(s1start=ifelse(strand=="+", pmax(X1, X2 - elen), pmax(X4, X5 - ilen)),
                      s1end=ifelse(strand=="+", pmin(X2 + ilen, X3), pmin(X5 + elen, X6)),
                      
                      s2start=ifelse(strand=="+", pmax(X2, X3 - ilen), pmax(X1, X4 - elen)),
                      s2end=ifelse(strand=="+", pmin(X3 + elen, X6), pmin(X4 + ilen, X5)),
                      
                      s3start=ifelse(strand=="+", pmax(X2, X5 - ilen), pmax(X1, X2 - elen)),
                      s3end=ifelse(strand=="+", pmin(X5 + elen, X6), pmin(X2 + ilen, X5)))
  
  ## get s1~s6 (A3, A5, ES, IR, ME), start & end
  a3 <- a3[,c("chr", "idMap", "strand",names(a3)[grepl("s[1-9]{1,2}[a-z]{3,5}",names(a3))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  a5 <- tcga[tcga$spliceType=="A5",]
  a5 <- a5 %>% mutate(s1start=ifelse(strand=="-", pmax(X2, X5 - ilen), pmax(X1, X2 - elen)),
                      s1end=ifelse(strand=="-", pmin(X5 + elen, X6), pmin(X2 + ilen, X5)),
                      
                      s2start=ifelse(strand=="-", pmax(X2, X3 - ilen), pmax(X1, X4 - elen)),
                      s2end=ifelse(strand=="-", pmin(X3 + elen, X6), pmin(X4 + ilen, X5)),
                      
                      s3start=ifelse(strand=="-", pmax(X1, X2 - elen), pmax(X4, X5 - ilen)),
                      s3end=ifelse(strand=="-", pmin(X2 + ilen, X3), pmin(X5 + elen, X6)))
  
  a5 <- a5[,c("chr", "idMap", "strand",names(a5)[grepl("s[1-9]{1,2}[a-z]{3,5}",names(a5))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  se <- tcga[tcga$spliceType=="ES",]
  se <- se %>% mutate(s1start=ifelse(strand=="+", pmax(X1, X2 - elen), pmax(X4, X5 - ilen)),
                      s1end=ifelse(strand=="+", pmin(X2 + ilen, X3), pmin(X5 + elen, X6)),
                      
                      s2start=ifelse(strand=="+", pmax(X2, X3 - ilen), pmax(X3, X4 - elen)),
                      s2end=ifelse(strand=="+", pmin(X3 + elen, X4), pmin(X4 + ilen, X5)),
                      
                      s3start=ifelse(strand=="+", pmax(X3, X4 - elen), pmax(X2, X3 - ilen)),
                      s3end=ifelse(strand=="+", pmin(X4 + ilen, X5), pmin(X3 + elen, X4)),
                      
                      s4start=ifelse(strand=="+", pmax(X4, X5 - ilen), pmax(X1, X2 - elen)),
                      s4end=ifelse(strand=="+", pmin(X5 + elen, X6), pmin(X2 + ilen, X3)))
  
  se <- se[,c("chr", "idMap", "strand",names(se)[grepl("s[1-9]{1,2}[a-z]{3,5}",names(se))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  
  ri <- tcga[tcga$spliceType=="IR",]
  ri <- ri %>% mutate(s1start=ifelse(strand=="+", pmax(X1, X2 - elen), pmax(X2, X5 - ilen)),
                      s1end=ifelse(strand=="+", pmin(X2 + ilen, X5), pmin(X5 + elen, X6)),
                      
                      s2start=ifelse(strand=="+", pmax(X2, X5 - ilen), pmax(X1, X2 - elen)),
                      s2end=ifelse(strand=="+", pmin(X5 + elen, X6), pmin(X2 + ilen, X5)))
  
  ri <- ri[,c("chr", "idMap", "strand",names(ri)[grepl("s[1-9]{1,2}[a-z]{3,5}",names(ri))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  me <- tcga[tcga$spliceType=="ME",]
  me <- me %>% mutate(s1start=ifelse(strand=="+", pmax(X1, X2 - elen), pmax(X6, X7 - ilen)),
                      s1end=ifelse(strand=="+", pmin(X2 + ilen, X3), pmin(X7 + elen, X8)),
                      
                      s2start=ifelse(strand=="+", pmax(X2, X3 - ilen), pmax(X5, X6 - elen)),
                      s2end=ifelse(strand=="+", pmin(X3 + elen, X4), pmin(X6 + ilen, X7)),
                      
                      s3start=ifelse(strand=="+", pmax(X3, X4 - elen), pmax(X4, X5 - ilen)),
                      s3end=ifelse(strand=="+", pmin(X4 + ilen, X5), pmin(X5 + elen, X6)),
                      
                      s4start=ifelse(strand=="+", pmax(X4, X5 - ilen), pmax(X3, X4 - elen)),
                      s4end=ifelse(strand=="+", pmin(X5 + elen, X6), pmin(X4 + ilen, X5)),
                      
                      s5start=ifelse(strand=="+", pmax(X5, X6 - elen), pmax(X2, X3 - ilen)),
                      s5end=ifelse(strand=="+", pmin(X6 + ilen, X7), pmin(X3 + elen, X4)),
                      
                      s6start=ifelse(strand=="+", pmax(X6, X7 - ilen), pmax(X1, X2 - elen)),
                      s6end=ifelse(strand=="+", pmin(X7 + elen, X8), pmin(X2 + ilen, X3)))
  
  me <- me[,c("chr", "idMap", "strand",names(me)[grepl("s[1-9]{1,2}[a-z]{3,5}",names(me))])] %>% 
    tidyr::gather(regionType, loc, -c("chr", "idMap", "strand")) %>% 
    dplyr::mutate(type=ifelse(grepl("start$", regionType), "start", "end"),
                  regionType = gsub("start", "", gsub("end", "", regionType)),
                  idMap=paste0(idMap, "::",regionType)) %>% 
    tidyr::spread(type, loc) %>% 
    dplyr::select(chr, start, end, idMap, regionType, strand)
  
  
  tcga <- rbind(a3, a5, se, ri, me) %>% select(chr, start, end, idMap, regionType, strand)
  tcga <- tcga[(tcga$end - tcga$start) > ilen,]
  
  options(scipen = 7)
  fileName <- paste0(outDir, "tcga_ase_rmaps_i", ilen, "_e", elen, "_site_separated.bed")
  write.table(tcga, file=fileName, sep="\t", row.names = F, col.names = F, quote = F)
  
  message("all done\n")
}

rmapsSite(tcga, 5, 5, outDir)


## 
regionBed <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_separated.bed"
siteBed <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i5_e5_site_separated.bed"

teFile <- "/home/u1357/RNAseq/repeats/repeat_processed.bed"
regionFile1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/ase_repeat_overlaps_region.txt"
regionFile2 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/ase_repeat_overlaps_region_antiS.txt"
siteFile1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/ase_repeat_overlaps_site.txt"
siteFile2 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/ase_repeat_overlaps_site_antiS.txt"

cmd1 <- paste0("/pub/anaconda3/bin/bedtools intersect -a ", regionBed, " -b ", teFile, " -wo -s >", regionFile1)
cmd2 <- paste0("/pub/anaconda3/bin/bedtools intersect -a ", regionBed, " -b ", teFile, " -wo -S >", regionFile2)
cmd3 <- paste0("/pub/anaconda3/bin/bedtools intersect -a ", siteBed, " -b ", teFile, " -wo -s >", siteFile1)
cmd4 <- paste0("/pub/anaconda3/bin/bedtools intersect -a ", siteBed, " -b ", teFile, " -wo -S >", siteFile2)

system(cmd1)
system(cmd2)
system(cmd3)
system(cmd4)



## load result
library(dplyr)

regionSense <- read.table(regionFile1, sep="\t", header = F)[,-c(7,12)]
regionAntisense <- read.table(regionFile2, sep="\t", header = F)[,-c(7,12)]
name <- c("chr","aseStart0","aseEnd","idMap","regionType","strand","repStart0","repEnd","repID","repName","ovlapLen")
colnames(regionSense) <- colnames(regionAntisense) <- name

## data process
### region
regionSense <- regionSense %>% 
  dplyr::mutate(senseType="sense")

regionAntisense <- regionAntisense %>% 
  dplyr::mutate(senseType="antisense")

regionRes <- rbind(regionSense, regionAntisense)

regionRes <- regionRes %>% 
  dplyr::mutate(len = aseEnd-aseStart0,
                spliceType = limma::strsplit2(idMap, split="\\|")[,2])

regionRes <- regionRes %>% 
  dplyr::mutate(exonType=ifelse(spliceType=="A3", ifelse(regionType %in% c("r1","r4","r6"), "exon", "intron"),
                                ifelse(spliceType=="A5", ifelse(regionType %in% c("r1", "r3", "r6"), "exon", "intron"),
                                       ifelse(spliceType=="ES", ifelse(regionType %in% c("r1","r4","r5","r8"), "exon", "intron"),
                                              ifelse(spliceType=="RI", ifelse(regionType %in% c("r1","r4"), "exon", "intron"), 
                                                     ifelse(regionType %in% c("r1","r4","r5","r8","r9","r12"), "exon", "intron"))))),
                binNum = ifelse(exonType=="exon", 5, 25))


## resign map star/end similar to str_detect 
regionRes1 <- regionRes %>% 
  dplyr::filter(strand=="+") %>% 
  dplyr::mutate(mapStart=ifelse(aseStart0 >= repStart0, 0, repStart0-aseStart0), 
                mapEnd=mapStart+ovlapLen,
                mapStart=ifelse(mapStart==0, 1, mapStart))

regionRes2 <- regionRes %>% 
  dplyr::filter(strand=="-") %>% 
  dplyr::mutate(mapStart=ifelse(aseEnd <= repEnd, 0, aseEnd-repEnd), 
                mapEnd=mapStart+ovlapLen,
                mapStart=ifelse(mapStart==0, 1, mapStart))

regionRes <- rbind(regionRes1, regionRes2)




### region bin process
regionRes <- regionRes %>%
  dplyr::mutate(idMap = limma::strsplit2(idMap, split="::")[,1],
                binStart = ceiling(mapStart*binNum/len),
                binEnd = ceiling(mapEnd*binNum/len))

binLocFun <- function(name, dfls){
  df <- dfls[[name]]
  for (i in 1:nrow(df)){
    # i <- 1
    mt <- df[i,]
    binEnd <- mt$binEnd
    binStart <- mt$binStart
    n <- binEnd-binStart+1
    mt <- mt[rep(1,n),]
    mt$binLoc <- seq(binStart, binEnd)
    
    if(i==1){
      res <- mt
    }else{
      res <- rbind(res, mt)
    }
  }
  
  return(res)
}


source("/home/u1357/RNAseq/pancan/oncosplicingv3/script/defined_function.R")
dfls <- splitBy(regionRes, 40)
mt <- parallel::mclapply(seq(40), binLocFun, dfls, mc.cores = 40)
res <- do.call(rbind, mt)


regionRes <- res %>% 
  dplyr::mutate(binLoc = stringr::str_pad(binLoc, width = 2, side = "left",pad = 0),
                regionType = stringr::str_pad(gsub("r","",regionType), width = 2, side = "left",pad = 0),
                regionBinLoc = paste0("r",regionType, "_b", binLoc)) %>% 
  dplyr::select(-binLoc,-regionType)
saveRDS(regionRes, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge.rds")


regionRes <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge.rds")


###site
siteFile1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/ase_repeat_overlaps_site.txt"
siteFile2 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/ase_repeat_overlaps_site_antiS.txt"

siteSense <- read.table(siteFile1, sep="\t", header = F)[,-c(7,12)]
siteAntisense <- read.table(siteFile2, sep="\t", header = F)[,-c(7,12)]

name <- c("chr","aseStart0","aseEnd","idMap","regionType","strand","repStart0","repEnd","repID","repName","ovlapLen")
colnames(siteSense) <- colnames(siteAntisense) <- name


siteAntisense <- siteAntisense %>% 
  dplyr::mutate(senseType="antisense") %>% 
  dplyr::filter(ovlapLen>=8)

siteSense <- siteSense %>% 
  dplyr::mutate(senseType="sense") %>% 
  dplyr::filter(ovlapLen>=8)

siteRes <- rbind(siteAntisense, siteSense)
siteRes <- siteRes %>% 
  dplyr::mutate(spliceType=limma::strsplit2(idMap, split="\\|")[,2])

##
a3 <- siteRes[siteRes$spliceType=="A3",]
a3$regionType <- gsub("s1", "r2", gsub("s2", "r3", gsub("s3", "r5", a3$regionType)))

a5 <- siteRes[siteRes$spliceType=="A5",]
a5$regionType <- gsub("s1", "r2", gsub("s2", "r4", gsub("s3", "r5", a5$regionType)))

es <- siteRes[siteRes$spliceType=="ES",]
es$regionType <- gsub("s1", "r2", gsub("s2", "r3", gsub("s3", "r6", gsub("s4", "r7", es$regionType))))

ir <- siteRes[siteRes$spliceType=="IR",]
ir$regionType <- gsub("s1", "r2", gsub("s2", "r3", ir$regionType))

me <- siteRes[siteRes$spliceType=="ME",]
me$regionType <- gsub("s1", "r2", gsub("s2", "r3", gsub("s3", "r6", gsub("s4", "r7", gsub("s5", "r10", gsub("s6", "r11", me$regionType))))))

siteRes <- rbind(a3, a5, es, ir, me)

siteRes <- siteRes %>%
  dplyr::mutate(binLoc = ifelse(spliceType=="A3", ifelse(regionType %in% c("r2"), "0", "99"),
                                  ifelse(spliceType=="A5", ifelse(regionType %in% c("r2", "r4"), "0", "99"),
                                         ifelse(spliceType=="ES", ifelse(regionType %in% c("r2","r6"), "0", "99"),
                                                ifelse(spliceType=="RI", ifelse(regionType %in% c("r2"), "0", "99"), 
                                                       ifelse(regionType %in% c("r2","r6","r10"), "0", "99"))))),
                idMap = limma::strsplit2(idMap, split="::")[,1]) %>%
  dplyr::select(-spliceType)


siteRes <- siteRes %>% 
  dplyr::mutate(binLoc = stringr::str_pad(binLoc, width = 2, side = "left", pad = 0),
                regionType = stringr::str_pad(gsub("r","",regionType), width = 2, side = "left",pad = 0),
                regionBinLoc = paste0("r",regionType, ifelse(binLoc=="99", "_c", "_a"), binLoc)) %>% 
  dplyr::select(-binLoc,-regionType)
saveRDS(siteRes, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_iX_eX_region_separated_result_infoMerge.rds")


resMergeRE <- rbind(regionRes[,names(siteRes)], siteRes)
resMergeRE <- resMergeRE %>% 
  dplyr::mutate(spliceType=limma::strsplit2(idMap, split="\\|")[,2])

teFilter <- readRDS(file="/home/u1357/RNAseq/refv19/te_filter_top6class.rds")
resMergeRE <- merge(resMergeRE, 
                    teFilter %>% dplyr::select(repID, repClass, repFamily), 
                    by="repID")

saveRDS(resMergeRE, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_i5_e5_site_separated_result.rds")
resMergeRE <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_i5_e5_site_separated_result.rds")
length(unique(resMergeRE$idMap))



#### make check file ####
library(dplyr)
source("/home/u1357/RNAseq/pancan/oncosplicingv3/script/defined_function.R")

teFilter <- readRDS(file="/home/u1357/RNAseq/refv19/te_filter_top6class.rds")
teFilter <- teFilter %>% 
  dplyr::mutate(repName=paste(repClass, repFamily, repName, sep="::"),
                repFamily=paste(repClass, repFamily, sep="::"))

checkX <- function(repLevelCheck, level){
  # level <- "repAll"
  checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
  checkAll <- pasteAll(repLevelCheck, checks, sep=".")
  checkAll <- pasteAll(checkAll, c("sense","antisense"), sep=".")
  checkAll <- pasteAll(checkAll, c("control","input"), sep=".")
  checkAll <- data.frame(limma::strsplit2(checkAll, split="\\.") %>% as.data.frame(), freq=0, checks=checkAll)
  names(checkAll)[1:5] <- c("repLevel", "spliceType", "regionBinLoc", "senseType", "dataType")
  checkAll$repLevelSenseSpliceType <- paste(checkAll$repLevel, checkAll$senseType, checkAll$spliceType, sep=".")
  saveRDS(checkAll, file=paste0("/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/check_", level, ".rds"))
}

repLevelCheck <- c(unique(teFilter$repName), unique(teFilter$repFamily), unique(teFilter$repClass))
# checkX(unique(teFilter$repName), "repName")
# checkX(unique(teFilter$repFamily), "repFamily")
# checkX(unique(teFilter$repClass), "repClass")
checkX(repLevelCheck, "repAll")



#### all in summary ####
resMergeStatAllFun <- function(resMergeAll, type, outDir){
  # type <- "all"
  message("befor rmdup events, nrow:", nrow(resMergeAll))
  outDir <- paste0(outDir, "/", type, "/")
  dir.create(outDir)
  
  tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all.rds")
  
  if (type=="spliceseq"){
    tcga <- tcga %>% 
      dplyr::filter(!is.na(idType))
  }else if(type=="spladder"){
    tcga <- tcga %>% 
      dplyr::filter(!is.na(Splice_Event))
  }else if(type=="all"){
    tcga <- tcga %>% 
      dplyr::mutate(projects=ifelse(!is.na(idType) & !is.na(Splice_Event), 2, 
                                    ifelse(!is.na(idType), 1, 3))) %>% 
      arrange(projects)
  }

  dup2 <- gsub("(.*\\|.*\\|).*?-|-.+$", "\\1", tcga$idMap)
  dup1 <- gsub("(:).*(:)", ":", tcga$idMap)
  ase <- tcga$idMap[!duplicated(dup1) & !duplicated(dup2)]
  
  resMergeAll.filter <- resMergeAll %>% 
    dplyr::filter(idMap %in% ase)
  message("after rmdup events, nrow:", nrow(resMergeAll.filter))
  
  resMergeStat.repname <- resMergeAll.filter %>% 
    group_by(repName, senseType, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=n()) 
  
  resMergeStat.repname <- resMergeStat.repname %>% 
    split(resMergeStat.repname$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  name <- names(resMergeStat.repname)
  resMergeStat.repname <- names(resMergeStat.repname) %>% purrr::map(function(x) checkRegionBin(x, resMergeStat.repname, 3))
  names(resMergeStat.repname) <- name
  saveRDS(resMergeStat.repname, file=paste0(outDir, "/", "resMergeStat.repname.rds"))
  
  ##
  resMergeStat.repfamily <- resMergeAll.filter %>% 
    group_by(repFamily, senseType, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=n())
  
  resMergeStat.repfamily <- resMergeStat.repfamily %>% 
    split(resMergeStat.repfamily$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  name <- names(resMergeStat.repfamily)
  resMergeStat.repfamily <- names(resMergeStat.repfamily) %>% purrr::map(function(x) checkRegionBin(x, resMergeStat.repfamily, 3))
  names(resMergeStat.repfamily) <- name
  saveRDS(resMergeStat.repfamily, file=paste0(outDir, "/", "resMergeStat.repfamily.rds"))
  
  
  ##
  resMergeStat.repclass <- resMergeAll.filter %>% 
    group_by(repClass, senseType, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=n())
  
  resMergeStat.repclass <- resMergeStat.repclass %>% 
    split(resMergeStat.repclass$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  name <- names(resMergeStat.repclass)
  resMergeStat.repclass <- names(resMergeStat.repclass) %>% purrr::map(function(x) checkRegionBin(x, resMergeStat.repclass, 3))
  names(resMergeStat.repclass) <- name
  saveRDS(resMergeStat.repclass, file=paste0(outDir, "/", "resMergeStat.repclass.rds"))
  
  
  datals <- list(resMergeStat.repname, resMergeStat.repfamily, resMergeStat.repclass)
  names(datals) <- c("statByRepName", "statByRepFamily", "statByRepClass")
  # saveRDS(datals, file=paste0(outDir, "/", "resMergeStat.datals.rds"))
  # return(datals)
  
  source("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/teEnrich.R")
  lapply(names(datals), plotFreqMapStatType, datals, outDir)
}

outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/repeatSummary"
dir.create(outDir)
resMergeStatAllFun(resMergeRE, "all", outDir)

resMergeStatAllFun(resMergeRE, "spliceseq", outDir)

resMergeStatAllFun(resMergeRE, "spladder", outDir)




#### enrichment ####
resMergeAll <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_i5_e5_site_separated_result.rds")

enrichRBPase <- function(name, inputASEls, resMergeAll, dir, controlASE=NULL, type="spliceseq"){
  message(name)

  inputASE <- inputASEls[[name]]
  outDir <- paste0(dir, "/", type, "/", name, "/")
  dir.create(outDir, recursive = T)
  
  source("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/teEnrich.R")
  regionBinCheckPval(resMergeAll, inputASE, controlASE, outDir, type)
}


##
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all_txAnno_repeatAnno.rds")
tcga <- tcga %>%
  dplyr::filter(!is.na(idType))
table(tcga$txAnno)
NMD <- tcga$idMap[tcga$txAnno=="PCD_NMD"]
PCD <- tcga$idMap[tcga$txAnno=="PCD_PCD"]
PCDNA <- tcga$idMap[tcga$txAnno=="PCD_NA"]
NAPCD <- tcga$idMap[tcga$txAnno=="NA_PCD"]
OTHER <-tcga$idMap[tcga$txAnno=="Other"]
NANA <- tcga$idMap[tcga$txAnno=="notAnno"]
NMDnoRE <- tcga$idMap[tcga$txAnno=="PCD_NMD" & tcga$repeatAnno=="no"]
aseTypels <- list("NMD"=NMD, "PCD"=PCD, "PCDNA"=PCDNA, "NAPCD"=NAPCD, "NANA"=NANA, "NMDnoRE"=NMDnoRE,"ALUNMD"=ALUNMD, "ALUnoNMD"=ALUnoNMD)


ALUNMD <- tcga$idMap[tcga$txAnno=="PCD_NMD" & tcga$aluAnno=="yes" & tcga$spliceType=="ES"]
ALUnoNMD <- tcga$idMap[tcga$txAnno!="PCD_NMD" & tcga$aluAnno=="yes" & tcga$spliceType=="ES"]
aseTypels <- list("ALUNMD"=ALUNMD, "ALUnoNMD"=ALUnoNMD)


##
dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/aseInDiffTypesByRBPsControlOTHER/"
dir.create(dir, recursive = T)
parallel::mclapply(names(aseTypels), enrichRBPase, aseTypels, resMergeRE, dir, controlASE=OTHER, type="spliceseq", mc.cores = 5)


##
dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/aseInDiffTypesByRBPsControlNA/"
dir.create(dir, recursive = T)
parallel::mclapply(names(aseTypels), enrichRBPase, aseTypels, resMergeRE, dir, controlASE=NULL, type="spliceseq", mc.cores = 5)







#### integration ####
# res <- read.table("/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/spliceseq/20231026182651NMD/resMergeStat.pval.filter.txt", sep="\t", header = T)
# res <- read.table("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/aseInDiffTypesByRBPsControlOTHER/spliceseq/NMD/resMergeStat.pval.filter.txt", sep="\t", header = T)
res <- read.table("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/aseInDiffTypesByRBPsControlNA/spliceseq/NMD/resMergeStat.pval.filter.txt", sep="\t", header = T)

res <- res  %>% 
  dplyr::mutate(FDR=-log10(FDR)) %>% 
  dplyr::filter(grepl("(::).*(::)", repLevel)) %>% 
  dplyr::filter(spliceType=="ES") %>% 
  dplyr::filter(abs(FDR) > -log(0.05, 10)) %>% 
  dplyr::select(seq(4), FDR) %>% 
  tidyr::spread(regionBinLoc, FDR, fill=0)

checkRegionBinS <- function(df, name, preColN){
  # name <- "ES"
  # dfls <- resMergeStat.repname
  checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
  checks <- limma::strsplit2(checks, split="\\.")
  checkls <- split(checks[,2], checks[,1])
  
  check <- checkls[[name]]
  check <- check[order(check)]
  checkNO <- check[!check%in%names(df)]
  df[,checkNO] <- 0
  df <- df[,c(names(df)[1:preColN], check)]
  
  return(df)
}


result <- checkRegionBinS(res, "ES" ,3)
row.names(result) <- paste0(result$repLevel, "_", result$senseType)
result <- result[,-seq(3)]

# filter <- apply(data.frame(result>-log10(0.05)), 1, sum)
# result <- result[filter>=3,]


plotPvalMap2 <- function(df, title, fillColor){
  # df <- result
  # split <- limma::strsplit2(colnames(df), split = "_")[,1]
  split <- limma::strsplit2(colnames(df), split = "[0-9]{2}$")[,1]
  # col_fun2 = colorRamp2(c(-300, -10, -4, log(0.05, 10), 0, -log(0.05, 10), 4, 10, 300), c("green4", "green", "springgreen","white","white","white","Salmon","red","red4"))
  col_fun2 = colorRamp2(c(0, -log(0.05, 10), 4, 20, 40), c("white","white","Salmon","red","red4"))
  p <- Heatmap(df, col = col_fun2, cluster_rows=T,cluster_columns=F,
               row_title = NULL, column_split = split, column_gap = unit(0, "mm"), border = T,
               column_title = title, column_title_gp = gpar(fill = fillColor, fontsize=15, fontface="italic"),
               column_title_side = "bottom", 
               show_row_names = T,rect_gp = gpar(col = NA), show_column_dend = F,
               show_row_dend = F, show_column_names = T,
               heatmap_legend_param = list(title = "-log10(FDR)"))
  ComplexHeatmapdraw(p)
}


columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))
# pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/spliceseq/20231026182651NMD/sig_repname2.pdf"
# pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/aseInDiffTypesByRBPsControlOTHER/spliceseq/NMD/sig_repname2.pdf"
pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/aseInDiffTypesByRBPsControlNA/spliceseq/NMD/sig_repname2.pdf"
pdf(file=pdfName, width=ncol(result)*0.1, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()







#### coupling with RBP enriched ASEs √ #####
library(dplyr)
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdByRBPs/", pattern = "*pval.result.txt$", recursive = T, full.names = T)
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*/(.*[n,p]cor)/.*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)
for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>%
    dplyr::mutate(RBP=names(res)[i])
}
res <- do.call(rbind, res)

result2 <- res %>%
  dplyr::filter(dataType=="input") %>%
  dplyr::filter(spliceType=="ES") %>%
  dplyr::mutate(RBP=limma::strsplit2(RBP, split="_")[,1],
                dup2=paste0(RBP, idType)) %>%
  dplyr::filter(!duplicated(dup2)) %>%
  group_by(RBP) %>%
  # summarise(freq=n()) %>%
  dplyr::filter(n()>=10)
inputASEls <- split(result2$idMap, result2$RBP)
dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/"
dir.create(dir)

resMergeRE <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_i5_e5_site_separated_result.rds")
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all_txAnno_repeatAnno.rds")
nmd <- tcga %>%
  dplyr::filter(!is.na(idType)) %>%
  dplyr::filter(txAnno=="PCD_NMD")

teEnrichFun <- function(name, inputASEls, resMergeRE, nmd, dir){
  message(name)
  ## name <- names(inputASEls)[1]
  outDir <- paste0(dir, "/", name, "/")
  dir.create(outDir)

  inputASE <- inputASEls[[name]]
  controlASE <- nmd$idMap[!nmd$idMap%in%inputASE]

  source("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/teEnrich.R")
  regionBinCheckPval(resMergeRE, inputASE, controlASE, outDir, "spliceseq", minInputASEs=5)
}

parallel::mclapply(names(inputASEls), teEnrichFun, inputASEls, resMergeRE, nmd, dir, mc.cores = 20)



### RBP-NMD-Events enrich to repeats
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/", pattern = "*pval.filter.txt$", recursive = T, full.names = T)
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*enrich2rbpMotifEnrichedASEs//(.*)/.*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(RBP=names(res)[i])
}


res <- do.call(rbind, res)
result <- res %>% 
  dplyr::filter(spliceType=="ES") %>%
  # dplyr::filter(senseType=="antisense") %>% 
  dplyr::filter(grepl("(::).*(::)", repLevel)) %>% 
  dplyr::mutate(repLevel=paste0(RBP, "_", repLevel, "_", senseType)) %>% 
  # dplyr::filter(!duplicated(repLevel)) %>% 
  dplyr::select(repLevel, regionBinLoc, pValue) %>% 
  tidyr::spread(regionBinLoc, pValue, fill=0) %>% 
  as.data.frame()
row.names(result) <- result$repLevel
result <- result[,-1]

checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
checks <- limma::strsplit2(checks, split="\\.")
checks <- checks[,2][checks[,1]=="ES"]
checks <- checks[order(checks)]
checkNO <- checks[!checks%in%names(result)]
result[,checkNO] <- 0
result <- result[,checks]

filter <- apply(data.frame(result>-log10(0.05)), 1, sum)
result <- result[filter>=5,]


columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/pan2_rbp.pdf"
pdf(file=pdfName, width=ncol(result)*0.2, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()


##2
filter <- apply(data.frame(result>-log10(0.001)), 1, sum)
result <- result[filter>=5,]

columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/pan2_rbp2.pdf"
pdf(file=pdfName, width=ncol(result)*0.2, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()



## RBP-to-repName (Alu family)
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/", pattern = "*pval.filter.txt$", recursive = T, full.names = T)
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*enrich2rbpMotifEnrichedASEs//(.*)/.*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(RBP=names(res)[i])
}
res <- do.call(rbind, res)
result.filter <- res %>% 
  dplyr::filter(spliceType=="ES") %>%
  # dplyr::filter(senseType=="antisense") %>% 
  dplyr::filter(grepl("(::).*(::)", repLevel)) %>% 
  dplyr::filter(pValue> -log10(0.001)) %>% 
  dplyr::mutate(term=paste(term, RBP, sep=".")) %>%  ## new 2023.12.19
  group_by(term) %>%
  summarise(pValue=max(pValue))
# result.filter <- unique(result.filter$term)


##
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/", pattern = "*pval.result.txt$", recursive = T, full.names = T)
# files <- files[grepl("_ncor", files)]
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*enrich2rbpMotifEnrichedASEs//(.*)/.*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(RBP=names(res)[i])
}

res <- do.call(rbind, res)
result <- res %>% 
  # dplyr::mutate(term=paste(repClassFamilyName, senseType, spliceType, regionBinLoc, sep=".")) %>% 
  dplyr::mutate(term=paste(repClassFamilyName, senseType, spliceType, regionBinLoc, RBP, sep=".")) ## new 2023.12.19
  
result <- merge(result, result.filter[,c("term", "pValue")], by="term")



result4sanky <- result %>% 
  dplyr::filter(spliceType=="ES" & dataType=="input") %>% 
  # dplyr::filter(senseType=="antisense") %>% 
  dplyr::mutate(dup=paste0(RBP, repName, idMap, senseType),
                rbpRepSense=paste(RBP, repName, senseType, sep=".")) %>% 
  dplyr::filter(!duplicated(dup)) %>%
  group_by(RBP, repName, senseType) %>% 
  summarise(freq=n(), pval=max(pValue)) %>% 
  dplyr::filter(freq>5) %>% 
  dplyr::filter(pval>0)


length(unique(result$RBP))
length(unique(result$repName))

## color:  Salmon, #D5A9E3
library(ggalluvial)
ggplot(data = result4sanky, aes(axis1 = RBP, axis2 = repName, y = freq)) +
  scale_x_discrete(limits = c("RBP","repName"), expand = c(.01, .05)) +
  scale_y_continuous(labels = scales::label_number())+
  ggalluvial::geom_alluvium(aes(color=pval, fill=pval)) +
  scale_fill_continuous(low="lightblue",high="#D5A9E3")+
  scale_color_continuous(low="lightblue",high="#D5A9E3")+
  geom_stratum() + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("AS-transcript switch") +
  theme(axis.ticks = element_line(color = "black"),
        panel.grid.major = element_line(color="white",size=0.2,linetype = "dotted"),
        panel.grid.minor = element_line(color="white",size=0.2,linetype = "dotted"))
fileName1 <- paste0("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/ase_te_switch_sanky4.pdf")
ggsave(file=fileName1,width=7,height=13)




## RBP-motif on structured Alu-ASEs
res <- read.table("/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/aseInDiffTypesByRBPsControlNULL/NMD/resMergeStat.pval.result.txt", sep="\t", header = T)

result2 <- res %>% 
  dplyr::filter(dataType=="input" & spliceType=="ES") %>% 
  dplyr::mutate(dup2=paste0(RBP, idType)) %>% 
  dplyr::filter(!duplicated(dup2)) %>% 
  group_by(RBP) %>%
  # summarise(freq=n()) %>%
  dplyr::filter(n()>=20)

result2 <- result2 %>% 
  dplyr::filter(idMap %in% res.motif.filter$idMap) %>% 
  dplyr::mutate(RBP=limma::strsplit2(term, split = "\\.")[,1],
                motif=limma::strsplit2(term, split = "\\.")[,2],
                rbpIdMapBinLoc=paste0(RBP, idMap, regionBinLoc))



##
res.motif <- result %>% 
  dplyr::filter(spliceType=="ES" & dataType=="input") %>% 
  # dplyr::filter(senseType=="antisense") %>% 
  dplyr::mutate(dup=paste0(RBP, repName, idMap, senseType),
                rbpRepSense=paste(RBP, repName, senseType, sep=".")) %>% 
  dplyr::filter(!duplicated(dup)) %>%
  group_by(RBP, repName, senseType) %>% 
  dplyr::filter(n()>5 & max(pValue)>0)

res.motif.filter <- result %>% 
  dplyr::filter(spliceType=="ES" & dataType=="input") %>% 
  dplyr::mutate(dup=paste0(RBP, repName, idMap, senseType),
                rbpRepSense=paste(RBP, repName, senseType, sep="."),
                rbpIdMapBinLoc=paste0(RBP, idMap, regionBinLoc)) %>% 
  dplyr::filter(rbpRepSense %in% res.motif$rbpRepSense)


## final result (motif hit on structured ASE)
res.motif.final <- merge(result2, res.motif.filter[,c("rbpIdMapBinLoc","repName", "pValue")],by="rbpIdMapBinLoc")
csv <- res.motif.final %>% 
  dplyr::select(idMap, idType, chr, strand, spliceType, regionBinLoc, RBP, motif, repName, pValue) %>% 
  arrange(repName, RBP)
write.csv(csv, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/ase_te_motif.csv")


res.motif.final <- res.motif.final %>% 
  dplyr::mutate(motifRep=paste0(repName, "_", motif),
                dup=paste0(motifRep, idType)) %>%
  group_by(motifRep, regionBinLoc) %>%
  summarise(pValue=max(pValue)) %>% 
  tidyr::spread(regionBinLoc, pValue,fill=0) %>% 
  as.data.frame()
row.names(res.motif.final) <- res.motif.final$motifRep
res.motif.final <- res.motif.final[,-1]

checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
checks <- limma::strsplit2(checks, split="\\.")
checks <- checks[,2][checks[,1]=="ES"]
checks <- checks[order(checks)]
checkNO <- checks[!checks%in%names(res.motif.final)]
res.motif.final[,checkNO] <- 0
res.motif.final <- res.motif.final[,checks]

columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifEnrichedASEs/ase_te_motif.pdf"
pdf(file=pdfName, width=ncol(res.motif.final)*0.2, height=nrow(res.motif.final)*0.2+1.5)
plotMotif(res.motif.final, columnTitle, fillColor)
dev.off()


plotMotif <- function(df, title, fillColor){
  # df <- result
  # split <- limma::strsplit2(colnames(df), split = "_")[,1]
  split <- limma::strsplit2(colnames(df), split = "[0-9]{2}$")[,1]
  # col_fun2 = colorRamp2(c(-300, -10, -4, log(0.05, 10), 0, -log(0.05, 10), 4, 10, 300), c("green4", "green", "springgreen","white","white","white","Salmon","red","red4"))
  col_fun2 = colorRamp2(c(0, -log(0.05, 10), 4, 20, 40), c("white","white","Salmon","red","red4"))
  p <- Heatmap(df, col = col_fun2, cluster_rows=F,cluster_columns=F,
               row_title = NULL, column_split = split, column_gap = unit(0, "mm"), border = T,
               column_title = title, column_title_gp = gpar(fill = fillColor, fontsize=15, fontface="italic"),
               column_title_side = "bottom", 
               show_row_names = T,rect_gp = gpar(col = NA), show_column_dend = F,
               show_row_dend = F, show_column_names = T,
               heatmap_legend_param = list(title = "-log10(FDR)"))
  draw(p)
}



#### coupling with RBP correlated ASEs √ ####
library(dplyr)
corr.ss.filter <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/correlation/mapas_spliceseq_correlation_filtered.rds")
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all_txAnno_repeatAnno.rds")
tcga <- tcga %>%
  dplyr::filter(!is.na(idType))
nmd <- tcga$idMap[tcga$txAnno=="PCD_NMD"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.NMD <- corr.ss.filter %>%
  dplyr::filter(Map_Event %in% nmd & RBP %in% rbpWithMotif) %>%
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>%
  group_by(RBP, Splice_Type, Corr_Type) %>%
  dplyr::filter(n()>50)
length(unique(corr.ss.NMD$RBP))

inputASEls <- split(corr.ss.NMD$Map_Event, corr.ss.NMD$rbpCorrType)
dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifASEs/"
dir.create(dir)

resMergeRE <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/tcga_ase_rmaps_i250_e50_region_i5_e5_site_separated_result.rds")

teEnrichFun <- function(name, inputASEls, resMergeRE, dir){
  message(name)
  ## name <- names(inputASEls)[1]
  outDir <- paste0(dir, "/", name, "/")
  dir.create(outDir)
  
  inputASE <- inputASEls[[name]]
  source("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/teEnrich.R")
  regionBinCheckPval(resMergeRE, inputASE, controlASE=NULL, outDir, "spliceseq")
}

parallel::mclapply(names(inputASEls), teEnrichFun, inputASEls, resMergeRE, dir, mc.cores = 20)




### RBP-NMD-Events enrich to repeats
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifASEs/", pattern = "*pval.filter.txt$", recursive = T, full.names = T)
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*enrich2rbpMotifASEs//(.*)/.*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(RBP=names(res)[i])
}


res <- do.call(rbind, res)
result <- res %>% 
  dplyr::filter(spliceType=="ES" & senseType=="antisense") %>% 
  dplyr::filter(grepl("(::).*(::)", repLevel)) %>% 
  dplyr::mutate(repLevel=paste0(RBP, "_", repLevel)) %>% 
  dplyr::select(repLevel, regionBinLoc, pValue) %>% 
  tidyr::spread(regionBinLoc, pValue, fill=0) %>% 
  as.data.frame()
row.names(result) <- result$repLevel
result <- result[,-1]

checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
checks <- limma::strsplit2(checks, split="\\.")
checks <- checks[,2][checks[,1]=="ES"]
checks <- checks[order(checks)]
checkNO <- checks[!checks%in%names(result)]
result[,checkNO] <- 0
result <- result[,checks]

filter <- apply(data.frame(result>-log10(0.05)), 1, sum)
result <- result[filter>=5,]


columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifASEs/pan_rbp.pdf"
pdf(file=pdfName, width=ncol(result)*0.2, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()


##2
filter <- apply(data.frame(result>-log10(0.001)), 1, sum)
result <- result[filter>=5,]

columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

pdfName <- "/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifASEs/pan_rbp2.pdf"
pdf(file=pdfName, width=ncol(result)*0.2, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()




## RBP-to-repName ()
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifASEs/", pattern = "*pval.filter.txt$", recursive = T, full.names = T)
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*enrich2rbpMotifASEs//(.*)/.*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(RBP=names(res)[i])
}
res <- do.call(rbind, res)
result.filter <- res %>% 
  dplyr::filter(spliceType=="ES" & senseType=="antisense") %>% 
  dplyr::filter(grepl("(::).*(::)", repLevel)) %>% 
  dplyr::filter(pValue> -log10(0.001)) %>% 
  group_by(term) %>% 
  summarise(pValue=max(pValue))


##
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifASEs/", pattern = "*pval.result.txt$", recursive = T, full.names = T)
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*enrich2rbpMotifASEs//(.*)/.*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(RBP=names(res)[i])
}

res <- do.call(rbind, res)
result <- res %>% 
  dplyr::mutate(term=paste(repClassFamilyName, senseType, spliceType, regionBinLoc, sep="."))
result <- merge(result, result.filter[,c("term", "pValue")], by="term")


result <- result %>% 
  dplyr::filter(spliceType=="ES" & senseType=="antisense" & dataType=="input") %>% 
  dplyr::mutate(RBP=limma::strsplit2(RBP, split="_")[,1],
                dup=paste0(RBP, repName, idMap)) %>% 
  dplyr::filter(!duplicated(dup)) %>%
  group_by(RBP, repName) %>% 
  summarise(freq=n(), pval=max(pValue)) %>% 
  dplyr::filter(freq>10) %>% 
  dplyr::filter(pval>0)


length(unique(result$RBP))
length(unique(result$repName))

## color:  Salmon, #D5A9E3
library(ggalluvial)
ggplot(data = result, aes(axis1 = RBP, axis2 = repName, y = freq)) +
  scale_x_discrete(limits = c("RBP","repName"), expand = c(.01, .05)) +
  ggalluvial::geom_alluvium(aes(color=pval, fill=pval)) +
  scale_fill_continuous(low="lightblue",high="#D5A9E3")+
  scale_color_continuous(low="lightblue",high="#D5A9E3")+
  geom_stratum() + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("AS-transcript switch") +
  theme(axis.ticks = element_line(color = "black"),
        panel.grid.major = element_line(color="white",size=0.2,linetype = "dotted"),
        panel.grid.minor = element_line(color="white",size=0.2,linetype = "dotted"))
fileName1 <- paste0("/home/u1357/RNAseq/pancan/oncosplicingv3/teEnrich/results/enrich2rbpMotifASEs/ase_te_switch_sanky4.pdf")
ggsave(file=fileName1,width=7,height=13)
