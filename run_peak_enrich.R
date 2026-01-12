#### peak enrichment ####

#### 1. ase intersect peak ####
peakDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/webdata/MapAS/peaktemp/"
files <- list.files(peakDir, pattern = ".*rep_merge.bed", full.names = T)

bed <- lapply(files, data.table::fread, header=F)
names(bed) <- limma::strsplit2(basename(files), split = "_")[,1]
bed <- do.call(rbind, bed)
bed <- bed %>% 
  dplyr::mutate(V4=limma::strsplit2(V4, split="_")[,1],
                V4=paste(V4, V1, V2, V3, sep="::"))
peakFile <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/pan_gene_peak_merge.bed"
write.table(bed, file=peakFile, sep = "\t", col.names = F, row.names = F, quote = F)


##
tcga <- read.table("/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all.bed", sep="\t")
colnames(tcga) <- c("chr", "start", "end", "idMap", "spliceType", "strand")
outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/"
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
regionBed <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_i250_e50_region_separated.bed"
siteBed <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_i5_e5_site_separated.bed"
peakBed <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/pan_gene_peak_merge.bed"

regionFile1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/ase_peak_overlaps_region.txt"
siteFile1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/ase_peak_overlaps_site.txt"

cmd1 <- paste0("/pub/anaconda3/bin/bedtools intersect -a ", regionBed, " -b ", peakBed, " -wo -s >", regionFile1)
cmd3 <- paste0("/pub/anaconda3/bin/bedtools intersect -a ", siteBed, " -b ", peakBed, " -wo -s >", siteFile1)

system(cmd1)
system(cmd3)




#### 2. data process ####
library(dplyr)
regionFile1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/ase_peak_overlaps_region.txt"
regionRes <- read.table(regionFile1, sep="\t", header = F)[,-c(7, 11, 12)]
name <- c("chr","aseStart0","aseEnd","idMap","regionType","strand","peakStart0","peakEnd","peakID","ovlapLen")
colnames(regionRes) <- name

## data process
### region
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
  dplyr::mutate(mapStart=ifelse(aseStart0 >= peakStart0, 0, peakStart0-aseStart0), 
                mapEnd=mapStart+ovlapLen,
                mapStart=ifelse(mapStart==0, 1, mapStart))

regionRes2 <- regionRes %>% 
  dplyr::filter(strand=="-") %>% 
  dplyr::mutate(mapStart=ifelse(aseEnd <= peakEnd, 0, aseEnd-peakEnd), 
                mapEnd=mapStart+ovlapLen,
                mapStart=ifelse(mapStart==0, 1, mapStart))

regionRes <- rbind(regionRes1, regionRes2)


### region bin process
regionRes <- regionRes %>%
  dplyr::mutate(idMap = limma::strsplit2(idMap, split="::")[,1],
                binStart = ceiling(mapStart*binNum/len),
                binEnd = ceiling(mapEnd*binNum/len))


source("/home/u1357/RNAseq/pancan/oncosplicingv3/script/defined_function.R")
dfls <- splitBy(regionRes, 300)

dir.create("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/temp/")
saveSplitRDS <- function(i, dfls){
  saveRDS(dfls[[i]], file = paste0("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/temp/regionRes_", i,".rds"))
}
lapply(seq(300), saveSplitRDS, dfls)

rm(regionRes, regionRes1, regionRes2, dfls)
gc()



library(dplyr)
binLocFun <- function(name){
  # name <- 1
  message(name)
  df <- readRDS(paste0("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/temp/regionRes_", name,".rds"))
  # df <- df[1:10000,]
  
  res <- lapply(1:nrow(df), function(i, df){
    # i <- 3
    mt <- df[i,]
    binEnd <- mt$binEnd
    binStart <- mt$binStart
    n <- binEnd-binStart+1
    mt <- mt[rep(1,n),]
    mt$binLoc <- seq(binStart, binEnd)
    return(mt)
  }, df)
  
  res <- do.call(rbind, res)
  return(res)
}


mt <- parallel::mclapply(seq(300), binLocFun, mc.cores = 30)
# mt <- lapply(seq(40), binLocFun, dfls)
regionRes <- do.call(rbind, mt)

regionRes <- regionRes %>% 
  dplyr::mutate(binLoc = stringr::str_pad(binLoc, width = 2, side = "left",pad = 0),
                regionType = stringr::str_pad(gsub("r","",regionType), width = 2, side = "left",pad = 0),
                regionBinLoc = paste0("r",regionType, "_b", binLoc)) %>% 
  dplyr::select(-binLoc,-regionType)
saveRDS(regionRes, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge.rds")



###site
siteFile1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/ase_peak_overlaps_site.txt"

siteRes <- read.table(siteFile1, sep="\t", header = F)[,-c(7,11,12)]
name <- c("chr","aseStart0","aseEnd","idMap","regionType","strand","peakStart0","peakEnd","peakID","ovlapLen")
colnames(siteRes) <- name

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


##
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
saveRDS(siteRes, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_i5_e5_site_separated_result_infoMerge.rds")


##
library(dplyr)
resMergePeak <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge.rds")
siteRes <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_i5_e5_site_separated_result_infoMerge.rds")
resMergePeak <- rbind(resMergePeak[,names(siteRes)], siteRes)

##
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all.rds")
ase.ss <- tcga[!is.na(tcga$idType), c("idMap", "idType")]
ase.sa <- tcga[!is.na(tcga$Splice_Event), c("idMap", "Splice_Event")] %>% 
  dplyr::rename(idType=Splice_Event)

ase.encd <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/asetotalhk_uniqueASE.rds")
ase.encd <- ase.encd[ase.encd%in%tcga$idMap]


##
resMerge.encd <- resMergePeak %>% 
  dplyr::select(idMap, regionBinLoc, peakID, chr, strand) %>% 
  dplyr::filter(idMap %in% ase.encd) %>% 
  dplyr::mutate(spliceType=limma::strsplit2(idMap, split="\\|")[,2],
                RBP=limma::strsplit2(peakID, split="::")[,1])

saveRDS(resMerge.encd, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_iX_eX_region_site_separated_result_infoMerge_encode.rds")
rm(resMerge.encd)
gc()


##
resMerge.ss <- resMergePeak %>% 
  dplyr::select(idMap, regionBinLoc, peakID, chr, strand) %>% 
  dplyr::filter(idMap %in% ase.ss$idMap) %>% 
  dplyr::mutate(spliceType=limma::strsplit2(idMap, split="\\|")[,2],
                RBP=limma::strsplit2(peakID, split="::")[,1])

resMerge.ss <- merge(ase.ss, resMerge.ss, by="idMap")
saveRDS(resMerge.ss, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_iX_eX_region_site_separated_result_infoMerge_spliceseq.rds")
rm(resMerge.ss)
gc()


##
resMerge.sa <- resMergePeak %>% 
  dplyr::select(idMap, regionBinLoc, peakID, chr, strand) %>% 
  dplyr::filter(idMap %in% ase.sa$idMap) %>% 
  dplyr::mutate(spliceType=limma::strsplit2(idMap, split="\\|")[,2],
                RBP=limma::strsplit2(peakID, split="::")[,1])

resMerge.sa <- merge(ase.sa, resMerge.sa, by="idMap")
saveRDS(resMerge.sa, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_iX_eX_region_site_separated_result_infoMerge_spladder.rds")
rm(resMerge.sa)
gc()





#### 3. enrichment ####
library(dplyr)
corr.ss.filter <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/correlation/mapas_spliceseq_correlation_filtered.rds")
table(corr.ss.filter$RBP, corr.ss.filter$Corr_Type)
ase.stat <- corr.ss.filter %>% 
  group_by(Map_Event, Corr_Type) %>% 
  summarise(freq=n()) %>% 
  tidyr::spread(Corr_Type, freq, fill=0) %>% 
  dplyr::mutate(SpliceType=limma::strsplit2(Map_Event, "\\|")[,2]) %>% 
  as.data.frame()

library(dplyr)
resMerge.ss <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_iX_eX_region_site_separated_result_infoMerge_spliceseq.rds")

resMerge.ss <- resMerge.ss %>% 
  dplyr::mutate(motifId=RBP)


##### 3.1. enrichment for RBP-correlated ASEs ##### 
rbpWithPeak <- unique(resMerge.ss$RBP)
rbpWithPeak <- rbpWithPeak[order(rbpWithPeak)]

corr.ss.filter.ls1 <- corr.ss.filter %>% 
  dplyr::filter(RBP %in% rbpWithPeak) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>20)

corr.ss.filter.ls2 <- corr.ss.filter %>% 
  dplyr::filter(RBP %in% rbpWithPeak) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_all")) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>20)

corr.ss.filter.ls <- rbind(corr.ss.filter.ls1, corr.ss.filter.ls2)
corr.ss.filter.ls <- split(corr.ss.filter.ls$Map_Event, corr.ss.filter.ls$rbpCorrType)


enrichRBPaseMotif <- function(name, inputASEls, resMergeAll, dir, rbp="", controlASE=NULL, binMotifUnique=F){
  message(name)
  ## inputASEls <- inputRBPasels
  ## name <- names(inputASEls)[1]
  # name <- "GRWD1_ncor"
  
  
  inputASE <- inputASEls[[name]]
  if (!is.null(rbp)){
    rbp <- limma::strsplit2(name, split="_")[,1]
  }
  outDir <- paste0(dir, "/", rbp, "/", name, "/")
  dir.create(outDir, recursive = T)
  
  source("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/peakEnrich.R")
  regionBinCheckPval(resMergeAll, inputASE, controlASE, rbp, outDir,binMotifUnique)
}


## ALL RBP-CORRTYPE
dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/aseByRBPs/"
dir.create(dir, recursive = T)

inputRBPasels <-  corr.ss.filter.ls
parallel::mclapply(names(inputRBPasels), enrichRBPaseMotif, inputRBPasels, resMerge.ss, 
                   dir=dir, rbp="", controlASE=NULL,
                   mc.cores = 20)

# resMergeAll.bak <- resMergeAll
# resMergeAll <- resMerge.ss

##### 3.2 enrichment for RBP-correlated NMD-ASEs ##### 
## controlASE should be non-correlated
corr.ss.filter <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/correlation/mapas_spliceseq_correlation_filtered.rds")
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all_txAnno_repeatAnno.rds")
library(dplyr)
tcga <- tcga %>%
  dplyr::filter(!is.na(idType))


##### 3.2.1 NMD  #####
nmd <- tcga$idMap[tcga$txAnno=="PCD_NMD"]

rbpWithPeak <- unique(resMerge.ss$RBP)
rbpWithPeak <- rbpWithPeak[order(rbpWithPeak)]

corr.ss.NMD <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% nmd & RBP %in% rbpWithPeak) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.NMD$RBP))

inputASEls <- split(corr.ss.NMD$Map_Event, corr.ss.NMD$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/nmdByRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMerge.ss, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)



## pval & plot
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/nmdByRBPs/", pattern = "*pval.filter.txt$", recursive = T, full.names = T)
# files <- files[grepl("_ncor", files)]
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*([n,p]cor).*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(corrType=names(res)[i])
}
res <- do.call(rbind, res)

termFilter <- unique(res$term[res$pValue> -log10(0.05)])
result <- res %>% 
  # dplyr::filter(spliceType=="ES") %>% 
  dplyr::filter(spliceType=="ES" & pValue> -log10(0.05)) %>% 
  dplyr::mutate(pValue=ifelse(corrType=="pcor", -pValue, pValue),
                dup=paste0(RBPmotifId, regionBinLoc)) %>% 
  dplyr::filter(!duplicated(dup)) %>% 
  dplyr::select(RBPmotifId, regionBinLoc, pValue) %>% 
  tidyr::spread(regionBinLoc, pValue, fill=0) %>% 
  as.data.frame()
row.names(result) <- result$RBPmotifId
result <- result[,-1]

checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
checks <- limma::strsplit2(checks, split="\\.")
checks <- checks[,2][checks[,1]=="ES"]
checks <- checks[order(checks)]
checkNO <- checks[!checks%in%names(result)]
result[,checkNO] <- 0
result <- result[,checks]


columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

# pdfFile <- "/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/motif/rbps/pan_rbp.pdf"
pdfFile <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/nmdByRBPs/pan_rbp.pdf"
pdf(file=pdfFile, width=ncol(result)*0.15, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()




##### 3.2.2 alu #####
alu <- tcga$idMap[tcga$txAnno=="PCD_NMD" & tcga$aluAnno=="yes"]

rbpWithPeak <- unique(resMerge.ss$RBP)
rbpWithPeak <- rbpWithPeak[order(rbpWithPeak)]

corr.ss.alu <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% alu & RBP %in% rbpWithPeak) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.alu$RBP))

inputASEls <- split(corr.ss.alu$Map_Event, corr.ss.alu$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/nmdAluByRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMerge.ss, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)


## pval & plot
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/nmdAluByRBPs/", pattern = "*pval.filter.txt$", recursive = T, full.names = T)
# files <- files[grepl("_ncor", files)]
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*([n,p]cor).*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(corrType=names(res)[i])
}
res <- do.call(rbind, res)

termFilter <- unique(res$term[res$pValue> -log10(0.05)])
# termFilter <- unique(res$term[res$pValue < log10(0.001)])

result <- res %>% 
  # dplyr::filter(spliceType=="ES") %>% 
  dplyr::filter(spliceType=="ES" & pValue> -log10(0.05)) %>% 
  dplyr::mutate(pValue=ifelse(corrType=="pcor", -pValue, pValue),
                dup=paste0(RBPmotifId, regionBinLoc)) %>% 
  dplyr::filter(!duplicated(dup)) %>% 
  dplyr::select(RBPmotifId, regionBinLoc, pValue) %>% 
  tidyr::spread(regionBinLoc, pValue, fill=0) %>% 
  as.data.frame()
row.names(result) <- result$RBPmotifId
result <- result[,-1]

checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
checks <- limma::strsplit2(checks, split="\\.")
checks <- checks[,2][checks[,1]=="ES"]
checks <- checks[order(checks)]
checkNO <- checks[!checks%in%names(result)]
result[,checkNO] <- 0
result <- result[,checks]


columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

pdfFile <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/nmdAluByRBPs/pan_rbp.pdf"
pdf(file=pdfFile, width=ncol(result)*0.15, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()




##### 3.2.3 aluNoNMD #####
aluNoNMD <- tcga$idMap[tcga$txAnno!="PCD_NMD" & tcga$aluAnno=="yes"]

rbpWithPeak <- unique(resMerge.ss$RBP)
rbpWithPeak <- rbpWithPeak[order(rbpWithPeak)]

corr.ss.aluNoNMD <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% aluNoNMD & RBP %in% rbpWithPeak) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.aluNoNMD$RBP))

inputASEls <- split(corr.ss.aluNoNMD$Map_Event, corr.ss.aluNoNMD$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/AluNoNMDbyRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMerge.ss, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)


## pval & plot
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/AluNoNMDbyRBPs/", pattern = "*pval.filter.txt$", recursive = T, full.names = T)
# files <- files[grepl("_ncor", files)]
res <- lapply(files, read.table, sep="\t", header=T)
names(res) <- gsub(".*([n,p]cor).*","\\1", files)
res <- res %>% purrr::keep(function(x) nrow(x)>0)

for (i in 1:length(res)){
  res[[i]] <- res[[i]] %>% 
    dplyr::mutate(corrType=names(res)[i])
}


res <- do.call(rbind, res)

termFilter <- unique(res$term[res$pValue> -log10(0.05)])
# termFilter <- unique(res$term[res$pValue < log10(0.001)])

result <- res %>% 
  # dplyr::filter(spliceType=="ES") %>% 
  dplyr::filter(spliceType=="ES" & pValue> -log10(0.05)) %>% 
  dplyr::mutate(pValue=ifelse(corrType=="pcor", -pValue, pValue),
                dup=paste0(RBPmotifId, regionBinLoc)) %>% 
  dplyr::filter(!duplicated(dup)) %>% 
  dplyr::select(RBPmotifId, regionBinLoc, pValue) %>% 
  tidyr::spread(regionBinLoc, pValue, fill=0) %>% 
  as.data.frame()
row.names(result) <- result$RBPmotifId
result <- result[,-1]

checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
checks <- limma::strsplit2(checks, split="\\.")
checks <- checks[,2][checks[,1]=="ES"]
checks <- checks[order(checks)]
checkNO <- checks[!checks%in%names(result)]
result[,checkNO] <- 0
result <- result[,checks]


columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))

pdfFile <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/AluNoNMDbyRBPs/pan_rbp.pdf"
pdf(file=pdfFile, width=ncol(result)*0.15, height=nrow(result)*0.2+1.5)
plotPvalMap2(result, columnTitle, fillColor)
dev.off()





#### STAT ####
##### stat by ase groups ##### 
prefix <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_iX_eX_region_site_separated_result_infoMerge_"
resMergeSummaryFun <- function(resMergeAll, prefix ,type){
  # type <- "spliceseq"
  suppressMessages(library(dplyr))
  message("befor rmdup events, nrow:", nrow(resMergeAll))
  tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all.rds")
  ase <- unique(tcga$idMap)
  dup1 <- gsub("(:).*(:)", ":", ase)
  dup2 <- gsub("(.*\\|.*\\|).*?-|-.+$", "\\1", ase)
  # ase <- ase[!duplicated(dup1)]
  ase <- ase[!duplicated(dup1) & !duplicated(dup2)]
  resMergeAll <- resMergeAll %>% 
    dplyr::filter(idMap %in% ase)
  message("after rmdup events, nrow:", nrow(resMergeAll))
  
  resMergeStat.motifid <- resMergeAll %>% 
    group_by(motifId, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=n()) 
  
  resMergeStat.motifid <- resMergeStat.motifid %>% 
    split(resMergeStat.motifid$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  fileName <- paste0(prefix, type, "_statByMotifid.rds")
  saveRDS(resMergeStat.motifid, file = fileName)
  
  
  ##
  resMergeStat.rbp <- resMergeAll %>% 
    group_by(RBP, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=n())
  
  resMergeStat.rbp <- resMergeStat.rbp %>% 
    split(resMergeStat.rbp$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  fileName <- paste0(prefix, type, "_statByRBP.rds")
  saveRDS(resMergeStat.rbp, file = fileName)
  
  
  ##
  resMergeStat.ase <- resMergeAll %>% 
    group_by(idMap, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=length(unique(RBP)))
  
  resMergeStat.ase <- resMergeStat.ase %>% 
    split(resMergeStat.ase$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  fileName <- paste0(prefix, type, "_statByASE.rds")
  saveRDS(resMergeStat.ase, file = fileName)
}
resMergeSummaryFun(resMerge.ss, prefix, "spliceseq")


statByASE <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/tcga_ase_rmaps_iX_eX_region_site_separated_result_infoMerge_spliceseq_statByASE.rds")
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all_txAnno_repeatAnno.rds")

statByASE <- statByASE %>% 
  purrr::map(function(x) merge(x[,-2], tcga[tcga$txAnno!="notAnno",c("idMap", "txAnno")], by="idMap"))

statByASE <- statByASE %>% 
  purrr::map(function(x) x[,-1] %>% group_by(txAnno) %>% summarise_all(mean))

datals <- list("statByASEmean" = statByASE)

outDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/results/motifSummaryMean/"
dir.create(outDir)
source("/home/u1357/RNAseq/pancan/oncosplicingv3/peakEnrich/peakEnrich.R")
lapply(names(datals), plotFreqMapStatType, datals, outDir)