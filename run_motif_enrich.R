#############################################
######
###### 1. precision locations of motifs in ASE sequecnes
######
#############################################

## global: aseFa, outfile
## NOTE that non-global variant need to be indicated in each defined function
## or it will be reported as: "In mclapply(): all scheduled cores encountered errors in user code"

## define functions
listframe <- function(x, listLoc){
  listLoc[[x]] <- data.frame(listLoc[[x]]) %>% 
    plyr::mutate(seqName=x)
}


strLocate <- function(motifId, i){
  cat(motifId, "\n")
  motif <- motifId
  # motif <- limma::strsplit2(motifId, "_")[,3]
  listLoc <- stringr::str_locate_all(aseFa$seq, motif)
  names(listLoc) <- aseFa$idMap
  
  listLoc <- listLoc %>% 
    purrr::keep(function(x) nrow(x) > 0)
  
  listLoc <- names(listLoc) %>% 
    purrr::map(function(x) listframe(x, listLoc))
  
  loc <- data.frame(do.call(rbind, listLoc)) %>%
    plyr::mutate(motifId=motifId)
  
  outtemp <- paste0(dirname(outfile), "/temp/", strsplit2(basename(outfile), "\\.")[,1], i, ".txt")
  write.table(loc, outtemp, sep = "\t", append = T, row.names = F, col.names = F, quote = F)
}


strLocateAll <- function(i, motifIdSplit){
  motifIds <- motifIdSplit[[i]]
  lapply(motifIds, strLocate, i)
}


findmotif <- function(motifIdsAll, cores=20, outfile){
  st <- rep(seq(cores), ceiling(length(motifIdsAll)/cores))[1:length(motifIdsAll)]
  motifIdSplit <- split(motifIdsAll, st[order(st)])
  
  set.seed(0123456789)
  parallel::mclapply(seq(cores), strLocateAll, motifIdSplit, mc.cores = cores)
  
  ## merge processed data
  files <- list.files(paste0(dirname(outfile), "/temp/"), pattern =  paste0("^", limma::strsplit2(basename(outfile), "\\.")[,1]), full.names = T)
  result = lapply(files, data.table::fread)
  result = data.frame(do.call(rbind, result))
  # file.remove(files)
  saveRDS(result, outfile)
}


library(dplyr)
library(stringr)
library(tidyr)
library(limma)
library(parallel)


## motifs
motifs <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
motifIdsAll <- unique(motifs$motif)
motifIdsAll <- sample(motifIdsAll, length(motifIdsAll))

dir.create("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/temp/")
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/", pattern = "_separated_fasta.rds$", full.names = T)
for (file in files){
  # file <- files[1]
  message(file)
  len <- limma::strsplit2(file, split="_")[,4]
  len <- as.numeric(gsub("i", "", len))
  
  if (len==250){
    motifIdsAll.filter <- motifIdsAll
  }else if(len+1==7){
    motifIdsAll.filter <- motifIdsAll[nchar(motifIdsAll)>=(len+1)]
  }else{
    motifIdsAll.filter <- motifIdsAll[nchar(motifIdsAll)==(len+1)]
  }
  
  aseFa <- readRDS(file)
  outfile <- gsub("_fasta", "_result", file)
  findmotif(motifIdsAll.filter, cores=60, outfile)
}


##
files <- list.files("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/", pattern = "_separated_antisense_fasta.rds$", full.names = T)
for (file in files){
  # file <- files[1]
  message(file)
  len <- limma::strsplit2(file, split="_")[,4]
  len <- as.numeric(gsub("i", "", len))
  
  if (len==250){
    motifIdsAll.filter <- motifIdsAll
  }else if(len+1==7){
    motifIdsAll.filter <- motifIdsAll[nchar(motifIdsAll)>=(len+1)]
  }else{
    motifIdsAll.filter <- motifIdsAll[nchar(motifIdsAll)==(len+1)]
  }
  
  aseFa <- readRDS(file)
  outfile <- gsub("_fasta", "_result", file)
  findmotif(motifIdsAll.filter, cores=60, outfile)
}
cat("done")


#############################################
######
###### 2. relative loacations of motifs in structured ASE events
######    +++++/---------...---------/+++++...+++++/---------...---------/+++++
######    +: 5 bins of exon; /:splice site; and -:20 bins of intron
#############################################
##### 2.1 region #####
library(dplyr)
aseFa <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_fasta.rds")
aseFa <- aseFa %>% 
  dplyr::mutate(len = end-start,
                spliceType = limma::strsplit2(idMap, split="\\|")[,2])

aseFa <- aseFa %>% 
  dplyr::mutate(exonType=ifelse(spliceType=="A3", ifelse(regionType %in% c("r1","r4","r6"), "exon", "intron"),
                                ifelse(spliceType=="A5", ifelse(regionType %in% c("r1", "r3", "r6"), "exon", "intron"),
                                       ifelse(spliceType=="ES", ifelse(regionType %in% c("r1","r4","r5","r8"), "exon", "intron"),
                                              ifelse(spliceType=="RI", ifelse(regionType %in% c("r1","r4"), "exon", "intron"), 
                                                     ifelse(regionType %in% c("r1","r4","r5","r8","r9","r12"), "exon", "intron"))))),
                binNum = ifelse(exonType=="exon", 5, 25))


##
message("load result & merge to aseFa")
result <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result.rds")
names(result) <- c("mapStart","mapEnd", "idMap", "motif")
# head(result)
# res <- merge(result[1:10000,], aseFa, by="idMap", all.x=T)
res <- merge(result, aseFa, by="idMap", all.x=T)


##
message("merge to motif")
motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
motif <- motif %>% 
  dplyr::mutate(motifId = paste0(geneName, "_", motif)) %>% 
  dplyr::select(geneName, motifId, motif) %>% 
  dplyr::rename(RBP=geneName)
resMerge <- merge(res, motif, by="motif", all.x=T)
rm(res, result, motif, aseFa)
gc()


## RIGHT 2023.10.31
message("process resMerge0-3")
resMerge <- resMerge %>%
  dplyr::mutate(idMap = limma::strsplit2(idMap, split="::")[,1],
                binStart = ceiling(mapStart*binNum/len),
                binEnd = ceiling(mapEnd*binNum/len))

## binEnd-binStart=0: 1 bin, use first bin
## binEnd-binStart>1: >2 bin, short region, use first bin
resMerge1 <- resMerge %>%
  dplyr::filter((binEnd-binStart)==0) %>%
  dplyr::mutate(binLoc=binStart) %>%
  dplyr::select(-binStart, -binEnd)

resMerge2 <- resMerge %>%
  dplyr::filter((binEnd-binStart)==1) %>%
  tidyr::gather(binType, binLoc, -c(seq(ncol(resMerge)-2))) %>%
  dplyr::select(-binType)

resMerge3 <- resMerge %>%
  dplyr::filter((binEnd-binStart)>1)

saveRDS(resMerge1, "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_temp1.rds")
saveRDS(resMerge2, "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_temp2.rds")
saveRDS(resMerge3, "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_temp3.rds")
rm(resMerge, resMerge1, resMerge2)
gc()


message("resMerge3 spread")
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
cors <- 20
dfls <- splitBy(resMerge3, cors)
cors <- min(length(dfls), cors)
resMerge3 <- parallel::mclapply(seq(cors), binLocFun, dfls, mc.cores = cors)
resMerge3 <- do.call(rbind, resMerge3)
resMerge3 <- resMerge3 %>% 
  dplyr::select(-binStart, -binEnd)
saveRDS(resMerge3, "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_temp33.rds")

resMerge1 <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_temp1.rds")
resMerge2 <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_temp2.rds")
resMerge <- rbind(resMerge1, resMerge2[,colnames(resMerge1)], resMerge3[,colnames(resMerge1)])

rm(resMerge1, resMerge2, resMerge3)
gc()

message("final process and save")
resMerge <- resMerge %>% 
  dplyr::mutate(binLoc = stringr::str_pad(binLoc, width = 2, side = "left",pad = 0),
                regionType = stringr::str_pad(gsub("r","",regionType), width = 2, side = "left",pad = 0),
                regionBinLoc = paste0("r",regionType, "_b", binLoc)) %>% 
  dplyr::select(-binLoc,-regionType)
saveRDS(resMerge, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge.rds")


#############################################
######
###### 3. merge data (structured AS events & splice site)
######    
#############################################

## merge region & site 
suppressMessages(library(dplyr))
suppressMessages(library(lubridate))
resMerge <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge.rds")
resMergeSite <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps2/tcga_ase_rmaps_iX_eX_region_separated_result_infoMerge.rds")
resMergeAll <- rbind(resMerge[,names(resMergeSite)], resMergeSite)
rm(resMerge, resMergeSite)
gc()


## filter by spliceSeq & spladder
message("filtering spliceseq and spladder")
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all.rds")
ase.ss <- tcga[!is.na(tcga$idType), c("idMap", "idType")]
ase.sa <- tcga[!is.na(tcga$Splice_Event), c("idMap", "Splice_Event")] %>% 
  dplyr::rename(idType=Splice_Event)

resMerge.ss <- resMergeAll %>% 
  dplyr::select(idMap, spliceType, regionBinLoc, RBP, motifId, chr, start, end, strand, mapStart, mapEnd) %>% 
  dplyr::filter(idMap %in% ase.ss$idMap)

resMerge.sa <- resMergeAll %>% 
  dplyr::select(idMap, spliceType, regionBinLoc, RBP, motifId, chr, start, end, strand, mapStart, mapEnd) %>% 
  dplyr::filter(idMap %in% ase.sa$idMap)


# saveRDS(resMergeAll, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_sssa_more.rds")
# resMergeAll <- resMergeAll %>% 
#   dplyr::select(idMap, idType, spliceType, regionBinLoc, RBP, motifId)
# saveRDS(resMergeAll, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_sssa.rds")
rm(resMergeAll)
gc()


##
message("merge and saving spliceseq")
resMerge.ss <- merge(ase.ss, resMerge.ss, by="idMap")
gc()
saveRDS(resMerge.ss, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_spliceseq_more.rds")
gc()
resMerge.ss <- resMerge.ss %>% 
  dplyr::select(idMap, idType, spliceType, regionBinLoc, RBP, motifId)
saveRDS(resMerge.ss, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_spliceseq.rds")
gc()


##
message("merge and saving spladder")
resMerge.sa <- merge(ase.sa, resMerge.sa, by="idMap")
saveRDS(resMerge.sa, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_spladder_more.rds")
gc()
resMerge.sa <- resMerge.sa %>% 
  dplyr::select(idMap, idType, spliceType, regionBinLoc, RBP, motifId)
saveRDS(resMerge.sa, file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_spladder.rds")
gc()



#############################################
######
###### 4. motif enrichment
######    
#############################################


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
resMergeAll <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_spliceseq_more.rds")
resMergeAll <- resMergeAll.bak

##### 1. enrichment for RBP-related ASEs ##### 
motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.filter.ls1 <- corr.ss.filter %>% 
  dplyr::filter(RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>20)

corr.ss.filter.ls2 <- corr.ss.filter %>% 
  dplyr::filter(RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_all")) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>20)

corr.ss.filter.ls <- rbind(corr.ss.filter.ls1, corr.ss.filter.ls2)
corr.ss.filter.ls <- split(corr.ss.filter.ls$Map_Event, corr.ss.filter.ls$rbpCorrType)


enrichRBPaseMotif <- function(name, inputASEls, resMergeAll, dir, rbp="", controlASE=NULL, binMotifUnique=F){
  message(name)
  ## name <- names(inputASEls)[1]
  ## name <- names(aseTypels)[3]
  ## inputASEls <- aseTypels
  
  inputASE <- inputASEls[[name]]
  if (!is.null(rbp)){
    rbp <- limma::strsplit2(name, split="_")[,1]
  }
  outDir <- paste0(dir, "/", rbp, "/", name, "/")
  dir.create(outDir, recursive = T)
  
  source("/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/motifEnrich.R")
  regionBinCheckPval(resMergeAll, inputASE, controlASE, rbp, outDir,binMotifUnique)
}


## ALL RBP-CORRTYPE
dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/aseByRBPs/"
dir.create(dir, recursive = T)

inputRBPasels <-  corr.ss.filter.ls
parallel::mclapply(names(inputRBPasels), enrichRBPaseMotif, inputRBPasels, resMergeAll, 
                   dir=dir, rbp="", controlASE=NULL,
                   mc.cores = 20)



##### 2. enrichment for RBP-related NMD-ASEs ##### 
## controlASE should be non-correlated
corr.ss.filter <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/correlation/mapas_spliceseq_correlation_filtered.rds")
# nmd <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/correlation/as_tr_loc_sanky_top0_InOutAllFALSE_new_NMD.rds")
tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all_txAnno_repeatAnno.rds")
library(dplyr)
tcga <- tcga %>%
  dplyr::filter(!is.na(idType))


###### 2.1 NMD ######
nmd <- tcga$idMap[tcga$txAnno=="PCD_NMD"]
# other <-tcga$idMap[tcga$txAnno=="Other"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.NMD <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% nmd & RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.NMD$RBP))

inputASEls <- split(corr.ss.NMD$Map_Event, corr.ss.NMD$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdByRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)

# dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdByRBPsControlledByOther/"
# dir.create(dir, recursive = T)
# parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, dir, controlASE=other, mc.cores = 20)


###### 2.2 aluNMD  ######
aluNMD <- tcga$idMap[tcga$txAnno=="PCD_NMD" & tcga$aluAnno=="yes"]
# other <-tcga$idMap[tcga$txAnno=="Other"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.aluNMD <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% aluNMD & RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.aluNMD$RBP))

inputASEls <- split(corr.ss.aluNMD$Map_Event, corr.ss.aluNMD$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/aluNMDbyRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)


resDir1 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/aluNMDbyRBPs/"
panRbpPlot(resDir1)




###### 2.3 aluNoNMD  ######

# other <-tcga$idMap[tcga$txAnno=="Other"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.aluNoNMD <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% aluNoNMD & RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.aluNoNMD$RBP))

inputASEls <- split(corr.ss.aluNoNMD$Map_Event, corr.ss.aluNoNMD$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/aluNoNMDbyRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)



## pval & plot
resDir <- resDir1
resDir2 <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/aluNoNMDbyRBPs/"
panRbpPlot <- function(resDir){
  files <- list.files(resDir, pattern = "*pval.filter.txt$", recursive = T, full.names = T)
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
  
  pdfFile <- paste0(resDir, "/pan_rbp.pdf")
  pdf(file=pdfFile, width=ncol(result)*0.15, height=nrow(result)*0.2+1.5)
  plotPvalMap2(result, columnTitle, fillColor)
  dev.off()
}


columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
columnTitle <- c(columnTitle, columnTitle)
fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))


result1 <- result

panRbpPlot(resDir2)



###
pdfFile <- paste0(resDir, "/pan_rbpXXX.pdf")
pdf(file=pdfFile, width=ncol(result)*0.15, height=nrow(result)*0.2+1.5)


colnames(result1) <- paste0("nmd",colnames(result1))
result <- merge(result1 %>% dplyr::mutate(motifId=row.names(result1)),
                result2 %>% dplyr::mutate(motifId=row.names(result2)),
                all=T)
row.names(result) <- result$motifId
result <- result[,-1]
result[is.na(result)] <- 0




###### 2.3 aluAll  ######
aluAll <- tcga$idMap[tcga$aluAnno=="yes"]
# other <-tcga$idMap[tcga$txAnno=="Other"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.aluAll <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% aluAll & RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.aluAll$RBP))

inputASEls <- split(corr.ss.aluAll$Map_Event, corr.ss.aluAll$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/aluAllByRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)





###### 2.4 nmd noAlu ######
nmdNoAlu <- tcga$idMap[tcga$txAnno=="PCD_NMD" & tcga$aluAnno=="no"]
# other <-tcga$idMap[tcga$txAnno=="Other"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.nmdNoAlu <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% nmdNoAlu & RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.nmdNoAlu$RBP))

inputASEls <- split(corr.ss.nmdNoAlu$Map_Event, corr.ss.nmdNoAlu$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdNoAluByRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)


resDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdNoAluByRBPs/"
panRbpPlot(resDir)




###### 2.5 nmd noRE ######
nmdNoRE <- tcga$idMap[tcga$txAnno=="PCD_NMD" & tcga$repeatAnno=="no"]
# other <-tcga$idMap[tcga$txAnno=="Other"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.nmdNoRE <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% nmdNoRE & RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.nmdNoRE$RBP))

inputASEls <- split(corr.ss.nmdNoRE$Map_Event, corr.ss.nmdNoRE$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdNoREbyRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)


resDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdNoREbyRBPs/"
panRbpPlot(resDir)




###### 2.6 nmd RE ######
nmdRE <- tcga$idMap[tcga$txAnno=="PCD_NMD" & tcga$repeatAnno=="yes"]
# other <-tcga$idMap[tcga$txAnno=="Other"]

motif <- readRDS("/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/rmaps_rbsn_motif.rds")
rbpWithMotif <- unique(motif$geneName)[order(unique(motif$geneName))]

corr.ss.nmdRE <- corr.ss.filter %>% 
  dplyr::filter(Map_Event %in% nmdRE & RBP %in% rbpWithMotif) %>% 
  dplyr::mutate(rbpCorrType=paste0(RBP, "_", Corr_Type)) %>% 
  group_by(RBP, Splice_Type, Corr_Type) %>% 
  dplyr::filter(n()>50)
length(unique(corr.ss.nmdRE$RBP))

inputASEls <- split(corr.ss.nmdRE$Map_Event, corr.ss.nmdRE$rbpCorrType)

dir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdREbyRBPs/"
dir.create(dir, recursive = T)
parallel::mclapply(names(inputASEls), enrichRBPaseMotif, inputASEls, resMergeAll, 
                   dir,  rbp="", controlASE=NULL, 
                   mc.cores = 20)


resDir <- "/home/u1357/RNAseq/pancan/oncosplicingv3/motifEnrich/results/nmdREbyRBPs/"
panRbpPlot(resDir)

