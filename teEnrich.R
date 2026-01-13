library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)


##### 1. get pval by phyperTest #####
#
phyperTest <- function(vects){
  inputHits <- as.numeric(vects[[1]])
  hits <- as.numeric(vects[[2]])
  all <- as.numeric(vects[[3]])
  inputAll <- as.numeric(vects[[4]])
  p = phyper(inputHits-1, hits, all, inputAll, lower.tail=F, log.p = F)
  return(p)
}

#
getPval <- function(df){
  # df <- resMergeStat
  dfls <- split(df[,c("input", "hits", "all", "inputAll")], seq(nrow(df)))
  pval <- unlist(lapply(dfls, phyperTest))
  fdr <- p.adjust(pval, method = "fdr", n=length(pval))
  logfdr <- ifelse(fdr==0 | fdr < 1e-300, 300, -log(fdr, 10))
  inputRatio <- df$input/df$inputAll
  controlRatio <- df$control/df$controlAll
  
  df <- df %>% 
    dplyr::mutate(pval=pval,
                  FDR=fdr,
                  pValue = ifelse(inputRatio>controlRatio, logfdr, -logfdr))
  
  return(df)
}



##### 2. process for summary stat #####
#
checkRegionBin <- function(name, dfls, preColN){
  # name <- "ES"
  # dfls <- resMergeStat.repname
  checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
  checks <- limma::strsplit2(checks, split="\\.")
  checkls <- split(checks[,2], checks[,1])
  
  df <- dfls[[name]]
  check <- checkls[[name]]
  check <- check[order(check)]
  checkNO <- check[!check%in%names(df)]
  df[,checkNO] <- 0
  df <- df[,c(names(df)[1:preColN], check)]
}

#
resMergeStatFun <- function(resMerge, outDir){
  ##
  resMergeStat.repname <- resMerge %>% 
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
  resMergeStat.repfamily <- resMerge %>% 
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
  resMergeStat.repclass <- resMerge %>% 
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
  
  return(datals)
}



##### 3. plot freq map #####
#
plotFreqMap <- function(df, title, fillColor, scaleType){
  if (scaleType=="row"){
    df <- data.frame(t(scale(t(df))))
    df[is.na(df)] <- 0
  }else if(scaleType=="col"){
    df <- scale(df)
    df[is.na(df)] <- 0
  }else if(scaleType%in%c("mean", "none")){
    df <- df
  }
  maxdf <- max(abs(min(df)), max(df))
  split <- limma::strsplit2(colnames(df), split = "[0-9]{2}$")[,1]
  if (scaleType%in%c("row", "col")){
    col_fun2 = colorRamp2(c(-maxdf, -maxdf/2, 0, maxdf/2, maxdf), c("navy","DodgerBlue","white","Salmon","red4"))
  }else{
    col_fun2 = colorRamp2(c(0, maxdf/4, 2*maxdf/4, 3*maxdf/4, maxdf), c("white","Salmon","red","red4", "black"))
  }
  
  if (nrow(df)>150){
    rowNameFont = 1
  }else{
    rowNameFont = NULL
  }
  
  
  ht_opt$message = FALSE
  p <- Heatmap(df, col = col_fun2, cluster_rows=T,cluster_columns=F, 
               row_title = NULL, column_split = split, column_gap = unit(0, "mm"), border = T,
               column_title = title, column_title_gp = gpar(fill = fillColor, fontsize=15, fontface="italic"),
               column_title_side = "bottom", row_names_gp = gpar(fontsize = rowNameFont),
               show_row_names = T, rect_gp = gpar(col = NA), show_column_dend = F,
               show_row_dend = F, show_column_names = T,
               heatmap_legend_param = list(title = paste0(scaleType,"Scale")))
  draw(p)
}

#
plotFreqMapSpliceType <- function(spliceType, data, name, scaleTypes, outPref){
  # spliceType <- "ES"
  if (spliceType=="A3"){
    columnTitle <- c("50", "S", "250", "250", "S", "50", "250", "S", "50")
  }else if(spliceType=="A5"){
    columnTitle <- c("50", "S", "250", "50", "S", "250", "250", "S", "50")
  }else if(spliceType=="ES"){
    columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
  }else if(spliceType=="IR"){
    columnTitle <- c("50", "S", "250", "250", "S", "50")
  }else if(spliceType=="ME"){
    columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
  }
  fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))
  
  df <- data[[spliceType]] %>% as.data.frame()
  row.names(df) <- paste0(df[,1], "_", df[,2])
  df <- df[, -c(seq(3))]
  
  
  for (scaleType in scaleTypes){
    # scaleType <- "col"
    heightPdf <- ifelse(nrow(df)>150, 26, nrow(df)*0.2+1.5)
    widthPdf <- ncol(df)*0.23
    
    pdf(file=paste0(outPref, name, "_", spliceType, "_", scaleType, "scale.pdf"), width = widthPdf, height = heightPdf)
    plotFreqMap(df, columnTitle, fillColor, scaleType)
    dev.off()
    
    png(file=paste0(outPref, name, "_", spliceType, "_", scaleType, "scale.png"), width = widthPdf*100, height = heightPdf*100)
    plotFreqMap(df, columnTitle, fillColor, scaleType)
    dev.off()
  }
}

#
plotFreqMapStatType <- function(name, datals, outDir){
  # name <- "statByRepFamily"
  outPref <- paste0(outDir, name, "/")
  dir.create(outPref)
  data <- datals[[name]]
  
  scaleTypes <- c("row", "col", "none")
  
  lapply(names(data), plotFreqMapSpliceType, data, name, scaleTypes, outPref)
}



##### 4. plot pval map #####
#
plotPvalMap <- function(df, title, fillColor){
  # split <- limma::strsplit2(colnames(df), split = "_")[,1]
  split <- limma::strsplit2(colnames(df), split = "[0-9]{2}$")[,1]
  col_fun2 = colorRamp2(c(min(min(df), -10), -4, log(0.05, 10), 0, -log(0.05, 10), 4, 10, max(max(df), 20)), c("darkgreen", "green", "DodgerBlue4","DodgerBlue4","white","Salmon","red","red4"))
  p <- Heatmap(df, col = col_fun2, cluster_rows=T,cluster_columns=F,
               row_title = NULL, column_split = split, column_gap = unit(0, "mm"), border = T,
               column_title = title, column_title_gp = gpar(fill = fillColor, fontsize=15, fontface="italic"),
               column_title_side = "bottom", 
               show_row_names = T,rect_gp = gpar(col = NA), show_column_dend = F,
               show_row_dend = F, show_column_names = T,
               heatmap_legend_param = list(title = "-log10(FDR)"))
  draw(p)
}

#
plotPvalMapGene <- function(gene, dfls, spliceType, columnTitle, outPref){
  # gene <- "allClass"
  fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))
  
  if (gene%in%c("allClass", "allFamily", "sigName")){
    outPref <- gsub(paste0("/",spliceType,"/"), paste0("/",gene,"/"), outPref)
    dir.create(outPref)
  }
  df <- dfls[[gene]]
  row.names(df) <- paste0(df$repLevel, "_", df$senseType)
  df <- df[, -c(seq(5))]
  
  pdf(file=paste0(outPref, spliceType, "_",  gene, ".pdf"), width=ncol(df)*0.23, height=nrow(df)*0.2+1.5)
  plotPvalMap(df, columnTitle, fillColor)
  dev.off()
  
  png(file=paste0(outPref, spliceType, "_",  gene, ".png"), width = ncol(df)*22, height = nrow(df)*20+150)
  plotPvalMap(df, columnTitle, fillColor)
  dev.off()
}

#
plotPvalMapSpliceType <- function(spliceType, resMergeStat.pval, outDir){
  # spliceType <- "IR"
  outPref <- paste0(outDir, spliceType, "/")
  dir.create(outPref)
  
  ## set title & color
  if (spliceType=="A3"){
    columnTitle <- c("50", "S", "250", "250", "S", "50", "250", "S", "50")
  }else if(spliceType=="A5"){
    columnTitle <- c("50", "S", "250", "50", "S", "250", "250", "S", "50")
  }else if(spliceType=="ES"){
    columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
  }else if(spliceType=="IR"){
    columnTitle <- c("50", "S", "250", "250", "S", "50")
  }else if(spliceType=="ME"){
    columnTitle <- c("50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50", "50", "S", "250", "250", "S", "50")
  }
  
  ##
  dfls <- resMergeStat.pval[[spliceType]]
  
  lapply(names(dfls), plotPvalMapGene, dfls, spliceType, columnTitle, outPref)
}





##### 5. main #####
## rm duplicated input/control ASEs
idMapRmdup <- function(ase){
  dup1 <- gsub("(:).*(:)", ":", ase)
  dup2 <- gsub("(.*\\|.*\\|).*?-|-.+$", "\\1", ase)
  ase <- ase[!duplicated(dup1) & !duplicated(dup2)]
  return(ase)
}


## main
regionBinCheckPval <- function(resMergeAll, inputASE, controlASE=NULL, outDir, type, minInputASEs=20){
  ## 0. ASE process
  message("processing input/control ASEs...")
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
  
  allASE <- tcga$idMap
  
  ##
  inputASE <- inputASE[inputASE%in%allASE]
  if (is.null(controlASE)){
    controlASE <- allASE[!allASE%in%inputASE]
  }else{
    controlASE <- controlASE[controlASE%in%allASE]
  }
  controlASEall <- controlASE
  controlASE <- idMapRmdup(controlASE)
  message("controlASEall:",length(controlASEall)," and controlASEfilter:", length(controlASE))
  controlAll <- data.frame(controlASE, spliceType=limma::strsplit2(controlASE, split="\\|")[,2]) %>% 
    dplyr::group_by(spliceType) %>% 
    summarise(controlAll=n())
  
  ##
  inputASEall <- inputASE
  inputASE <- idMapRmdup(inputASE)
  message("inputASEall:",length(inputASEall)," and inputASEfilter:", length(inputASE))
  inputAll <- data.frame(inputASE, spliceType=limma::strsplit2(inputASE, split="\\|")[,2]) %>% 
    dplyr::group_by(spliceType) %>% 
    summarise(inputAll=n()) %>%
    dplyr::filter(inputAll>=minInputASEs)
  
  
  
  
  
  ## 1. data preprocess
  message("processing mapping data...")
  dataTypeFilter <- c("input")
  
  
  resMerge <- resMergeAll %>% 
    dplyr::mutate(repClassFamilyName=paste(repClass, repFamily, repName, sep="::"),
                  repClassFamily=paste(repClass, repFamily, sep="::"))
  
  
  resMerge <- resMerge %>% 
    dplyr::mutate(dataType=ifelse(idMap%in%inputASE, "input", 
                                  ifelse(idMap%in%controlASE, "control", "other"))) %>% 
    dplyr::filter(dataType %in% c("input", "control"))
  
  resMergeStat <- resMerge %>%
    dplyr::select(spliceType, regionBinLoc, dataType, senseType, repClassFamilyName, repClassFamily, repClass) %>% 
    tidyr::gather(type, repLevel,-c(seq(4))) %>% 
    dplyr::select(-type) %>% 
    group_by(repLevel, senseType, spliceType, regionBinLoc, dataType) %>% 
    dplyr::summarise(freq=n())
  
  
  ## 2. check stat result
  message("checking data for statistic...")
  # checkAll <- readRDS(file=paste0("/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/check_", repLevel, ".rds"))
  checkAll <- readRDS(file=paste0("/home/u1357/RNAseq/pancan/oncosplicingv3/analysis/te/maps/check_repAll.rds"))
  resMergeIn <- unique(paste(resMergeStat$repLevel, resMergeStat$senseType, resMergeStat$spliceType, sep="."))
  
  resMergeStatCheck <- paste(resMergeStat$repLevel, resMergeStat$spliceType, resMergeStat$regionBinLoc, resMergeStat$senseType, resMergeStat$dataType, sep=".")
  checks <- checkAll %>% 
    dplyr::filter(repLevelSenseSpliceType%in%resMergeIn & !checks%in%resMergeStatCheck)
  resMergeStat <- rbind(resMergeStat, checks[,names(resMergeStat)])
  
  
  ## 3. get pval
  message("geting p value...")
  resMergeStat <- merge(resMergeStat, inputAll, by="spliceType")
  resMergeStat <- merge(resMergeStat, controlAll, by="spliceType")
  resMergeStat <- resMergeStat %>% 
    tidyr::spread(dataType, freq, fill=0) %>% 
    as.data.frame() %>% 
    dplyr::mutate(all=inputAll+controlAll,
                  hits=input+control)
  
  ## by getPval function
  resMergeStat.pval <- getPval(resMergeStat)
  resMergeStat.pval <- resMergeStat.pval %>% 
    arrange(repLevel, spliceType, regionBinLoc)
  
  
  ## 4. result report, set cutoff
  ## export using un-rmdup ASEs, inputASEall & controlASEall
  message("generate report data...")
  
  resMergeStat.pval.filter <- resMergeStat.pval %>% 
    # dplyr::filter(abs(pValue) > -log(0.0001, 10) & RBP!=RBPmotifId) %>% 
    dplyr::filter(abs(pValue) > -log(0.05, 10)) %>%
    dplyr::mutate(term=paste(repLevel, senseType, spliceType, regionBinLoc, sep="."))
  write.table(resMergeStat.pval.filter, file=paste0(outDir, "resMergeStat.pval.filter.txt"), sep="\t", row.names = F)

  resMergeStat.pval.filter <- gsub("(::).*(\\..*sense.)","\\2",resMergeStat.pval.filter$term)
  resMergeStat.pval.result <- resMerge %>% 
    dplyr::filter(paste(repClass, senseType, spliceType, regionBinLoc, sep=".")%in%resMergeStat.pval.filter) %>%
    dplyr::mutate(dataType=ifelse(idMap%in%inputASEall, "input", 
                                  ifelse(idMap%in%controlASEall, "control", "other"))) %>% 
    dplyr::filter(dataType %in% dataTypeFilter) 
  write.table(resMergeStat.pval.result, file=paste0(outDir, "resMergeStat.pval.result.txt"), sep="\t", row.names = F)
  
  
  ##  5. processing data for plot
  message("processing data for plot: freq map")
  resMergeFreq <- resMergeAll %>% 
    dplyr::filter(idMap %in% inputASE)
  datals <- resMergeStatFun(resMergeFreq, outDir)
  
  
  message("processing data for plot: pval map")
  ## pValue = -log10(FDR)
  resMergeStat.pval <- resMergeStat.pval %>% 
    dplyr::select(repLevel, senseType, spliceType, regionBinLoc, pValue) %>% 
    dplyr::mutate(repClass=limma::strsplit2(repLevel, split="::")[,1],
                  repFamily=limma::strsplit2(repLevel, split = "::")[,2])
  
  ##
  class <- resMergeStat.pval[resMergeStat.pval$repFamily=="",]
  class <- split(class, class$spliceType)
  class <- class %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, pValue, fill=0))

  
  repname <- limma::strsplit2(resMergeStat.pval$repLevel, split="::")[,3]
  family <- resMergeStat.pval[resMergeStat.pval$repFamily!="" & repname=="",]
  family <- split(family, family$spliceType)
  family <- family %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, pValue, fill=0))
  
  
  name <- resMergeStat.pval[resMergeStat.pval$repFamily!="" & repname!="",] %>% 
    dplyr::filter(abs(pValue) > -log(0.05, 10))
  name <- split(name, name$spliceType)
  name <- name %>% 
    purrr::keep(function(x) nrow(x)>=1) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, pValue, fill=0))
  namename <- names(name)
  name <- names(name) %>% purrr::map(function(x) checkRegionBin(x, name, 5))
  names(name) <- namename
  
  ## split
  resMergeStat.pval <- resMergeStat.pval[resMergeStat.pval$repFamily!="",]
  resMergeStat.pval <- split(resMergeStat.pval, resMergeStat.pval$spliceType)
  resMergeStat.pval <- resMergeStat.pval %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, pValue))
  
  resMergeStat.pval <- resMergeStat.pval %>%
    purrr::map(function(x) split(x, x$repFamily))
  
  # xx <-resMergeStat.pval<-xx
  
  for (i in names(resMergeStat.pval)){
    resMergeStat.pval[[i]][["allClass"]] <- class[[i]]
    resMergeStat.pval[[i]][["allFamily"]] <- family[[i]]
    resMergeStat.pval[[i]][["sigName"]] <- name[[i]]
  }
  
  
  ## 6. plot freq map
  message("plotting freq map...")
  lapply(names(datals), plotFreqMapStatType, datals, outDir)
  
  
  ## 7. plot pval map
  message("plotting pval map...")
  lapply(names(resMergeStat.pval), plotPvalMapSpliceType, resMergeStat.pval, outDir)
  
  message("all done!")
}

