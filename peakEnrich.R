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
                  pValue = ifelse(inputRatio>controlRatio, logfdr, -logfdr), 
                  RBP=limma::strsplit2(RBPmotifId, split="_")[,1])
  
  return(df)
}



##### 2. process for summary stat #####
#
checkRegionBin <- function(name, dfls){
  # name <- "A3"
  # dfls <- resMergeStat.motifid
  checks <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_region_site.rds")
  checks <- limma::strsplit2(checks, split="\\.")
  checkls <- split(checks[,2], checks[,1])
  
  df <- dfls[[name]]
  check <- checkls[[name]]
  check <- check[order(check)]
  checkNO <- check[!check%in%names(df)]
  df[,checkNO] <- 0
  df <- df[,c(names(df)[1:2], check)]
}

#
resMergeStatFun <- function(resMerge, outDir){
  ##
  resMergeStat.motifid <- resMerge %>% 
    group_by(motifId, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=n()) 
  
  resMergeStat.motifid <- resMergeStat.motifid %>% 
    split(resMergeStat.motifid$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  name <- names(resMergeStat.motifid)
  resMergeStat.motifid <- names(resMergeStat.motifid) %>% purrr::map(function(x) checkRegionBin(x, resMergeStat.motifid))
  names(resMergeStat.motifid) <- name
  saveRDS(resMergeStat.motifid, file=paste0(outDir, "/", "resMergeStat.motifid.rds"))
  
  
  ##
  resMergeStat.rbp <- resMerge %>% 
    group_by(RBP, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=n())
  
  resMergeStat.rbp <- resMergeStat.rbp %>% 
    split(resMergeStat.rbp$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  name <- names(resMergeStat.rbp)
  resMergeStat.rbp <- names(resMergeStat.rbp) %>% purrr::map(function(x) checkRegionBin(x, resMergeStat.rbp))
  names(resMergeStat.rbp) <- name
  saveRDS(resMergeStat.rbp, file=paste0(outDir, "/", "resMergeStat.rbp.rds"))
  
  ##
  resMergeStat.ase <- resMerge %>% 
    group_by(idMap, spliceType, regionBinLoc) %>% 
    dplyr::summarise(freq=length(unique(RBP)))
  
  resMergeStat.ase <- resMergeStat.ase %>% 
    split(resMergeStat.ase$spliceType) %>% 
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, freq, fill=0))
  name <- names(resMergeStat.ase)
  resMergeStat.ase <- names(resMergeStat.ase) %>% purrr::map(function(x) checkRegionBin(x, resMergeStat.ase))
  names(resMergeStat.ase) <- name
  saveRDS(resMergeStat.ase, file=paste0(outDir, "/", "resMergeStat.ase.rds"))
  
  datals <- list(resMergeStat.motifid, resMergeStat.rbp, resMergeStat.ase)
  names(datals) <- c("statByMotifId", "statByRBP", "statByASE")
  
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
  if (maxdf==0){
    maxdf <- 1
    col_fun2 <- colorRamp2(c(-maxdf, -maxdf/2, 0, maxdf/2, maxdf), c("navy","DodgerBlue","white","Salmon","red4"))
  }else if (scaleType%in%c("row", "col")){
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
  names(df)[1] <- "RBPmotifId"
  row.names(df) <- df$RBPmotifId
  df <- df[, -c(seq(2))]
  
  
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
  # name <- "statByRBP"
  # datals <- freqls
  # name <- "input"
  outPref <- paste0(outDir, name, "/")
  dir.create(outPref)
  data <- datals[[name]]
  
  if (name=="statByASE"){
    data <- data %>% 
      purrr::map(function(x) x[,-1] %>% group_by(spliceType) %>% summarise_all(mean))
    data <- data %>% 
      purrr::map(function(x) cbind(x[,"spliceType", drop=F] %>% dplyr::rename(type=spliceType), x))
    
    scaleTypes <- c("mean")
  }else{
    scaleTypes <- c("row", "col", "none")
  }
  
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
  # gene <- names(dfls)[1]
  fillColor <- gsub("50","green4", gsub("250", "grey90", gsub("S", "red",columnTitle)))
  
  if (gene=="allRBP"){
    outPref <- gsub(spliceType, "", outPref)
  }
  df <- dfls[[gene]]
  row.names(df) <- df$RBPmotifId
  df <- df[, -c(seq(3))]
  
  pdf(file=paste0(outPref, spliceType, "_",  gene, ".pdf"), width=ncol(df)*0.23, height=nrow(df)*0.2+1.5)
  plotPvalMap(df, columnTitle, fillColor)
  dev.off()
  
  png(file=paste0(outPref, spliceType, "_",  gene, ".png"), width = ncol(df)*22, height = nrow(df)*20+150)
  plotPvalMap(df, columnTitle, fillColor)
  dev.off()
}

#
plotPvalMapSpliceType <- function(spliceType, resMergeStat.pval, rbp, outDir){
  # spliceType <- "ES"
  if (is.null(rbp)){
    outPref <- paste0(outDir, spliceType, "/")
    dir.create(outPref)
  }else{
    outPref <- outDir
  }
  
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
regionBinCheckPval <- function(resMergeAll, inputASE, controlASE=NULL, rbp=NULL, outDir, binMotifUnique=F){
  
  ## 0. ASE process
  message("processing input/control ASEs...")
  tcga <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/result/tcga_ase_all.rds")
  if (T){
    ## project = spliceseq
    allASE <- tcga$idMap[!is.na(tcga$idType)]
  }else{
    allASE <- tcga$idMap[!is.na(tcga$Splice_Event)]
  }
  

  if (is.null(controlASE)){
    controlASE <- allASE[!allASE%in%inputASE]
  }
  controlASEall <- controlASE
  controlASE <- idMapRmdup(controlASE)
  message("controlASEall:",length(controlASEall)," and controlASEfilter:", length(controlASE))
  controlAll <- data.frame(controlASE, spliceType=limma::strsplit2(controlASE, split="\\|")[,2]) %>% 
    dplyr::group_by(spliceType) %>% 
    summarise(controlAll=n())
  
  inputASEall <- inputASE
  inputASE <- idMapRmdup(inputASE)
  message("inputASEall:",length(inputASEall)," and inputASEfilter:", length(inputASE))
  inputAll <- data.frame(inputASE, spliceType=limma::strsplit2(inputASE, split="\\|")[,2]) %>% 
    dplyr::group_by(spliceType) %>% 
    summarise(inputAll=n()) %>% 
    dplyr::filter(inputAll>20)
  
  
  
  ## 1. data preprocess
  message("processing mapping data...")
  if (T){
    # resMergeAll <- project
    # project=="spliceseq"
    # resMergeAll <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_spliceseq.rds")
  }else{
    # project=="spladder"
    # resMerge <- readRDS(file = "/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/tcga_ase_rmaps_i250_e50_region_separated_result_infoMerge_spladder.rds")
  }
  
  if (!is.null(rbp)){
    resMerge <- resMergeAll %>% 
      dplyr::filter(RBP %in% rbp)
    dataTypeFilter <- c("input", "control")
  }else{
    resMerge <- resMergeAll
    dataTypeFilter <- c("input")
  }
  
  resMerge <- resMerge %>% 
    dplyr::mutate(dataType=ifelse(idMap%in%inputASE, "input", 
                                  ifelse(idMap%in%controlASE, "control", "other"))) %>% 
    dplyr::filter(dataType %in% c("input", "control"))
  
  
  # if (binMotifUnique==T){
  #   resMerge <- resMerge %>%
  #     dplyr::mutate(binDup=paste0(motifId, idType, regionBinLoc)) %>% 
  #     dplyr::filter(!duplicated(binDup))
  # }


  resMergeStat <- resMerge %>%
    dplyr::select(spliceType, regionBinLoc, dataType, RBP, motifId) %>% 
    # tidyr::gather(type, RBPmotifId,-c(seq(3))) %>% 
    # dplyr::select(-type) %>% 
    dplyr::mutate(RBPmotifId=RBP) %>% 
    group_by(RBPmotifId, spliceType, regionBinLoc, dataType) %>% 
    dplyr::summarise(freq=n())
  
  
  ## 2. check stat result
  message("checking data for statistic...")
  checkAll <- readRDS(file="/home/u1357/RNAseq/pancan/oncosplicingv3/rmaps3/check_RBPmotifId_spliceType_regionBins_dataType.rds")
  resMergeIn <- unique(paste(resMergeStat$RBPmotifId, resMergeStat$spliceType, sep="."))
  
  resMergeStatCheck <- paste(resMergeStat$RBPmotifId, resMergeStat$spliceType, resMergeStat$regionBinLoc, resMergeStat$dataType, sep=".")
  checks <- checkAll %>% 
    dplyr::filter(motifSpliceType%in%resMergeIn & !checks%in%resMergeStatCheck)
  resMergeStat <- rbind(resMergeStat, checks[,1:5])
  
  
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
    arrange(RBP, RBPmotifId, spliceType, regionBinLoc)
  
  
  ## 4. result report, set cutoff
  ## export using un-rmdup ASEs, inputASEall & controlASEall
  message("generate report data...")
  
  saveRDS(resMergeStat.pval, file=paste0(outDir, "resMergeStat.pval.rds"))
  
  resMergeStat.pval.filter <- resMergeStat.pval %>%
    # dplyr::filter(abs(pValue) > -log(0.0001, 10) & RBP!=RBPmotifId) %>%
    # dplyr::filter(abs(pValue) > -log(0.05, 10) & RBP!=RBPmotifId) %>%
    dplyr::filter(abs(pValue) > -log(0.05, 10)) %>%
    dplyr::mutate(term=paste(RBP, RBPmotifId, spliceType, regionBinLoc, sep="."))
  write.table(resMergeStat.pval.filter, file=paste0(outDir, "resMergeStat.pval.filter.txt"), sep="\t", row.names = F)

  resMergeStat.pval.result <- resMerge %>%
    dplyr::mutate(term=paste(RBP, motifId, spliceType, regionBinLoc, sep=".")) %>% 
    dplyr::filter(term%in%resMergeStat.pval.filter$term) %>%
    dplyr::mutate(dataType=ifelse(idMap%in%inputASEall, "input",
                                  ifelse(idMap%in%controlASEall, "control", "other"))) %>%
    dplyr::filter(dataType %in% dataTypeFilter)
  write.table(resMergeStat.pval.result, file=paste0(outDir, "resMergeStat.pval.result.txt"), sep="\t", row.names = F)


  ##  5. processing data for plot
  if (is.null(rbp)){
    message("processing data for plot: freq map")
    resMergeFreq <- resMergeAll %>%
      dplyr::filter(idMap %in% inputASE)
    datals <- resMergeStatFun(resMergeFreq, outDir)
  }

  
  ##
  message("processing data for plot: RBP freq")
  if (!is.null(rbp)){
    message("plotting freq map...")
    resMergeStat.input <- resMergeStat.pval %>%
      dplyr::select(RBP, RBPmotifId, spliceType, regionBinLoc, input)
    resMergeStat.input <- split(resMergeStat.input %>% dplyr::select(-RBP), resMergeStat.input$spliceType)
    resMergeStat.input <- resMergeStat.input %>%
      purrr::map(function(x) tidyr::spread(x, regionBinLoc, input))
    
    resMergeStat.control <- resMergeStat.pval %>%
      dplyr::select(RBP, RBPmotifId, spliceType, regionBinLoc, control)
    resMergeStat.control <- split(resMergeStat.control %>% dplyr::select(-RBP), resMergeStat.control$spliceType)
    resMergeStat.control <- resMergeStat.control %>%
      purrr::map(function(x) tidyr::spread(x, regionBinLoc, control))
    freqls <- list(input=resMergeStat.input, control=resMergeStat.control)
  }
 
  
  
  ##
  message("processing data for plot: pval map")
  resMergeStat.pval <- resMergeStat.pval %>%
    dplyr::select(RBP, RBPmotifId, spliceType, regionBinLoc, pValue)

  ##
  if (is.null(rbp)){
    rbps <- resMergeStat.pval[!grepl("_",resMergeStat.pval$RBPmotifId),]
    rbps <- split(rbps, rbps$spliceType)
    rbps <- rbps %>%
      purrr::map(function(x) tidyr::spread(x, regionBinLoc, pValue, fill=0))
  }

  ## split
  resMergeStat.pval <- split(resMergeStat.pval, resMergeStat.pval$spliceType)
  resMergeStat.pval <- resMergeStat.pval %>%
    purrr::map(function(x) tidyr::spread(x, regionBinLoc, pValue))

  resMergeStat.pval <- resMergeStat.pval %>%
    purrr::map(function(x) split(x, x$RBP))
  

  if (is.null(rbp)){
    for (i in names(resMergeStat.pval)){
      resMergeStat.pval[[i]][["allRBP"]] <- rbps[[i]]
    }
  }


  ## 6. plot freq map
  if (is.null(rbp)){
    message("plotting freq map...")
    lapply(names(datals), plotFreqMapStatType, datals, outDir)
  }
  
  
  ## 7. plot pval map
  message("plotting pval map...")
  lapply(names(resMergeStat.pval), plotPvalMapSpliceType, resMergeStat.pval, rbp, outDir)
  
  
  ## 8. plot RBP freq
  if (!is.null(rbp)){
    message("plotting RBP freq...")
    lapply(names(freqls), plotFreqMapStatType, freqls, outDir)
  }
  

  message("all done!")
}

