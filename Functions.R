
GetDGEObject <- function(cts,sampleInfo,geneInfo){

  
  # Reorder genes (rows) based on gene_info
  cts <- cts[match(geneInfo$ensembl_gene_id, rownames(cts)), ]
  
  # Reorder samples (columns) based on sample_info
  cts <- cts[, match(sampleInfo$sampleid, colnames(cts))]
  
  
  smallestSampleNumber <- nrow(sampleInfo)
  for(sample in unique(sampleInfo$samplename)){
    
    if(nrow(sampleInfo[sampleInfo$samplename==sample,]) < smallestSampleNumber){
      smallestSampleNumber <- nrow(sampleInfo[sampleInfo$samplename==sample,])
    }
    
    
  }
  
  
  #Filtering out genes with low counts in smallest sample group number
  flag <- rowSums(cts >= 10) >= smallestSampleNumber
  cts <- cts[flag,]
  
  #Filter gene information
  geneInfo <- geneInfo[flag,]
  
  
  # Reorder genes (rows) based on gene_info
  cts <- cts[match(geneInfo$ensembl_gene_id, rownames(cts)), ]
  
  # Reorder samples (columns) based on sample_info
  cts <- cts[, match(sampleInfo$sampleid, colnames(cts))]
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  #Obtain log CPM counts
  dge <- cpm(dge,log = TRUE, prior.count = 1)
  
  
  return(dge)
  
  
}

DOPCA <- function(cts,sampleInfo,geneInfo,tissueOfInterest,samplesToExclude = c(),internal = FALSE, pcaPrior = NULL, filteredSamples = c()){
  
  #Filter tissue data
  if(!internal){
    
    if(tissueOfInterest == "Small Bowel Mucosa"){
      
      sampleInfo <- sampleInfo[sampleInfo$tissuegeneral==tissueOfInterest,]
      samplesOfInterest <- sampleInfo$sampleid
      cts <- cts[,samplesOfInterest]
      plotWidth <- 24
      plotHeight <- 20
      titleSize <- 25
      labelSize <- 8
      
    }else{
      
      sampleInfo <- sampleInfo[sampleInfo$tissuedetail==tissueOfInterest,]
      samplesOfInterest <- sampleInfo$sampleid
      cts <- cts[,samplesOfInterest]
      plotWidth <- 12
      plotHeight <- 10
      titleSize <- 15
      labelSize <- 5
      
    }
    
    
    
  }
  
  
  # Reorder genes (rows) based on gene_info
  cts <- cts[match(geneInfo$ensembl_gene_id, rownames(cts)), ]
  
  # Reorder samples (columns) based on sample_info
  cts <- cts[, match(sampleInfo$sampleid, colnames(cts))]
  
  
  smallestSampleNumber <- nrow(sampleInfo)
  for(sample in unique(sampleInfo$samplename)){
    
    if(nrow(sampleInfo[sampleInfo$samplename==sample,]) < smallestSampleNumber){
      smallestSampleNumber <- nrow(sampleInfo[sampleInfo$samplename==sample,])
    }
    
    
  }
  
  
  #Filtering out genes with low counts in smallest sample group number
  flag <- rowSums(cts >= 10) >= smallestSampleNumber
  cts <- cts[flag,]
  
  #Filter gene information
  geneInfo <- geneInfo[flag,]
  
  
  # Reorder genes (rows) based on gene_info
  cts <- cts[match(geneInfo$ensembl_gene_id, rownames(cts)), ]
  
  # Reorder samples (columns) based on sample_info
  cts <- cts[, match(sampleInfo$sampleid, colnames(cts))]
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  
  #Library size plot to understand the variability between samples
  # ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=samplename)) +
  #   geom_bar(stat = "identity") +
  #   labs(y = "Library size (total number of mapped and quantified reads)",
  #        x = "Sample", fill = "samplename") +
  #   coord_flip()
  
  
  #MDS dimension reduction plot
  # glMDSPlot(dge,labels = dge$samples$sampleid, groups = dge$samples$samplename)
  
  
  
  
  #PCA Analysis
  pcaAnalysis <- prcomp(t(cpm(dge,log = TRUE,prior.count = 1)))
  pcaAnalysisSummary <- summary(pcaAnalysis)
  varianceExplained <- pcaAnalysisSummary$importance[2,]
  
  screeData <- data.frame(
    PC = factor(1:length(varianceExplained), levels = 1:length(varianceExplained)),
    Variance = varianceExplained
  )
  
  dirPath <- paste0("data/PCAplots")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  
  dirPath <- paste0("data/PCAplots/",tissueOfInterest)
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  
  
  #PCA plot
  if(internal){
    
    pcaPlot <- ggbiplot::ggbiplot(pcaAnalysis,
                                  groups = dge$samples$samplename,
                                  # labels = dge$samples$sampleid,
                                  var.axes = FALSE) +
      geom_text_repel(aes(label = dge$samples$sampleid), size = labelSize,max.overlaps = Inf) +
      ggtitle("After filtering") +   
      theme(plot.title = element_text(hjust = 0.5, size = titleSize),
            axis.title = element_text(size = titleSize),
            axis.text  = element_text(size = titleSize),
            legend.title = element_text(size = titleSize),
            legend.text = element_text(size = titleSize))
    
    plotCaption <- paste("Samples removed: ",paste(filteredSamples,collapse = ", "))
    
    pcaPrior <- pcaPrior + theme(plot.title = element_blank())
    pcaPrior <- pcaPrior + ggtitle("Before filtering") + theme(plot.title = element_text(hjust = 0.5))
    
    dualPlot <- pcaPrior + pcaPlot +
      plot_annotation(
        title = paste0("PCA for ",tissueOfInterest," tissue group before and after removing certain samples"),
        caption = plotCaption
        )& theme(
          plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5)
        )
    
    ggsave(filename = paste0(dirPath,"/",tissueOfInterest,"_filtered.jpeg"),dualPlot,width = plotWidth, height = plotHeight, dpi = 600)
    
  }else{
    
    
    
    pcaPlot <- ggbiplot::ggbiplot(pcaAnalysis,
                                  groups = dge$samples$samplename,
                                  # labels = dge$samples$sampleid,
                                  var.axes = FALSE) +
      geom_text_repel(aes(label = dge$samples$sampleid), size = labelSize,max.overlaps = Inf) +
      ggtitle(paste0("PCA Plot of RNA-seq Data for ",tissueOfInterest," tissue")) +   
      theme(plot.title = element_text(hjust = 0.5, size = titleSize),
            axis.title = element_text(size = titleSize),
            axis.text  = element_text(size = titleSize),
            legend.title = element_text(size = titleSize),
            legend.text = element_text(size = titleSize)
            )
    
    ggsave(filename = paste0(dirPath,"/PCA.jpeg"),pcaPlot,width = plotWidth, height = plotHeight, dpi = 600)
    
    screePlot <- ggplot(screeData, aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", fill = "pink", color = "black", alpha = 0.7) +
      labs(
        title = paste0("Scree Plot of PCA for ",tissueOfInterest," tissue"),
        x = "Principal Component",
        y = "Proportion of Variance Explained"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
    
    ggsave(filename = paste0(dirPath,"/Scree.jpeg"),screePlot,width = plotWidth, height = plotHeight, dpi = 600)
    
  }
  

  #Redo PCA if there are any samples to be removed
  samplesRemoved <- colnames(cts)[colnames(cts) %in% samplesToExclude]
  if(!length(samplesRemoved) == 0){
    
    #Remove certain samples to redo PCA
    cts <- cts[, !colnames(cts) %in% samplesToExclude]
    sampleInfo <- sampleInfo[!sampleInfo$sampleid %in% samplesToExclude,]
    
    #Call DoPCA
    DOPCA(cts = cts, sampleInfo = sampleInfo,
          geneInfo = geneInfo,
          tissueOfInterest = tissueOfInterest,
          internal = TRUE,
          pcaPrior = pcaPlot,
          filteredSamples = samplesRemoved
          )
    
  }
  
  
  
}

DoGroupCorrelationPlots <- function(cts,sampleInfo,geneInfo,group){
  
  groupSampleInfo <- sampleInfo[sampleInfo$samplename == group,]
  
  groupCts <- cts[,colnames(cts) %in% groupSampleInfo$sampleid]
  
  dge <- GetDGEObject(groupCts,groupSampleInfo,geneInfo)
  
  corr <- round(cor(dge),1)
  
  corplot <-  ggcorrplot(corr, hc.order = TRUE, outline.color = "white",lab = TRUE) +
    ggtitle(paste0("Correlation Plot of RNA-seq Data for ",group," sample group"))+
    theme(plot.title = element_text(hjust = 0.5))
  
  dirPath <- paste0("data/GroupCorPlots")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  
  dirPath <- paste0("data/GroupCorPlots/",group)
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  ggsave(filename = paste0(dirPath,"/",group,".jpeg"),corplot,width = 12,height = 10,dpi = 600)
  
  
}

DoComprehensiveCorrelationPlot <- function(cts,sampleInfo,geneInfo){
  
  
  dge <- GetDGEObject(cts,sampleInfo,geneInfo)
  
  corr <- round(cor(dge),1)
  
  corplot <-  ggcorrplot(corr, hc.order = TRUE, outline.color = "white",lab = TRUE) +
    ggtitle(paste0("Correlation Plot of RNA-seq Data for all sample groups"))+
    theme(plot.title = element_text(hjust = 0.5))
  
  dirPath <- paste0("data/GroupCorPlots")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  ggsave(filename = paste0(dirPath,"/","ComprehensiveCorrelationPlot.jpeg"),corplot,width = 42,height = 30,dpi = 600)
  
  
}

DoCrossCorrelationPlot <- function(cts,sampleInfo,geneInfo,group1,group2){
  
  sampleIdsGroup1 <- sampleInfo[sampleInfo$samplename == group1,]$sampleid
  sampleIdsGroup2 <- sampleInfo[sampleInfo$samplename == group2,]$sampleid
  
  sampleIds <- c(sampleIdsGroup1,sampleIdsGroup2)
  
  groupSampleInfo <- sampleInfo[sampleInfo$sampleid %in% sampleIds,]
  
  groupCts <- cts[,colnames(cts) %in% sampleIds]
  
  
  dge <- GetDGEObject(groupCts,groupSampleInfo,geneInfo)
  
  corr <- round(cor(dge),1)
  
  corplot <-  ggcorrplot(corr, hc.order = TRUE, outline.color = "white",lab = TRUE) +
    ggtitle(paste0("Cross correlation Plot of RNA-seq Data between groups ",group1," and ",group2)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  dirPath <- paste0("data/CrossCorPlots")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  ggsave(filename = paste0(dirPath,"/",group1,"_",group2,".jpeg"),corplot,width = 12,height = 10,dpi = 600)
  
  
  
  
}

FilterLowQualitySamples <- function(cts,sampleInfo,geneInfo){
  
  sampleSymbols <- c("HF_C","LF_C","HF_S","LF_S")
  sampleGroups <- c("Highfat_Control","Lowfat_Control","Highfat_Surgery","Lowfat_Surgery")
  tissueGroups <- c("Duodenum","Jejunum","Ileum","Ileum interposition","Epididymal white adipose","Liver","Subcutaneous brown adipose","Subcutaneous white adipose")
  tissueSymbols <- c("_D_","_J_","_I_","_II_","_EA_","_L_","_SBr_","_SW_")
  
  lowQualitySamples <- c()
  
  sampleSummary <- as.data.frame(matrix(ncol = 5,nrow = 0))
  colnames(sampleSummary) <- c("Tissue",sampleGroups)
  
  for (symbol in tissueSymbols) {
    
    sampleSummaryRow <- as.data.frame(matrix(ncol = 5))
    colnames(sampleSummaryRow) <- c("Tissue",sampleGroups)
    sampleSummaryRow$Tissue <- tissueGroups[which(tissueSymbols==symbol)]
    
    for (group in sampleSymbols) {
      
      
      groupSamples <- grep(group,colnames(cts),value = TRUE)
      groupSamples <- grep(symbol,groupSamples,value = TRUE)
      
      if(length(groupSamples) == 0) next
      
      groupSampleInfo <- sampleInfo[sampleInfo$sampleid %in% groupSamples,]
      
      groupCts <- cts[,colnames(cts) %in% groupSampleInfo$sampleid]
      
      dge <- GetDGEObject(groupCts,groupSampleInfo,geneInfo)
      
      corr <- round(cor(dge),1)
      
      corrdf <- as.data.frame(corr)
      
      samplesWithLowValues <- colnames(corrdf)[colSums(corrdf < 0.8) >= 2]
      
      if( length(samplesWithLowValues)!=0){
        plotCaption <- paste0("Samples with low correlation: ",paste(samplesWithLowValues,collapse = ", "))
        lowQualitySamples <- c(lowQualitySamples,samplesWithLowValues)
      }else{
        plotCaption <- "All samples have good correlation with each other"
      } 
      
      
      corplot <-  ggcorrplot(corr, hc.order = TRUE, outline.color = "white",lab = TRUE) +
        plot_annotation(
          title = paste0(paste0("Correlation Plot of RNA-seq Data for ",tissueGroups[which(tissueSymbols==symbol)]," tissue group among ",sampleGroups[which(sampleSymbols==group)]," samples")),
          caption = plotCaption
        )& theme(
          plot.title = element_text(hjust = 0.5),
          plot.caption = element_text(hjust = 0.5)
        )
      
      dirPath <- paste0("data/TissueConditionCor")
      # Check if the directory exists to save plot
      if (!dir.exists(dirPath)) {
        # Create the directory
        dir.create(dirPath)
      }
      
      
      
      dirPath <- paste0("data/TissueConditionCor/",tissueGroups[which(tissueSymbols==symbol)],"_",sampleGroups[which(sampleSymbols==group)])
      # Check if the directory exists to save plot
      if (!dir.exists(dirPath)) {
        # Create the directory
        dir.create(dirPath)
      }
      
      print(corplot)
      
      ggsave(filename = paste0(dirPath,"/",tissueGroups[which(tissueSymbols==symbol)],"_",sampleGroups[which(sampleSymbols==group)],".jpeg"),corplot,width = 12,height = 10,dpi = 600)
      
      
      sampleSummaryRow[1,sampleGroups[which(sampleSymbols==group)]] <- length(groupSamples) - length(samplesWithLowValues)
      
    }
    
    sampleSummary <- rbind(sampleSummary,sampleSummaryRow)
    
  }
  
  write.csv(sampleSummary,file = "data/SampleSummary.csv",row.names = FALSE )
  
  return(lowQualitySamples)
  
}

Generate_Density_N_Library_Plot_Grid <- function(allCts,allSampleInfo,geneInfo){
  
  tissues <- unique(allSampleInfo$tissuedetail)
  
  librarySizePlotList <- list()
  
  densityPlotList <- list()
  
  for (i in 1:length(tissues)) {
    
    
    sampleInfo <- allSampleInfo[allSampleInfo$tissuedetail == tissues[i],]
    
    cts <- allCts[,colnames(allCts) %in% sampleInfo$sampleid]
    
    # Reorder genes (rows) based on gene_info
    cts <- cts[match(geneInfo$ensembl_gene_id, rownames(cts)), ]
    
    # Reorder samples (columns) based on sample_info
    cts <- cts[, match(sampleInfo$sampleid, colnames(cts))]
    
    #Create an interaction of diet and treatment
    sampleInfo$diet_treatment <- interaction(sampleInfo$diet,sampleInfo$treatment)
    
    
    smallestSampleNumber <- nrow(sampleInfo)
    for(sample in unique(sampleInfo$samplename)){
      
      if(nrow(sampleInfo[sampleInfo$samplename==sample,]) < smallestSampleNumber){
        smallestSampleNumber <- nrow(sampleInfo[sampleInfo$samplename==sample,])
      }
      
      
    }
    
    
    
    #Create a DGE object using the counts, sampleInfo and geneInfo
    dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
    
    #Normalization of the counts using TMM method
    dge <- calcNormFactors(dge,method = "TMM")
    
    #Library size plot to understand the variability between samples
    libSizePlot <- ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=diet_treatment)) +
      geom_bar(stat = "identity") +
      labs(y = "Library size (total number of mapped and quantified reads)",
           x = "Sample", fill = "Group") +
      ggtitle(tissues[i])+
      theme(plot.title = element_text(hjust = 0.5))+
      coord_flip() 
    
    librarySizePlotList[[i]] <- libSizePlot
    
    
    #Variability between samples
    max(dge$samples$lib.size)/min(dge$samples$lib.size)
    
    #Set a threshold for filtering genes
    flag <- (rowSums(cpm(dge,log = TRUE, prior.count = 1) > 1) >= smallestSampleNumber)
    
    
    #Generate a density plot before filtering
    beforeFiltering_plot <- dge %>% 
      cpm(log = TRUE, prior.count = 1) %>% 
      melt %>% 
      dplyr::filter(is.finite(value)) %>% 
      ggplot(aes(x = value, colour = Var2)) +
      geom_density() + 
      guides(colour = FALSE,scale = 'none') +
      ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes"))+
      labs(x = "logCPM", y = "Density")
    
    #Generate a Density plot after filtering
    afterFiltering_plot <- dge %>% 
      cpm(log = TRUE, prior.count = 1) %>% 
      magrittr::extract(flag,) %>%
      melt %>% 
      dplyr::filter(is.finite(value)) %>% 
      ggplot(aes(x = value, colour = Var2)) +
      geom_density() + 
      guides(colour = FALSE,scale = 'none') +
      ggtitle("B. After filtering", subtitle = paste0(table(flag)[[2]], " genes"))+
      labs(x = "logCPM", y = "Density")
    
    #plot both the density plots together
    densityPlot <- cowplot::plot_grid(beforeFiltering_plot, afterFiltering_plot) +
      ggtitle(tissues[i])+
      theme(plot.title = element_text(hjust = 0.5))
    
    densityPlotList[[i]] <- densityPlot
  
  }
  
  librarySizePlotGrid <- plot_grid(plotlist =  librarySizePlotList) + 
    ggtitle("Tissue wise library sizes of samples")+
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave("data/LibrarySizePlotGrid.jpeg", librarySizePlotGrid, width = 18, height = 18,dpi = 300)
  
  densityPlotGrid <- plot_grid(plotlist = densityPlotList)+ 
    ggtitle("Tissue wise density plots")+
    theme(plot.title = element_text(hjust = 0.5,  size = 20))
  ggsave("data/DensityPlotGrid.jpeg", densityPlotGrid, width = 18, height = 18,dpi = 300)
  
}

SaveLogCPMCounts <- function(logCPM,dirPath){
  
  write.csv(logCPM,file = paste0(dirPath,"/logCPMCounts.csv"),row.names = TRUE)
  
}

SaveAllNSigGenes <- function(fit,contrast,geneInfo,dirPath){
  
  allDEGenes <- topTable(fit,coef = contrast,number = Inf, adjust.method = "fdr")
  
  allDEGenes$EnsemblID <-  rownames(allDEGenes)
  ensemblIds <- allDEGenes$EnsemblID
  
  geneNames <- geneInfo$gene_name[match(ensemblIds,geneInfo$ensembl_gene_id)]
  geneNames[is.na(geneNames)] <- ensemblIds[is.na(geneNames)]
  allDEGenes$GeneSymbol <- geneNames

  rowNamesallDEGenes <- paste0(allDEGenes$EnsemblID,"_",allDEGenes$GeneSymbol)
  rownames(allDEGenes) <- rowNamesallDEGenes

  allDEGenes <- allDEGenes[,c("EnsemblID","GeneSymbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")]

  allDEGenes <- allDEGenes %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE,
      TRUE ~ FALSE
    ))
  
  file <- paste0(dirPath,"/AllGenes.xlsx")
  
  if (file.exists(file)) {
    # Load the existing workbook
    wb <- loadWorkbook(file)
    
    # Add a new sheet (check if the sheet already exists first)
    if (!(contrast %in% names(wb))) {
      addWorksheet(wb, contrast)
    }
    
    # Optionally, write some data to the new sheet
    writeData(wb, contrast,allDEGenes, startCol = 1, startRow = 1,rowNames = TRUE)
    
    # Save the workbook
    saveWorkbook(wb, file, overwrite = TRUE)
    
  } else {
    
    # Create a new workbook if the file doesn't exist
    wb <- createWorkbook()
    
    # Add a new sheet
    addWorksheet(wb, contrast)
    
    # Optionally, write some data to the new sheet
    writeData(wb, contrast,allDEGenes, startCol = 1, startRow = 1, rowNames = TRUE)
    
    # Save the new workbook
    saveWorkbook(wb, file, overwrite = TRUE)
  }
  
  
  sigDEResults <- allDEGenes %>%
    dplyr::filter(isSignificant == TRUE)
  
  if(nrow(sigDEResults)!=0){
    
    file <- paste0(dirPath,"/SigGenes.xlsx")
    
    if (file.exists(file)) {
      # Load the existing workbook
      wb <- loadWorkbook(file)
      
      # Add a new sheet (check if the sheet already exists first)
      if (!(contrast %in% names(wb))) {
        addWorksheet(wb, contrast)
      }
      
      # Optionally, write some data to the new sheet
      writeData(wb, contrast,sigDEResults, startCol = 1, startRow = 1,rowNames = TRUE)
      
      # Save the workbook
      saveWorkbook(wb, file, overwrite = TRUE)
      
    } else {
      
      # Create a new workbook if the file doesn't exist
      wb <- createWorkbook()
      
      # Add a new sheet
      addWorksheet(wb, contrast)
      
      # Optionally, write some data to the new sheet
      writeData(wb, contrast,sigDEResults, startCol = 1, startRow = 1, rowNames = TRUE)
      
      # Save the new workbook
      saveWorkbook(wb, file, overwrite = TRUE)
    }
    
  } 
  
}

PlotTrueNullFractions <- function(dirPath,contrasts){
  
  trueNULLDF <- as.data.frame(matrix(ncol = 2, nrow = 0))
  colnames(trueNULLDF) <- c("Contrast","TrueNullProportion")
  for (contrast in contrasts) {
    
    allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
    
    trueNull <- propTrueNull(allGenes$P.Value, method="lfdr", nbins=20)
    
    trueNULLDFRow <- data.frame(contrast,trueNull)
    colnames(trueNULLDFRow) <- c("Contrast","TrueNullProportion")
    
    trueNULLDF <- rbind(trueNULLDF,trueNULLDFRow)
    
  }
  
  jpeg(filename = paste0(dirPath,"TrueNullCol.jpeg"),width = 1000,height = 800,quality = 100)
  
  trueNullPlot <- ggplot(data = trueNULLDF,aes(x=Contrast, y=TrueNullProportion))+
    geom_col()
  
  print(trueNullPlot)
  
  dev.off()
  
}

PlotBioNTechCorr <- function(fit,contrasts,dirPath){
  
  contrastCombinations <- combn(contrasts, 2)
  
  contrastCombinations <- as.data.frame(t(contrastCombinations))
  
  if(nrow(contrastCombinations)%%2 == 0) rows <- nrow(contrastCombinations)/2
  if(nrow(contrastCombinations)%%2 != 0) rows <- (nrow(contrastCombinations)/2) + 1
  cols <- 2
  par(mfrow = c(rows, cols))
  
  jpeg(filename = paste0(dirPath,"/Genas.jpeg"),width = 1000,height = 500 * rows,quality = 100)
  
  for (i in 1:nrow(contrastCombinations)) {
    
    genas(fit,
          coef = c(contrastCombinations[i,1],
                   contrastCombinations[i,2]
          ),
          subset = "all",
          plot = TRUE,
          alpha = 0.4
    )
    
  }
  
  dev.off()
  
}

#Venn with adjusted p value
# PlotVennDiagrams <- function(fit,contrasts,geneInfo,tissueOfInterest,dirPath){
#   
#   
#   contrasts <- contrasts[contrasts!="Surgery_Diet"]
#   
#   decidedTests <- decideTests(fit,method = "separate",adjust.method = "fdr")
#   
#   decidedTests <- decidedTests[,colnames(decidedTests) != "Surgery_Diet"]
#   
#   vennCounts <- vennCounts(decidedTests)
#   
#   vennCounts <- as.data.frame(as.table(vennCounts))
#   
#   vennCounts <- pivot_wider(vennCounts,names_from = Var2, values_from = Freq)
#   
#   vennCounts <- vennCounts[,-1]
#   
#   write.csv(vennCounts,paste0(dirPath,"VennCounts.csv"),row.names = FALSE)
#   
#   vennCounts <- vennCounts[vennCounts$Counts!=0,]
#   
#   vennCounts <- vennCounts[!(vennCounts$HF_Surgery == 0 &
#                                vennCounts$LF_Surgery == 0 &
#                                vennCounts$Diet == 0 &
#                                vennCounts$Surgery == 0),]
#   
#   for (i in 1:nrow(vennCounts)) {
#     
#     selectedContrasts <- colnames(vennCounts)[which(vennCounts[i,(1:4)] == 1)]
#     nonSelectedContrasts <- setdiff(colnames(decidedTests), selectedContrasts)
#     
#     filteredGenes <- decidedTests[
#           rowSums(abs(decidedTests[, selectedContrasts, drop = FALSE]) == 1) == length(selectedContrasts) &  
#             rowSums(abs(decidedTests[, nonSelectedContrasts, drop = FALSE]) == 1) == 0, 
#         ]
#      
#     commonGenes <- rownames(filteredGenes)
#     
#     commonGeneSymbols <- geneInfo[geneInfo$ensembl_gene_id %in% commonGenes,]$gene_name
#     
#     if(i==1){
#       # Create a new workbook if the file doesn't exist
#       wb <- createWorkbook()
#       
#       # Add a new sheet
#       firstCharacters <- substr(selectedContrasts,1,1)
#       sheetName <- paste(firstCharacters,collapse = "_")
#       addWorksheet(wb, sheetName)
#       
#       # Optionally, write some data to the new sheet
#       writeData(wb, sheetName,commonGeneSymbols, startCol = 1, startRow = 1)
#       
#       # Save the new workbook
#       saveWorkbook(wb, paste0(dirPath,"VennGeneList.xlsx"), overwrite = TRUE)
#     }else{
#       
#       # Load the existing workbook
#       wb <- loadWorkbook(paste0(dirPath,"VennGeneList.xlsx"))
#       
#       # Add a new sheet (check if the sheet already exists first)
#       firstCharacters <- substr(selectedContrasts,1,1)
#       sheetName <- paste(firstCharacters,collapse = "_")
#       if (!(sheetName %in% names(wb))) {
#         addWorksheet(wb, sheetName)
#       }
#       
#       # Optionally, write some data to the new sheet
#       writeData(wb, sheetName,commonGeneSymbols, startCol = 1, startRow = 1)
#       
#       # Save the workbook
#       saveWorkbook(wb, paste0(dirPath,"VennGeneList.xlsx"), overwrite = TRUE)
#       
#     }
#     
#   }
#   
#   baseColors <- c("red","blue","green")
#   
#   custom_palette <- colorRampPalette(baseColors)
#   
#   colors <- custom_palette(length(contrasts))
#   
#   jpeg(filename = paste0(dirPath,"Venn.jpeg"),width = 1000,height = 800,quality = 100)
#   
#   
#   vennDiagram(decidedTests,
#               lwd = 1,
#               circle.col = colors,
#               names = contrasts,
#               main = paste0("Venn diagram for various contrasts in ",tissueOfInterest," tissue")
#   )
#   
#   
#   dev.off()
#   
# }


#Venn with nominal p value
PlotVennDiagrams <- function(fit, contrasts, geneInfo, tissueOfInterest, dirPath) {
  
  # Remove "Surgery_Diet" from contrasts
  contrasts <- contrasts[contrasts != "Surgery_Diet"]
  
  # Nominal p-value filtering (0.05) without adjustment
  decidedTests <- ifelse(fit$p.value < 0.05, 
                         ifelse(fit$coef > 0, 1, -1), 
                         0)
  decidedTests <- as.matrix(decidedTests)
  colnames(decidedTests) <- colnames(fit$p.value)
  rownames(decidedTests) <- rownames(fit$p.value)
  
  # Remove "Surgery_Diet" column
  decidedTests <- decidedTests[, colnames(decidedTests) != "Surgery_Diet"]
  
  # Venn counts
  vennCounts <- vennCounts(decidedTests)
  vennCounts <- as.data.frame(as.table(vennCounts))
  vennCounts <- tidyr::pivot_wider(vennCounts, names_from = Var2, values_from = Freq)
  vennCounts <- vennCounts[, -1]
  
  # Save Venn counts
  write.csv(vennCounts, paste0(dirPath, "VennCounts.csv"), row.names = FALSE)
  
  # Filter out zero counts and cases where all are zero
  vennCounts <- vennCounts[vennCounts$Counts != 0, ]
  vennCounts <- vennCounts[!(vennCounts$HF_Surgery == 0 &
                               vennCounts$LF_Surgery == 0 &
                               vennCounts$Diet == 0 &
                               vennCounts$Surgery == 0), ]
  
  # Loop through Venn combinations and save gene lists
  for (i in 1:nrow(vennCounts)) {
    
    selectedContrasts <- colnames(vennCounts)[which(vennCounts[i, (1:4)] == 1)]
    nonSelectedContrasts <- setdiff(colnames(decidedTests), selectedContrasts)
    
    filteredGenes <- decidedTests[
      rowSums(abs(decidedTests[, selectedContrasts, drop = FALSE]) == 1) == length(selectedContrasts) &
        rowSums(abs(decidedTests[, nonSelectedContrasts, drop = FALSE]) == 1) == 0, 
    ]
    
    commonGenes <- rownames(filteredGenes)
    commonGeneSymbols <- geneInfo[geneInfo$ensembl_gene_id %in% commonGenes, ]$gene_name
    
    if (i == 1) {
      wb <- openxlsx::createWorkbook()
      firstCharacters <- substr(selectedContrasts, 1, 1)
      sheetName <- paste(firstCharacters, collapse = "_")
      openxlsx::addWorksheet(wb, sheetName)
      openxlsx::writeData(wb, sheetName, commonGeneSymbols, startCol = 1, startRow = 1)
      openxlsx::saveWorkbook(wb, paste0(dirPath, "VennGeneList.xlsx"), overwrite = TRUE)
    } else {
      wb <- openxlsx::loadWorkbook(paste0(dirPath, "VennGeneList.xlsx"))
      firstCharacters <- substr(selectedContrasts, 1, 1)
      sheetName <- paste(firstCharacters, collapse = "_")
      if (!(sheetName %in% names(wb))) {
        openxlsx::addWorksheet(wb, sheetName)
      }
      openxlsx::writeData(wb, sheetName, commonGeneSymbols, startCol = 1, startRow = 1)
      openxlsx::saveWorkbook(wb, paste0(dirPath, "VennGeneList.xlsx"), overwrite = TRUE)
    }
    
  }
  
  # Color palette
  baseColors <- c("red", "blue", "green")
  custom_palette <- colorRampPalette(baseColors)
  colors <- custom_palette(length(contrasts))
  
  # # Save Venn diagram image
  # jpeg(filename = paste0(dirPath, "Venn.jpeg"), width = 1000, height = 900, quality = 100)
  # vennDiagram(decidedTests,
  #             lwd = 1,
  #             circle.col = colors,
  #             names = contrasts,
  #             main = paste0("Venn diagram for various contrasts in ", tissueOfInterest, " tissue. Pval < 0.05")
  # )
  # dev.off()
  
  
  # Save Venn diagram image
  jpeg(filename = paste0(dirPath, "Venn.jpeg"), width = 1000, height = 800, quality = 100)
  
  
  # Increase top margin to give room for big title
  old_par <- par(mar = c(5, 5, 6, 5))  # bottom, left, top, right
  
  vennDiagram(decidedTests,
              lwd = 1,
              circle.col = colors,
              names = contrasts,
              cex = 2.0,         # Increase set label text size
              cex.main = 1.5,    # Increase title size
              main = paste0("Venn diagram for various contrasts in ", tissueOfInterest, " tissue. Pval < 0.05")
  )
  
  par(old_par) # restore original margins
  
  dev.off()
  
  
}






#Venn with nominal p value
PlotVennDiagramsIleum <- function(fit, contrasts, geneInfo, tissueOfInterest, dirPath) {
  
  browser()
  
  # Nominal p-value filtering (0.05) without adjustment
  decidedTests <- ifelse(fit$p.value < 0.05, 
                         ifelse(fit$coef > 0, 1, -1), 
                         0)
  decidedTests <- as.matrix(decidedTests)
  colnames(decidedTests) <- colnames(fit$p.value)
  rownames(decidedTests) <- rownames(fit$p.value)
  
  # Venn counts
  vennCounts <- vennCounts(decidedTests)
  vennCounts <- as.data.frame(as.table(vennCounts))
  vennCounts <- tidyr::pivot_wider(vennCounts, names_from = Var2, values_from = Freq)
  vennCounts <- vennCounts[, -1]
  
  # Save Venn counts
  write.csv(vennCounts, paste0(dirPath, "VennCounts.csv"), row.names = FALSE)
  
  # Filter out zero counts and cases where all are zero
  vennCounts <- vennCounts[vennCounts$Counts != 0, ]
  # vennCounts <- vennCounts[!(vennCounts$HF_Surgery == 0 &
  #                              vennCounts$LF_Surgery == 0 &
  #                              vennCounts$Diet == 0 &
  #                              vennCounts$Surgery == 0), ]
  
  # Loop through Venn combinations and save gene lists
  for (i in 1:nrow(vennCounts)) {
    
    selectedContrasts <- colnames(vennCounts)[which(vennCounts[i, (1:4)] == 1)]
    nonSelectedContrasts <- setdiff(colnames(decidedTests), selectedContrasts)
    
    filteredGenes <- decidedTests[
      rowSums(abs(decidedTests[, selectedContrasts, drop = FALSE]) == 1) == length(selectedContrasts) &
        rowSums(abs(decidedTests[, nonSelectedContrasts, drop = FALSE]) == 1) == 0, 
    ]
    
    commonGenes <- rownames(filteredGenes)
    commonGeneSymbols <- geneInfo[geneInfo$ensembl_gene_id %in% commonGenes, ]$gene_name
    
    if (i == 1) {
      wb <- openxlsx::createWorkbook()
      firstCharacters <- substr(selectedContrasts, 1, 1)
      sheetName <- paste(firstCharacters, collapse = "_")
      openxlsx::addWorksheet(wb, sheetName)
      openxlsx::writeData(wb, sheetName, commonGeneSymbols, startCol = 1, startRow = 1)
      openxlsx::saveWorkbook(wb, paste0(dirPath, "VennGeneList.xlsx"), overwrite = TRUE)
    } else {
      wb <- openxlsx::loadWorkbook(paste0(dirPath, "VennGeneList.xlsx"))
      firstCharacters <- substr(selectedContrasts, 1, 1)
      sheetName <- paste(firstCharacters, collapse = "_")
      if (!(sheetName %in% names(wb))) {
        openxlsx::addWorksheet(wb, sheetName)
      }
      openxlsx::writeData(wb, sheetName, commonGeneSymbols, startCol = 1, startRow = 1)
      openxlsx::saveWorkbook(wb, paste0(dirPath, "VennGeneList.xlsx"), overwrite = TRUE)
    }
    
  }
  
  # Color palette
  baseColors <- c("red", "blue", "green")
  custom_palette <- colorRampPalette(baseColors)
  colors <- custom_palette(length(contrasts))
  
  # # Save Venn diagram image
  # jpeg(filename = paste0(dirPath, "Venn.jpeg"), width = 1000, height = 900, quality = 100)
  # vennDiagram(decidedTests,
  #             lwd = 1,
  #             circle.col = colors,
  #             names = contrasts,
  #             main = paste0("Venn diagram for various contrasts in ", tissueOfInterest, " tissue. Pval < 0.05")
  # )
  # dev.off()
  
  
  # Save Venn diagram image
  jpeg(filename = paste0(dirPath, "Venn.jpeg"), width = 1000, height = 800, quality = 100)
  
  
  # Increase top margin to give room for big title
  old_par <- par(mar = c(5, 5, 6, 5))  # bottom, left, top, right
  
  vennDiagram(decidedTests,
              lwd = 1,
              circle.col = colors,
              names = contrasts,
              cex = 2.0,         # Increase set label text size
              cex.main = 1.5,    # Increase title size
              main = paste0("Venn diagram for various contrasts in ", tissueOfInterest, " tissue. Pval < 0.05")
  )
  
  par(old_par) # restore original margins
  
  dev.off()
  
  
}









#Version1 - Works fine but names not clearly visible
# PlotVolcano <- function(contrasts,tissueOfInterest,dirPath){
#   
#   volcanoPlotList <- list()
#   
#   for (contrast in contrasts) {
#     allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast,rowNames = TRUE, colNames = TRUE)
#     #allGenes <- read.csv(paste0(dirPath,"/",contrast,"_AllGenes.csv"))
#     
#     sheetNames <- getSheetNames(paste0(dirPath,"SigGenes.xlsx"))
#     
#     
#     if(contrast %in% sheetNames){
#       
#       sigGenes <- read.xlsx(paste0(dirPath,"SigGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
#       
#       top10 <- head(sigGenes,n=10)
#       
#       #Volcano plot to visualize the results
#       volcanoPlot <- allGenes %>%
#         ggplot(aes(x = logFC,
#                    y = -log10(adj.P.Val),
#                    colour = isSignificant)) +
#         geom_vline(xintercept = 1, linetype = "dotted") +
#         geom_vline(xintercept = -1, linetype = "dotted") +
#         geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
#         geom_point(size = 1, alpha = 0.5) +
#         scale_colour_manual(values = c("grey", "red")) +
#         #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
#         geom_text_repel(data = top10 ,aes(label = GeneSymbol)) +
#         ggtitle(paste0(contrast,", Thresholds: FDR=0.05, logFC=1")) +
#         theme(plot.title = element_text(hjust = 0.5))
#       
#       volcanoPlotList <- append(volcanoPlotList,list(volcanoPlot))
#       
#     }else{
# 
#       
#       top10 <- as.data.frame(matrix(ncol = ncol(allGenes),nrow = 0))
#       p <- 0.001
#       while (TRUE) {
#         
#         allGenes <- allGenes %>%
#           dplyr::mutate(isPValSignificant = case_when(
#             P.Value < p & abs(logFC) > 1 ~ TRUE,
#             TRUE ~ FALSE
#           ))
#         
#         sigGenes <- allGenes %>%
#           dplyr::filter(isPValSignificant == TRUE)
#         
#         top10 <- head(sigGenes,n=10)
#         
#         if(nrow(top10) > 0) break
#         
#         p <- p + 0.001
#         
#       }
#       
#       #Volcano plot to visualize the results
#       volcanoPlot <- allGenes %>%
#         ggplot(aes(x = logFC,
#                    y = -log10(P.Value),
#                    colour = isPValSignificant)) +
#         geom_vline(xintercept = 1, linetype = "dotted") +
#         geom_vline(xintercept = -1, linetype = "dotted") +
#         geom_hline(yintercept = -log10(p), linetype = "dotted") +
#         geom_point(size = 1, alpha = 0.5) +
#         scale_colour_manual(values = c("grey", "red")) +
#         #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
#         geom_text_repel(data = top10 ,aes(label = GeneSymbol)) +
#         ggtitle(paste0(contrast,", Thresholds: Pval=",p,", logFC=1")) +
#         theme(plot.title = element_text(hjust = 0.5))
#       
#       volcanoPlotList <- append(volcanoPlotList,list(volcanoPlot))
#       
#     } 
#     
#     
#     
#     
#   }
#   
#   volcanoPlotGrid <- plot_grid(plotlist =  volcanoPlotList) + 
#     ggtitle(paste0("Volcano plots of different contrasts for ",tissueOfInterest," tissue"))+
#     theme(plot.title = element_text(hjust = 0.5, size = 20))
#   
#   ggsave(paste0(dirPath,"VolcanoPlot.jpeg"), volcanoPlotGrid, width = 18, height = 18,dpi = 300)
#   
# }

#Version2 - Corrected for unclear labels - Tested manually but not in loop
PlotVolcano <- function(contrasts,tissueOfInterest,dirPath){
  
  volcanoPlotList <- list()
  
  for (contrast in contrasts) {
    allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast,rowNames = TRUE, colNames = TRUE)
    #allGenes <- read.csv(paste0(dirPath,"/",contrast,"_AllGenes.csv"))
    
    sheetNames <- getSheetNames(paste0(dirPath,"SigGenes.xlsx"))
    
    
    if(contrast %in% sheetNames){
      
      sigGenes <- read.xlsx(paste0(dirPath,"SigGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
      
      top10 <- head(sigGenes,n=10)
      
      #Volcano plot to visualize the results
      volcanoPlot <- allGenes %>%
        ggplot(aes(x = logFC,
                   y = -log10(adj.P.Val),
                   colour = isSignificant)) +
        geom_vline(xintercept = 1, linetype = "dotted") +
        geom_vline(xintercept = -1, linetype = "dotted") +
        geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey", "red")) +
        #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
        geom_text_repel(data = top10 ,aes(label = GeneSymbol),size = 5) +
        ggtitle(paste0(contrast,", Thresholds: FDR<0.05, logFC>1")) +
        guides(color = guide_legend(title = "IsSignificant")) +
        theme(plot.title = element_text(hjust = 0.5,size = 15),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 10)
        )
      
      if(!contrast == contrasts[length(contrasts)]) volcanoPlot <- volcanoPlot + theme(legend.position = "none")
      
      volcanoPlotList <- append(volcanoPlotList,list(volcanoPlot))
      
    }else{
      
      
      top10 <- as.data.frame(matrix(ncol = ncol(allGenes),nrow = 0))
      p <- 0.001
      while (TRUE) {
        
        allGenes <- allGenes %>%
          dplyr::mutate(isPValSignificant = case_when(
            P.Value < p & abs(logFC) > 1 ~ TRUE,
            TRUE ~ FALSE
          ))
        
        sigGenes <- allGenes %>%
          dplyr::filter(isPValSignificant == TRUE)
        
        top10 <- head(sigGenes,n=10)
        
        if(nrow(top10) > 0) break
        
        p <- p + 0.001
        
      }
      
      #Volcano plot to visualize the results
      volcanoPlot <- allGenes %>%
        ggplot(aes(x = logFC,
                   y = -log10(P.Value),
                   colour = isPValSignificant)) +
        geom_vline(xintercept = 1, linetype = "dotted") +
        geom_vline(xintercept = -1, linetype = "dotted") +
        geom_hline(yintercept = -log10(p), linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey", "red")) +
        #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
        geom_text_repel(data = top10 ,aes(label = GeneSymbol), size = 5) +
        ggtitle(paste0(contrast,", Thresholds: Pval<",p,", logFC>1")) +
        guides(color = guide_legend(title = "IsSignificant")) +
        theme(plot.title = element_text(hjust = 0.5,size = 15),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 10))
      
      if(!contrast == contrasts[length(contrasts)]) volcanoPlot <- volcanoPlot + theme(legend.position = "none")
      
      volcanoPlotList <- append(volcanoPlotList,list(volcanoPlot))
      
    } 
    
    
    
    
  }
  
  volcanoPlotGrid <- plot_grid(plotlist =  volcanoPlotList,nrow = 1,ncol = 5) + 
    ggtitle(paste0("Volcano plots of different contrasts for ",tissueOfInterest," tissue"))+
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0(dirPath,"VolcanoPlot.jpeg"), volcanoPlotGrid, width = 35, height = 9,dpi = 600)
  
}

#Version2 - Corrected for unclear labels - Tested manually but not in loop
PlotVolcanoIleum <- function(contrasts,tissueOfInterest,dirPath){
  
  for (contrast in contrasts) {
    allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast,rowNames = TRUE, colNames = TRUE)
    #allGenes <- read.csv(paste0(dirPath,"/",contrast,"_AllGenes.csv"))
    
    sheetNames <- getSheetNames(paste0(dirPath,"SigGenes.xlsx"))
    
    
    if(contrast %in% sheetNames){
      
      sigGenes <- read.xlsx(paste0(dirPath,"SigGenes.xlsx"), sheet = contrast, rowNames = TRUE)
      
      top10 <- head(sigGenes,n=10)
      
      #Volcano plot to visualize the results
      volcanoPlot <- allGenes %>%
        ggplot(aes(x = logFC,
                   y = -log10(adj.P.Val),
                   colour = isSignificant)) +
        geom_vline(xintercept = 1, linetype = "dotted") +
        geom_vline(xintercept = -1, linetype = "dotted") +
        geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey", "red")) +
        #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
        geom_text_repel(data = top10 ,aes(label = GeneSymbol),size = 5) +
        ggtitle(paste0("Volcano plot\n",contrast,", Thresholds: FDR<0.05, logFC>1")) +
        guides(color = guide_legend(title = "IsSignificant")) +
        theme(plot.title = element_text(hjust = 0.5,size = 15),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 10)
        )
      
    }else{
      
      
      top10 <- as.data.frame(matrix(ncol = ncol(allGenes),nrow = 0))
      p <- 0.001
      while (TRUE) {
        
        allGenes <- allGenes %>%
          dplyr::mutate(isPValSignificant = case_when(
            P.Value < p & abs(logFC) > 1 ~ TRUE,
            TRUE ~ FALSE
          ))
        
        sigGenes <- allGenes %>%
          dplyr::filter(isPValSignificant == TRUE)
        
        top10 <- head(sigGenes,n=10)
        
        if(nrow(top10) > 0) break
        
        p <- p + 0.001
        
      }
      
      #Volcano plot to visualize the results
      volcanoPlot <- allGenes %>%
        ggplot(aes(x = logFC,
                   y = -log10(P.Value),
                   colour = isPValSignificant)) +
        geom_vline(xintercept = 1, linetype = "dotted") +
        geom_vline(xintercept = -1, linetype = "dotted") +
        geom_hline(yintercept = -log10(p), linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey", "red")) +
        #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
        geom_text_repel(data = top10 ,aes(label = GeneSymbol), size = 5) +
        ggtitle(paste0("Volcano plot\n",contrast,", Thresholds: Pval<",p,", logFC>1")) +
        guides(color = guide_legend(title = "IsSignificant")) +
        theme(plot.title = element_text(hjust = 0.5,size = 15),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 10))
      
    } 
    
    
    
    
  }
  
  
  
  ggsave(paste0(dirPath,"VolcanoPlot.jpeg"), volcanoPlot, width =15 , height = 10,dpi = 600)
  
}

PlotPValueHistogram <- function(contrasts,tissueOfInterest,dirPath){
  
  pValueHistogramList <- list()
  
  for (contrast in contrasts) {
    
    allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast,rowNames = TRUE, colNames = TRUE)
    
    pValueHistogramPlot <- ggplot(allGenes, aes(x = P.Value)) +
      geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.7) +
      labs(
        title = paste0(contrast," contrast"),
        x = "P-Value",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 25)
      )
    
    pValueHistogramList <- append(pValueHistogramList,list(pValueHistogramPlot))
    
  }
  
  
  pValueHistogramGrid <- plot_grid(plotlist =  pValueHistogramList,nrow = 5,ncol = 1) + 
    ggtitle(paste0("pValue histograms\n",tissueOfInterest," tissue"))+
    theme(plot.title = element_text(hjust = 0.5, size = 45,face = "bold"))
  
  ggsave(paste0(dirPath,"pValueHistogram.jpeg"), pValueHistogramGrid, width = 20, height = 40,dpi = 600)
  
}


PlotPValueHistogramIleum <- function(contrasts,tissueOfInterest,dirPath){
  
  
  for (contrast in contrasts) {
    
    allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast,rowNames = TRUE, colNames = TRUE)
    
    pValueHistogramPlot <- ggplot(allGenes, aes(x = P.Value)) +
      geom_histogram(binwidth = 0.05, color = "black", fill = "lightblue", alpha = 0.7) +
      labs(
        title = paste0("pValue histogram\n",contrast, " contrast."),
        x = "P-Value",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
        axis.title = element_text(size = 30),
        axis.text = element_text(size = 25)
      )
    
  }
  
  ggsave(paste0(dirPath,"pValueHistogram.jpeg"), pValueHistogramPlot, width = 20, height = 20,dpi = 600)
  
}


PlotMA <- function(contrasts,tissueOfInterest,dirPath){
  
  MAPlotList <- list()
  
  for (contrast in contrasts) {
    allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
    #allGenes <- read.csv(paste0(dirPath,"/",contrast,"_AllGenes.csv"))
    
    
    sheetNames <- getSheetNames(paste0(dirPath,"SigGenes.xlsx"))
    
    if(contrast %in% sheetNames){
      
      sigGenes <- read.xlsx(paste0(dirPath,"SigGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
      
      top10 <- head(sigGenes,n=10)
      
      #Volcano plot to visualize the results
      MAPlot <- allGenes %>%
        ggplot(aes(x = AveExpr,
                   y = logFC,
                   colour = isSignificant)) +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey", "red")) +
        geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
        ggtitle(paste0(contrast)) +
        theme(plot.title = element_text(hjust = 0.5))
      
      MAPlotList <- append(MAPlotList,list(MAPlot))
      
    }else{
      
      #Volcano plot to visualize the results
      MAPlot <- allGenes %>%
        ggplot(aes(x = AveExpr,
                   y = logFC,
                   colour = isSignificant)) +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey","red")) +
        ggtitle(paste0(contrast)) +
        theme(plot.title = element_text(hjust = 0.5))
      
      MAPlotList <- append(MAPlotList,list(MAPlot))
      
    } 
    
    
    
    
  }
  
  MAPlotGrid <- plot_grid(plotlist =  MAPlotList) + 
    ggtitle(paste0("MA plots of different contrasts for ",tissueOfInterest," tissue"))+
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0(dirPath,"MAPlot.jpeg"), MAPlotGrid, width = 18, height = 18,dpi = 300)
  
}

PlotMAIleum <- function(contrasts,tissueOfInterest,dirPath){
  
  for (contrast in contrasts) {
    allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
    #allGenes <- read.csv(paste0(dirPath,"/",contrast,"_AllGenes.csv"))
    
    
    sheetNames <- getSheetNames(paste0(dirPath,"SigGenes.xlsx"))
      
      sigGenes <- read.xlsx(paste0(dirPath,"SigGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
      
      top10 <- head(sigGenes,n=10)
      
      #Volcano plot to visualize the results
      MAPlot <- allGenes %>%
        ggplot(aes(x = AveExpr,
                   y = logFC,
                   colour = isSignificant)) +
        geom_hline(yintercept = 0, linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey", "red")) +
        #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
        geom_text_repel(data = top10 ,aes(label = GeneSymbol), size = 5) +
        ggtitle(paste0("MA plot\n",contrast,", Thresholds: FDR<0.05, logFC>1")) +
        theme(plot.title = element_text(hjust = 0.5))
    
    
    
    
  }
  
  ggsave(paste0(dirPath,"MAPlot.jpeg"), MAPlot, width = 9, height = 9,dpi = 600)
  
}



PlotHeatmaps <- function(ctsAll,contrasts,tissueOfInterest,dirPath){
  
  for (contrast in contrasts) {
    
    if(contrast == "LF_Surgery"){
      
      requiredColumns <- grep("^LF",colnames(ctsAll))
      cts <- ctsAll[,requiredColumns]
      
      sampleNames <- colnames(cts)
      
      lf_C_SamplesCount <- length(grep(paste0("^", "LF_C"), sampleNames))
      lf_S_SamplesCount <- length(grep(paste0("^", "LF_S"), sampleNames))
      
      sampleGroups <- factor(c(rep("Low fat Control",lf_C_SamplesCount),
                               rep("Low fat Surgery",lf_S_SamplesCount)
      ))
      
      groupColors <- c("Low fat Control" = "blue",
                       "Low fat Surgery" = "red"
      )
      
      desiredOrder <- c("LF_C", "LF_S")
      
    }else if(contrast == "HF_Surgery"){
      
      requiredColumns <- grep("^HF",colnames(ctsAll))
      cts <- ctsAll[,requiredColumns]
      
      sampleNames <- colnames(cts)
      
      hf_C_SamplesCount <- length(grep(paste0("^", "HF_C"), sampleNames))
      hf_S_SamplesCount <- length(grep(paste0("^", "HF_S"), sampleNames))
      
      sampleGroups <- factor(c(rep("High fat Control",hf_C_SamplesCount),
                               rep("High fat Surgery",hf_S_SamplesCount)
      ))
      
      groupColors <- c("High fat Control" = "blue",
                       "High fat Surgery" = "red"
      )
      
      desiredOrder <- c("HF_C", "HF_S")
      
    }else if(contrast == "Surgery") {
      
      
      cts <- ctsAll
      sampleNames <- colnames(cts)
      
      hf_C_SamplesCount <- length(grep(paste0("^", "HF_C"), sampleNames))
      hf_S_SamplesCount <- length(grep(paste0("^", "HF_S"), sampleNames))
      lf_C_SamplesCount <- length(grep(paste0("^", "LF_C"), sampleNames))
      lf_S_SamplesCount <- length(grep(paste0("^", "LF_S"), sampleNames))
      
      controlSamplesCount <- hf_C_SamplesCount + lf_C_SamplesCount
      surgerySamplesCount <- hf_S_SamplesCount + lf_S_SamplesCount
      
      sampleGroups <- factor(c(rep("Control",controlSamplesCount),
                               rep("Surgery",surgerySamplesCount)
      ))
      
      groupColors <- c("Control" = "blue",
                       "Surgery" = "red"
      )
      
      desiredOrder <- c("HF_C","LF_C", "HF_S","LF_S")
      
    }else if(contrast == "Diet"){
      
      cts <- ctsAll
      sampleNames <- colnames(cts)
      
      hf_C_SamplesCount <- length(grep(paste0("^", "HF_C"), sampleNames))
      hf_S_SamplesCount <- length(grep(paste0("^", "HF_S"), sampleNames))
      lf_C_SamplesCount <- length(grep(paste0("^", "LF_C"), sampleNames))
      lf_S_SamplesCount <- length(grep(paste0("^", "LF_S"), sampleNames))
      
      hfSamplesCount <- hf_C_SamplesCount + hf_S_SamplesCount
      lfSamplesCount <- lf_C_SamplesCount + lf_S_SamplesCount
      
      sampleGroups <- factor(c(rep("Highfat",hfSamplesCount),
                               rep("Lowfat",lfSamplesCount)
      ))
      
      groupColors <- c("Highfat" = "blue",
                       "Lowfat" = "red"
      )
      
      desiredOrder <- c("HF_C","HF_S", "LF_C","LF_S")
      
    }else {
      
      cts <- ctsAll
      sampleNames <- colnames(cts)
      
      hf_C_SamplesCount <- length(grep(paste0("^", "HF_C"), sampleNames))
      hf_S_SamplesCount <- length(grep(paste0("^", "HF_S"), sampleNames))
      lf_C_SamplesCount <- length(grep(paste0("^", "LF_C"), sampleNames))
      lf_S_SamplesCount <- length(grep(paste0("^", "LF_S"), sampleNames))
      
      sampleGroups <- factor(c(rep("High fat Control",hf_C_SamplesCount),
                               rep("High fat Surgery",hf_S_SamplesCount),
                               rep("Low fat Control",lf_C_SamplesCount),
                               rep("Low fat Surgery",lf_S_SamplesCount)
      ))
      
      groupColors <- c("High fat Control" = "blue",
                       "High fat Surgery" = "red",
                       "Low fat Control" = "orange",
                       "Low fat Surgery" = "green"
      )
      
      desiredOrder <- c("HF_C","HF_S", "LF_C","LF_S")
      
    }
    
    sheetNames <- getSheetNames(paste0(dirPath,"SigGenes.xlsx"))
    
    if(!contrast %in% sheetNames) next
    
    sigDEresults <- read.xlsx(paste0(dirPath,"SigGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
    
    sigDEresults <- head(sigDEresults[order(sigDEresults$adj.P.Val),],n=50)
    
    for (i in 1:nrow(sigDEresults)) {
      
      if(!nzchar(sigDEresults$GeneSymbol[i])) sigDEresults$GeneSymbol[i] <- sigDEresults$EnsemblID[i]
      
    }
    
    sigGenes <- sigDEresults$EnsemblID
    
    sigGeneExpressionData <- cts[rownames(cts) %in% sigGenes,]
    
    sampleNames <- colnames(sigGeneExpressionData)
    
    orderSamples <- unlist(lapply(desiredOrder, function(prefix) {
      sampleNames[grep(paste0("^", prefix), sampleNames)]
    }))
    
    sigGeneExpressionData <- sigGeneExpressionData[,orderSamples]
    
    rownames(sigGeneExpressionData) <- sigDEresults$GeneSymbol
    
    scaledSigExpData <- t(scale(t(sigGeneExpressionData)))
    
    topAnnotation <- HeatmapAnnotation(
      Group = sampleGroups,
      col = list(Group = groupColors)
    )
    
    jpeg(paste0(dirPath,contrast,"_TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          cluster_columns = FALSE,
                          top_annotation = topAnnotation,
                          column_title = paste0("Heatmap of significant genes (FDR<0.05) for ",contrast," contrast of ",tissueOfInterest," tissue"),
                          column_title_gp = gpar(fontsize = 20, fontface = "bold")
    )
    
    print(topHeatmap)
    
    dev.off()
    
    
    
  }
  
}


PlotHeatmapsIleum <- function(ctsAll,contrasts,tissueOfInterest,dirPath){

  for (contrast in contrasts) {
    
    cts <- ctsAll
    sampleNames <- colnames(cts)
    
    II_SamplesCount <- length(grep(paste0("_II_"), sampleNames))
    I_SamplesCount <- length(grep(paste0("_I_"), sampleNames))
    
    if(grepl("Control",contrast)){
      sampleGroups <- factor(c(rep("Ileum Interposition",II_SamplesCount),
                               rep("Ileum Control",I_SamplesCount)
      ))
      
      groupColors <- c("Ileum Interposition" = "blue",
                       "Ileum Control" = "red"
      )
    }else{
      sampleGroups <- factor(c(rep("Ileum Interposition",II_SamplesCount),
                               rep("Ileum Surgery",I_SamplesCount)
      ))
      
      groupColors <- c("Ileum Interposition" = "blue",
                       "Ileum Surgery" = "red"
      )
    }
    
    
    
    desiredOrder <- c("_II_","_I_")
    
    sheetNames <- getSheetNames(paste0(dirPath,"SigGenes.xlsx"))
    
    if(!contrast %in% sheetNames) next
    
    sigDEresults <- read.xlsx(paste0(dirPath,"SigGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
    
    sigDEresults <- head(sigDEresults[order(sigDEresults$adj.P.Val),],n=50)
    
    for (i in 1:nrow(sigDEresults)) {
      
      if(!nzchar(sigDEresults$GeneSymbol[i])) sigDEresults$GeneSymbol[i] <- sigDEresults$EnsemblID[i]
      
    }
    
    sigGenes <- sigDEresults$EnsemblID
    
    sigGeneExpressionData <- cts[rownames(cts) %in% sigGenes,]
    
    sampleNames <- colnames(sigGeneExpressionData)
    
    orderSamples <- unlist(lapply(desiredOrder, function(prefix) {
      sampleNames[grep(paste0(prefix), sampleNames)]
    }))
    
    sigGeneExpressionData <- sigGeneExpressionData[,orderSamples]
    
    rownames(sigGeneExpressionData) <- sigDEresults$GeneSymbol
    
    scaledSigExpData <- t(scale(t(sigGeneExpressionData)))
    
    topAnnotation <- HeatmapAnnotation(
      Group = sampleGroups,
      col = list(Group = groupColors),
      annotation_legend_param = list(
        Group = list(
          title_gp = gpar(fontsize = 16, fontface = "bold"),  # legend title
          labels_gp = gpar(fontsize = 14)                     # legend labels
        )
    ))
    
    jpeg(paste0(dirPath,contrast,"_TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          cluster_columns = FALSE,
                          top_annotation = topAnnotation,
                          column_title = paste0("Heatmap of significant genes (FDR<0.05) for ",contrast," contrast."),
                          column_title_gp = gpar(fontsize = 20, fontface = "bold"),
                          heatmap_legend_param = list(
                            title_gp = gpar(fontsize = 16, fontface = "bold"),  # legend title font
                            labels_gp = gpar(fontsize = 14)                     # legend labels font
                          )
    )
    
    print(topHeatmap)
    
    dev.off()
    
    
    
  }
  
}


SummarizeDEGenesFDR <- function(tissues,fdr,logFC){
  
  summaryDf <- data.frame(matrix(ncol = 5,nrow = 0))
  colnames(summaryDf) <- c("HF_Surgery","LF_Surgery","Diet","Surgery_Diet","Surgery")
  
  for (tissue in tissues) {
    
    dirPath <- paste0("data/LimmaResults/",tissue,"/")
    
    contrasts <- readLines(paste0(dirPath,"ContrastsUsed.txt"))
    
    for (i in 1:length(contrasts)) {
      
      contrasts[i] <- sub(" =.*","",contrasts[i])
      
    }
    
    tissueDf <- data.frame(matrix(ncol = length(contrasts),nrow = 1))
    colnames(tissueDf) <- contrasts
    tissueDf[1,] <- 0
    rownames(tissueDf) <- tissue
    
    for (contrast in contrasts) {
      
      allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
      
      allGenes <- allGenes %>%
        dplyr::mutate(isFDRSignificant = case_when(
          adj.P.Val < fdr & abs(logFC) > logFC ~ TRUE,
          TRUE ~ FALSE
        ))
      
      tissueDf[][[contrast]] <- nrow(allGenes[allGenes$isFDRSignificant==TRUE,])
      
    }
    
    summaryDf <- rbind(summaryDf,tissueDf)
    
  }
  
  summaryPlotDf <- summaryDf
  summaryPlotDf$Tissue <- rownames(summaryPlotDf)
  summaryPlotList <- list()
  
  for (contrast in contrasts){
    
    
    
    summaryPlot <- ggplot(summaryPlotDf,aes(x=Tissue,y=.data[[contrast]]))+
      geom_col() +
      labs(y="No of significant DE genes")+
      ggtitle(paste0("Contrast : ",contrast,", FDR cutoff: ",fdr,", logFC cutoff: ",logFC))+
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1,size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15)
      )
    
    
    summaryPlotList <- append(summaryPlotList,list(summaryPlot))
    
    
    
  }
  
  
  summaryPlotGrid <- plot_grid(plotlist =  summaryPlotList) + 
    ggtitle(paste0("Total significant DE genes of different contrasts Vs tissues"))+
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0("data/LimmaResults/GeneFDRSummary.jpeg"), summaryPlotGrid, width = 20, height = 18,dpi = 300)
  
  summaryPlotDf <- summaryPlotDf %>%
    dplyr::select(last_col(), everything())
  
  writeLines(paste0("FDR cutoff: ",fdr," and LogFC cutoff: ",logFC), "data/LimmaResults/GeneFDRSummary.csv")
  
  write.table(summaryPlotDf,file = "data/LimmaResults/GeneFDRSummary.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
  
}

SummarizeDEGenesPValue <- function(tissues,pValue,logFC){
  
  summaryDf <- data.frame(matrix(ncol = 5,nrow = 0))
  colnames(summaryDf) <- c("HF_Surgery","LF_Surgery","Diet","Surgery_Diet","Surgery")
  
  for (tissue in tissues) {
    
    dirPath <- paste0("data/LimmaResults/",tissue,"/")
    
    contrasts <- readLines(paste0(dirPath,"ContrastsUsed.txt"))
    
    for (i in 1:length(contrasts)) {
      
      contrasts[i] <- sub(" =.*","",contrasts[i])
      
    }
    
    tissueDf <- data.frame(matrix(ncol = length(contrasts),nrow = 1))
    colnames(tissueDf) <- contrasts
    tissueDf[1,] <- 0
    rownames(tissueDf) <- tissue
    
    for (contrast in contrasts) {
      
      allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
      
      allGenes <- allGenes %>%
        dplyr::mutate(isPVSignificant = case_when(
          P.Value < pValue & abs(logFC) > logFC ~ TRUE,
          TRUE ~ FALSE
        ))
      
      tissueDf[][[contrast]] <- nrow(allGenes[allGenes$isPVSignificant==TRUE,])
      
    }
    
    summaryDf <- rbind(summaryDf,tissueDf)
    
    
  }
  
  
  summaryPlotDf <- summaryDf
  summaryPlotDf$Tissue <- rownames(summaryPlotDf)
  summaryPlotList <- list()
  
  for (contrast in contrasts){
    
    
    
    summaryPlot <- ggplot(summaryPlotDf,aes(x=Tissue,y=.data[[contrast]]))+
      geom_col() +
      labs(y="No of significant DE genes")+
      ggtitle(paste0("Contrast : ",contrast,", Pvalue cutoff: ",pValue,", logFC cutoff: ",logFC))+
      theme(plot.title = element_text(hjust = 0.5, size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1,size = 15),
            axis.text.y = element_text(size = 15),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15)
      )
    
    
    summaryPlotList <- append(summaryPlotList,list(summaryPlot))
    
    
  }
  
  summaryPlotGrid <- plot_grid(plotlist =  summaryPlotList) + 
    ggtitle(paste0("Total significant DE genes of different contrasts Vs tissues"))+
    theme(plot.title = element_text(hjust = 0.5, size = 20))
  
  ggsave(paste0("data/LimmaResults/GenePValueSummary.jpeg"), summaryPlotGrid, width = 20, height = 18,dpi = 300)
  
  summaryPlotDf <- summaryPlotDf %>%
    dplyr::select(last_col(), everything())
  
  writeLines(paste0("Pvalue cutoff: ",pValue," and LogFC cutoff: ",logFC), "data/LimmaResults/GenePValueSummary.csv")
  
  write.table(summaryPlotDf,file = "data/LimmaResults/GenePValueSummary.csv",append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
  
}

SummarizePathwayFDR <- function(tissues,fdr,pathwayNames=c()){
  
  for (pathwayName in pathwayNames) {
    
    summaryDf <- data.frame(matrix(ncol = 5,nrow = 0))
    colnames(summaryDf) <- c("HF_Surgery","LF_Surgery","Diet","Surgery_Diet","Surgery")
    
    for (tissue in tissues) {
      
      readDirPath <- paste0("data/LimmaResults/",tissue,"/")
      writeDirPath <- paste0("data/FGSEAResults/",tissue,"/")
      
      contrasts <- readLines(paste0(readDirPath,"ContrastsUsed.txt"))
      
      for (i in 1:length(contrasts)) {
        
        contrasts[i] <- sub(" =.*","",contrasts[i])
        
      }
      
      tissueDf <- data.frame(matrix(ncol = length(contrasts),nrow = 1))
      colnames(tissueDf) <- contrasts
      tissueDf[1,] <- 0
      rownames(tissueDf) <- tissue
      
      for (contrast in contrasts) {
        
    
        
        sheetName <- paste0(pathwayName,"_",contrast)
        allPathways <- read.xlsx(paste0(writeDirPath,"AllPathways.xlsx"), sheet = sheetName, colNames = TRUE)
        
        allPathways <- allPathways %>%
          dplyr::mutate(isFDRSignificant = case_when(
            padj < fdr ~ TRUE,
            TRUE ~ FALSE
          ))
        
        tissueDf[][[contrast]] <- nrow(allPathways[allPathways$isFDRSignificant==TRUE,])
        
      }
      
      summaryDf <- rbind(summaryDf,tissueDf)
      
    }
    
    summaryPlotDf <- summaryDf
    summaryPlotDf$Tissue <- rownames(summaryPlotDf)
    summaryPlotList <- list()
    
    for (contrast in contrasts){
      
      
      
      summaryPlot <- ggplot(summaryPlotDf,aes(x=Tissue,y=.data[[contrast]]))+
        geom_col() +
        labs(y="No of significant pathways")+
        ggtitle(paste0("Contrast : ",contrast,", FDR cutoff: ",fdr))+
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1,size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15)
        )
      
      
      summaryPlotList <- append(summaryPlotList,list(summaryPlot))
      
      
      
    }
    
    
    summaryPlotGrid <- plot_grid(plotlist =  summaryPlotList) + 
      ggtitle(paste0("Total significant pathways of different contrasts Vs tissues"))+
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    
    ggsave(paste0("data/FGSEAResults/",pathwayName,"_PathwayFDRSummary.jpeg"), summaryPlotGrid, width = 20, height = 18,dpi = 300)
    
    summaryPlotDf <- summaryPlotDf %>%
      dplyr::select(last_col(), everything())
    
    writeLines(paste0("FDR cutoff: ",fdr), paste0("data/FGSEAResults/",pathwayName,"_PathwayFDRSummary.csv"))
    
    write.table(summaryPlotDf,file = paste0("data/FGSEAResults/",pathwayName,"_PathwayFDRSummary.csv"),append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
    
    
    
  }
  
}

SummarizePathwayPValue <- function(tissues,pValue,pathwayNames=c()){
  
  for (pathwayName in pathwayNames) {
    
    summaryDf <- data.frame(matrix(ncol = 5,nrow = 0))
    colnames(summaryDf) <- c("HF_Surgery","LF_Surgery","Diet","Surgery_Diet","Surgery")
    
    for (tissue in tissues) {
      
      readDirPath <- paste0("data/LimmaResults/",tissue,"/")
      writeDirPath <- paste0("data/FGSEAResults/",tissue,"/")
      
      contrasts <- readLines(paste0(readDirPath,"ContrastsUsed.txt"))
      
      for (i in 1:length(contrasts)) {
        
        contrasts[i] <- sub(" =.*","",contrasts[i])
        
      }
      
      tissueDf <- data.frame(matrix(ncol = length(contrasts),nrow = 1))
      colnames(tissueDf) <- contrasts
      tissueDf[1,] <- 0
      rownames(tissueDf) <- tissue
      
      for (contrast in contrasts) {
        
        sheetName <- paste0(pathwayName,"_",contrast)
        allPathways <- read.xlsx(paste0(writeDirPath,"AllPathways.xlsx"), sheet = sheetName, colNames = TRUE)
        
        allPathways <- allPathways %>%
          dplyr::mutate(isFDRSignificant = case_when(
            pval < pValue ~ TRUE,
            TRUE ~ FALSE
          ))
        
        tissueDf[][[contrast]] <- nrow(allPathways[allPathways$isFDRSignificant==TRUE,])
        
      }
      
      summaryDf <- rbind(summaryDf,tissueDf)
      
    }
    
    summaryPlotDf <- summaryDf
    summaryPlotDf$Tissue <- rownames(summaryPlotDf)
    summaryPlotList <- list()
    
    for (contrast in contrasts){
      
      
      
      summaryPlot <- ggplot(summaryPlotDf,aes(x=Tissue,y=.data[[contrast]]))+
        geom_col() +
        labs(y="No of significant pathways")+
        ggtitle(paste0("Contrast : ",contrast,", Pvalue cutoff: ",pValue))+
        theme(plot.title = element_text(hjust = 0.5, size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1,size = 15),
              axis.text.y = element_text(size = 15),
              axis.title.x = element_text(size = 15),
              axis.title.y = element_text(size = 15))
      
      
      summaryPlotList <- append(summaryPlotList,list(summaryPlot))
      
      
      
    }
    
    
    summaryPlotGrid <- plot_grid(plotlist =  summaryPlotList) + 
      ggtitle(paste0("Total significant pathways of different contrasts Vs tissues"))+
      theme(plot.title = element_text(hjust = 0.5, size = 20))
    
    ggsave(paste0("data/FGSEAResults/",pathwayName,"_PathwayPValueSummary.jpeg"), summaryPlotGrid, width = 20, height = 18,dpi = 300)
    
    summaryPlotDf <- summaryPlotDf %>%
      dplyr::select(last_col(), everything())
    
    writeLines(paste0("PValue cutoff: ",pValue), paste0("data/FGSEAResults/",pathwayName,"_PathwayPValueSummary.csv"))
    
    write.table(summaryPlotDf,file = paste0("data/FGSEAResults/",pathwayName,"_PathwayPValueSummary.csv"),append = TRUE,sep = ",",row.names = FALSE,col.names = TRUE)
    
    
    
  }
  
}

PlotBoxPlotsForTOPGenes <- function(tissueOfInterest,contrasts,noOfGenes,testStatistic,logFC,sampleInfo){
  
  dirPath <- "data/LimmaResults/"
  
  for (tissue in tissueOfInterest) {
    
    dirPath <- "data/LimmaResults/"
    dirPath <- paste0(dirPath,tissue,"/")
    
    for (contrast in contrasts) {
      
      allGenes <- read.xlsx(paste0(dirPath,"AllGenes.xlsx"), sheet = contrast,rowNames = TRUE, colNames = TRUE)

      if(testStatistic == "pvalue"){
        allGenes <- allGenes[abs(allGenes$logFC) > logFC,]
        allGenes <- allGenes[order(allGenes$P.Value), ]
      }else if(testStatistic == "fdr"){
        allGenes <- allGenes[abs(allGenes$logFC) > logFC,]
        allGenes <- allGenes[order(allGenes$adj.P.Val), ]
      }
      
      allGenes <- allGenes[1:noOfGenes,]
      
      topGenesEnsemblIds <- allGenes$EnsemblID
      topGenesSymbols <- allGenes$GeneSymbol
      
      cts <- read.csv(paste0(dirPath,"logCPMCounts.csv"),row.names = 1)
      
      filteresCts <- cts[rownames(cts) %in% topGenesEnsemblIds,]
      
      filteresCts <- filteresCts[match(topGenesEnsemblIds, rownames(filteresCts)), ]
      
      rownames(filteresCts) <- topGenesSymbols
      
      ctsLong <- pivot_longer(filteresCts,cols = everything(),names_to = "SampleID",values_to = "Counts")
      
      ctsLong$Gene <- rep(rownames(filteresCts), each = ncol(filteresCts))
      
      
      if(tissue != "Small Bowel Mucosa"){
        ctsLong <- ctsLong %>%
          dplyr::mutate(Group = sampleInfo$samplename[match(SampleID,sampleInfo$sampleid)])
      }else{
        ctsLong <- ctsLong %>%
          dplyr::mutate(Group = sampleInfo$samplenamegeneral[match(SampleID,sampleInfo$sampleid)])
      }
      
      
      boxPlot <- ggplot(ctsLong,aes(x=Gene,y=Counts,fill = Group))+
        geom_boxplot()+
        labs(x="Gene",
             y="logCPM",
             title = paste0("Box plot of top ",noOfGenes," genes in ",contrast," contrast of  ",tissue," tissue"), 
             caption = paste0("Genes ordered based on : ", testStatistic, " and logFC threshold: ",logFC))+
        theme(plot.title = element_text(hjust = 0.5, size = 15), 
              plot.caption = element_text(hjust = 0.5, size = 15),
              axis.title = element_text(size = 15),
              axis.text.x = element_text(size = 15),
              axis.text.y = element_text(size = 10),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 15)
              )
      
      ggsave(paste0(dirPath,contrast,"_BoxPlot.jpeg"),plot = boxPlot,height = 10, width = 10,dpi = 600)
      
    }
    
    
  }
  
}

AppendTissueNamesToAllFiles <- function(tissueOfInterest){
  
  dgeDir <- "data/LimmaResults/"
  fgseaDir <- "data/FGSEAResults/"
  
  #Append to DGE files
  for(tissue in tissueOfInterest){
    
    filesDir <- paste0(dgeDir,tissue)
    
    files <- list.files(filesDir,full.names = TRUE)
    
    for (file in files) {
      
      # Get the file name
      file_name <- basename(file)
      
      # Create the new file name by prepending the folder name
      new_name <- paste0(filesDir, "/", tissue, "_", file_name)
      
      # Rename the file
      file.rename(file, new_name)
      
      
    }
    
  }
  
  #Append to FGSEA files
  for(tissue in tissueOfInterest){
    
    filesDir <- paste0(fgseaDir,tissue)
    
    files <- list.files(filesDir,full.names = TRUE, recursive = TRUE)
    
    for (file in files) {
      
      # Get the file name and directory path
      file_name <- basename(file)
      file_dir <- dirname(file)
      
      # Create the new file name by prepending the tissue name
      new_name <- paste0(file_dir, "/", tissue, "_", file_name)
      
      # Rename the file
      file.rename(file, new_name)
      
      
    }
    
  }
  
  
}

GenerateHeatmapForSummaryValues <- function(inputFiles,OutputFiles){
  
  for (i in 1:length(inputFiles)) {
    
    fDRSummary <- read.csv(inputFiles[i],header = FALSE)
    threshold <- fDRSummary[1,1]
    colnames(fDRSummary) <- fDRSummary[2,]
    fDRSummary <- fDRSummary[-(1:2),]
    rownames(fDRSummary) <- fDRSummary$Tissue
    rowNames <- rownames(fDRSummary)
    fDRSummary <- fDRSummary[,-1]
    fDRSummary <- as.matrix(data.frame(lapply(fDRSummary, as.numeric)))
    rownames(fDRSummary) <- rowNames
    
    jpeg(filename = OutputFiles[i],quality = 100)
    
    
    heatMap <- Heatmap(fDRSummary,
            row_names_side = "left",
            column_title = paste0("Heatmap of significant pathways, ",threshold),
            column_names_side = "top",
            col = colorRamp2(c(min(fDRSummary),max(fDRSummary)),c("white","red")),
            heatmap_legend_param = list(title = "Number of Pathways"),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(round(fDRSummary[i, j], 2), x, y, gp = gpar(fontsize = 10))
            },
            cluster_rows = TRUE,
            cluster_columns = FALSE,
            
    )
    
    print(heatMap)
    
    
    
    dev.off()
    
    
    
  }
  
  
}

GenerateBarPlotsForPathways <- function(tissueOfInterest,pathwaysListNames,contrasts){
  
  
  for (tissue in tissueOfInterest) {
    
    #Use this before appending files
    wb <- loadWorkbook(paste0("data/FGSEAResults/",tissue,"/AllPathways.xlsx"))
    
    #Use this after appending files
    #wb <- loadWorkbook(paste0("data/FGSEAResults/",tissue,"/",tissue,".xlsx"))
    
    for (pathwaysList in pathwaysListNames) {
      
      
      
      for (contrast in contrasts) {
        
        sheetName <- paste0(pathwaysList,"_",contrast)
        
        sheet <- read.xlsx(wb,sheet = sheetName)
        
        sheet <- sheet[sheet$padj<0.05,]
        
        outputFile <- paste0("data/FGSEAResults/",tissue,"/",pathwaysList,"_Plots/",pathwaysList,"_",contrast,"_BarPlot.jpeg")
        
        plot <- ggplot(sheet, aes(x = reorder(pathway,NES), y = NES, fill = padj)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_gradient(low = "blue", high = "red", name = "padj",transform = "reverse") +
          coord_flip() + 
          labs(title = paste0("Up- and Down-Regulated Pathways, FDR < 0.05 \n",tissue," tissue, ",contrast," contrast"), x = "Pathway", y = "NES") +
          theme_minimal() +
          theme(legend.position = "right",
                plot.title = element_text(hjust = 0.5,size = 15),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                axis.text.y = element_text(size = 15,colour = "black",),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 15))
          
        
        ggsave(filename = outputFile,height = 10, width = 12,  dpi = 600)
        
      }
      
    }
    
  }
  
  
  
}

GenerateBarPlotsForPathwaysIleum <- function(tissueOfInterest,pathwaysListNames,contrasts){
  
  
  for (tissue in tissueOfInterest) {
    
    #Use this before appending files
    wb <- loadWorkbook(paste0("data/FGSEAResults/",tissue,"/AllPathways.xlsx"))
    
    #Use this after appending files
    #wb <- loadWorkbook(paste0("data/FGSEAResults/",tissue,"/",tissue,".xlsx"))
    
    for (pathwaysList in pathwaysListNames) {
      
      
      
      for (contrast in contrasts) {
        
        sheetName <- paste0(pathwaysList)
        
        sheet <- read.xlsx(wb,sheet = sheetName)
        
        sheet <- sheet[sheet$padj<0.05,]
        
        outputFile <- paste0("data/FGSEAResults/",tissue,"/",pathwaysList,"_Plots/",pathwaysList,"_",contrast,"_BarPlot.jpeg")
        
        plot <- ggplot(sheet, aes(x = reorder(pathway,NES), y = NES, fill = padj)) +
          geom_bar(stat = "identity", position = "stack") +
          scale_fill_gradient(low = "blue", high = "red", name = "padj",transform = "reverse") +
          coord_flip() + 
          labs(title = paste0("Up- and Down-Regulated Pathways, FDR < 0.05 \n",contrast," contrast"), x = "Pathway", y = "NES") +
          theme_minimal() +
          theme(legend.position = "right",
                plot.title = element_text(hjust = 0.5,size = 15),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                axis.text.y = element_text(size = 15,colour = "black",),
                legend.text = element_text(size = 15),
                legend.title = element_text(size = 15))
        
        
        ggsave(filename = outputFile,height = 10, width = 12,  dpi = 600)
        
      }
      
    }
    
  }
  
  
  
}

GenerateSignificanceHeatmapForPathways <- function(tissueOfInterest,pathwaysListNames){
  
  
  for (tissue in tissueOfInterest) {
    
    for (pathwaysList in pathwaysListNames) {
      
      sigMatrix <- read.csv(paste0("data/FGSEAResults/",tissue,"/",pathwaysList,"_Plots/",pathwaysList,"_SignificanceMatrix.csv"),row.names = 1)
      sigMatrix <- sigMatrix[!apply(sigMatrix, 1, function(row) all(row == 0)), ]
      sigMatrix <- as.matrix(sigMatrix)
      
      colFun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
      
      
      hm <- Heatmap(sigMatrix,
                    column_names_side = "top",
                    col = colFun,
                    border = TRUE,
                    row_names_gp = gpar(fontsize = 10),
                    heatmap_width = unit(17,"cm"),
                    heatmap_height = unit(22,"cm"),
                    show_heatmap_legend = FALSE,
                    rect_gp = gpar(col = "black", lwd = 1),
                    column_title =  paste0("Heatmap of significant pathways for ",tissue," tissue"),
                    column_title_gp = gpar(fontsize = 10, fontface = "bold"))
      
      lgd = Legend(
        labels = c("Up", "NotSig", "Down"),
        title = "Significance",
        legend_gp = gpar(fill = c("red", "white", "blue")),
        border = TRUE,
        direction = "horizontal", 
        title_gp = gpar(fontsize = 14),
        labels_gp = gpar(fontsize = 12)
      )
      
      outputFile <- paste0("data/FGSEAResults/",tissue,"/",pathwaysList,"_Plots/",pathwaysList,"_SignificanceHeatmap.jpeg")
      
      jpeg(filename = outputFile,quality = 100,height = 25,width = 32,units = "cm",res = 300)
      
      draw(hm,annotation_legend_list = list(lgd),
           annotation_legend_side = "left"
      )
      
      dev.off()
      
    }
    
  }
  
  
}

ConvertVennGeneListToSM <- function(tissueOfInterest){
  
  
  for (tissue in tissueOfInterest) {
    
    dataPath <- paste0("data/LimmaResults/",tissue,"/VennGeneList.xlsx")
    
    wb <- loadWorkbook(dataPath)
    
    sheetNames <- names(wb)
    
    genes <- c()
    
    for (sheet in sheetNames) {
      
      geneData <- read.xlsx(wb,sheet = sheet,colNames = FALSE)
      genes <- c(genes,geneData$X1)
      
    }
    
    genes <- unique(genes)
    
    vennGeneListMatrix <- as.data.frame(matrix(ncol = 4,nrow = length(genes),data = 0))
    colnames(vennGeneListMatrix) <- c("HF_Surgery","LF_Surgery","Surgery","Diet")
    rownames(vennGeneListMatrix) <- genes
    
    for (sheet in sheetNames) {
      
      combinations <- strsplit(sheet,"_")
      combinations <- combinations[[1]]
      
      for (i in 1:length(combinations)) {
        
        if(combinations[i]=="H") combinations[i] <- "HF_Surgery"
        if(combinations[i]=="L") combinations[i] <- "LF_Surgery"
        if(combinations[i]=="S") combinations[i] <- "Surgery"
        if(combinations[i]=="D") combinations[i] <- "Diet"
        
      }
      
      geneData <- read.xlsx(wb,sheet = sheet,colNames = FALSE)
      geneData <- geneData$X1
      
      for (gene in geneData) {
        
        
        vennGeneListMatrix[gene,colnames(vennGeneListMatrix) %in% combinations] <- 1
        
      }
      
    }
    
    writePath <- paste0("data/LimmaResults/",tissue,"/VennGeneListMatrix.csv")
    
    
    write.csv(vennGeneListMatrix,file = writePath,row.names = TRUE,col.names = TRUE)
    
    
  }
  
  
  
}

AppendRelevantFilesToAllGenes <- function(tissueOfInterest){
  
  
  for (tissue in tissueOfInterest) {
    
    dataPath <- paste0("data/LimmaResults/",tissue,"/")
    allGenes <- paste0(dataPath,"AllGenes.xlsx")
    logCPM <- paste0(dataPath,"logCPMCounts.csv")
    vennCounts <- paste0(dataPath,"VennCounts.csv")
    vennGeneListMatrix <- paste0(dataPath,"VennGeneListMatrix.csv")
    
    allGenes <- loadWorkbook(allGenes)
    logCPM <- read.csv(logCPM,row.names = 1)
    vennCounts <- read.csv(vennCounts)
    vennGeneListMatrix <- read.csv(vennGeneListMatrix,row.names = 1)
    
    addWorksheet(allGenes,"logCPM")
    writeData(allGenes,"logCPM",logCPM,startCol = 1,startRow = 1,rowNames = TRUE)
    
    addWorksheet(allGenes,"VennCounts")
    writeData(allGenes,"VennCounts",vennCounts,startCol = 1,startRow = 1)
    
    addWorksheet(allGenes,"VennGeneListMatrix")
    writeData(allGenes,"VennGeneListMatrix",vennGeneListMatrix,startCol = 1,startRow = 1,rowNames = TRUE)
    
    addWorksheet(allGenes,"Notes")
    
    saveWorkbook(allGenes,file = paste0("data/LimmaResults/",tissue,"/",tissue,".xlsx"),overwrite = TRUE)
    
  }
  
  
}

AppendRelevantFilesToAllPathways <- function(tissueOfInterest,pathwaysListNames){
  
  
  for (tissue in tissueOfInterest) {
    
    
    allPathways <- loadWorkbook(paste0("data/FGSEAResults/",tissue,"/AllPathways.xlsx"))
    
    
    for (pathwayList in pathwaysListNames) {
      
      pathwaySignificanceMatrix <- read.csv(paste0("data/FGSEAResults/",tissue,"/",pathwayList,"_Plots/",pathwayList,"_SignificanceMatrix.csv"),row.names = 1)
      
      addWorksheet(allPathways,sheetName = paste0(pathwayList,"_SignificanceMatrix"))
      writeData(allPathways,paste0(pathwayList,"_SignificanceMatrix"),pathwaySignificanceMatrix, startCol = 1, startRow = 1,rowNames = TRUE)
      
      
    }
    
    addWorksheet(allPathways,sheetName = "Notes")
    
    saveWorkbook(allPathways,file = paste0("data/FGSEAResults/",tissue,"/",tissue,".xlsx"),overwrite = TRUE )
    
    
  }
  
  
}

AddNotesToGeneFiles <- function(tissueOfInterest){
  
  GeneNotes <- c(
    "HF_Surgery" = "Differential Gene Expression Results for Highfat Control Vs Highfat Surgery samples within the respective tissue. The first column in the dataframe is a combination of EnsemblID and Genesymbol for a respective gene and it can be used as a rowname when reading into R or Python.",
    "LF_Surgery" = "Differential Gene Expression Results for Lowfat Control Vs Lowfat Surgery samples within the respective tissue. The first column in the dataframe is a combination of EnsemblID and Genesymbol for a respective gene and it can be used as a rowname when reading into R or Python.",
    "Surgery" = "Differential Gene Expression Results for all Control samples Vs Surgery samples, irrespective of diet, within the respective tissue. The first column in the dataframe is a combination of EnsemblID and Genesymbol for a respective gene and it can be used as a rowname when reading into R or Python.",
    "Diet" = "Differential Gene Expression Results for all Highfat samples Vs Lowfat samples, irrespective of surgical treatment, within the respective tissue. The first column in the dataframe is a combination of EnsemblID and Genesymbol for a respective gene and it can be used as a rowname when reading into R or Python.",
    "Surgery_Diet" = "Differential Gene Expression Results for interaction effect of diet and surgical treatment among all samples within the respective tissue. The first column in the dataframe is a combination of EnsemblID and Genesymbol for a respective gene and it can be used as a rowname when reading into R or Python.",
    "logCPM" = "Log normalizes CPM counts data used for the Differential Gene Expression Analysis",
    "VennCounts" = "The table shows the number of common genes that were differentially expressed across contrasts (except interaction contrast Surgery_Diet)",
    "VennGeneListMatrix" = "The matrix shows the actual genes that were significantly differentially expressed across contrasts ((except interaction contrast Surgery_Diet)). 1 indicates that a respective gene is significantly differentially expressed(p.adj < 0.05) in that particular contrast analysis "
  )
  
  for (tissue in tissueOfInterest) {
    
    wb <- loadWorkbook(paste0("data/LimmaResults/",tissue,"/",tissue,".xlsx"))
    
    sheetNames <- names(wb)
    
    sheetNames <- sheetNames[sheetNames!="Notes"]
    
    notesDf <- as.data.frame(matrix(ncol = 2,nrow = length(sheetNames),data = NA))
    colnames(notesDf) <- c("SheetName","Description")
    
    for (i in 1:length(sheetNames)) {
      
      notesDf$SheetName[i] <- sheetNames[i]
      notesDf$Description[i] <- GeneNotes[sheetNames[i]]
      
    }
    
    writeData(wb,"Notes",notesDf,startCol = 1, startRow = 1)
    
    saveWorkbook(wb,file = paste0("data/LimmaResults/",tissue,"/",tissue,".xlsx"),overwrite = TRUE)
    
  }
  
}

AddNotesToPathwayFiles <- function(tissueOfInterest){
  
  pathwayNotes <- c(
    "HallMark_HF_Surgery" = "Geneset Enrichment Analysis Results for HF_Surgery contrast using HallMark pathways",
    "HallMark_LF_Surgery" = "Geneset Enrichment Analysis Results for LF_Surgery contrast using HallMark pathways",
    "HallMark_Surgery" = "Geneset Enrichment Analysis Results for Surgery contrast using HallMark pathways",
    "HallMark_Diet" = "Geneset Enrichment Analysis Results for Diet contrast using HallMark pathways",
    "HallMark_Surgery_Diet" = "Geneset Enrichment Analysis Results for Surgery_Diet contrast using HallMark pathways",
    "Wiki_HF_Surgery" = "Geneset Enrichment Analysis Results for HF_Surgery contrast using Wiki pathways",
    "Wiki_LF_Surgery" = "Geneset Enrichment Analysis Results for LF_Surgery contrast using Wiki pathways",
    "Wiki_Surgery" = "Geneset Enrichment Analysis Results for Surgery contrast using Wiki pathways",
    "Wiki_Diet" = "Geneset Enrichment Analysis Results for Diet contrast using Wiki pathways",
    "Wiki_Surgery_Diet" = "Geneset Enrichment Analysis Results for Surgery_Diet contrast using Wiki pathways",
    "HallMark_SignificanceMatrix" = "The matrix shows up regulated, down regulated and not significant HallMark pathways across all contrasts within the respective tissue",
    "Wiki_SignificanceMatrix" = "The matrix shows up regulated, down regulated and not significant Wiki pathways across all contrasts within the respective tissue"
  )
  
  for (tissue in tissueOfInterest) {
    
    wb <- loadWorkbook(paste0("data/FGSEAResults/",tissue,"/",tissue,".xlsx"))
    
    sheetNames <- names(wb)
    
    sheetNames <- sheetNames[sheetNames!="Notes"]
    
    notesDf <- as.data.frame(matrix(ncol = 2,nrow = length(sheetNames),data = NA))
    colnames(notesDf) <- c("SheetName","Description")
    
    for (i in 1:length(sheetNames)) {
      
      notesDf$SheetName[i] <- sheetNames[i]
      notesDf$Description[i] <- pathwayNotes[sheetNames[i]]
      
    }
    
    writeData(wb,"Notes",notesDf,startCol = 1, startRow = 1)
    
    saveWorkbook(wb,file = paste0("data/FGSEAResults/",tissue,"/",tissue,".xlsx"),overwrite = TRUE)
    
  }
  
}

DeleteUnnecessaryGeneFiles <- function(tissueOfInterest){

  
  dataPath <- "data/Share/LimmaResults/"
  
  geneFilesToDeleteWithInTissue <- c("AllGenes.xlsx",
                                     "ContrastsUsed.txt",
                                     "DensityPlot.jpeg",
                                     "LibrarySizePlot.jpeg",
                                     "logCPMCounts.csv",
                                     "MAPlot.jpeg",
                                     "SigGenes.xlsx",
                                     "TrueNullCol.jpeg",
                                     "VennCounts.csv",
                                     "VennGeneList.xlsx",
                                     "VennGeneListMatrix.csv")
  
  for (tissue in tissueOfInterest) {
    
    geneFilesToDeleteWithInTissueFullPaths <- paste0(dataPath,tissue,"/",geneFilesToDeleteWithInTissue)
    
    file.remove(geneFilesToDeleteWithInTissueFullPaths)
    
  }
  
  
  
}

DeleteUnnecessaryPathwayFiles <- function(tissueOfInterest){
  
  dataPath <- "data/Share/FGSEAResults/"
  
  pathwayFilesToDeleteWithInTissue <- c("MainPathways.xlsx","AllPathways.xlsx")
  
  pathwayFilesToDeleteWithInTissuePlots <- c("MainPathways.jpeg",
                                                    "TopPathways.jpeg",
                                                    "UpVenn.jpeg",
                                                    "DownVenn.jpeg",
                                                    "UpVennPathways.xlsx",
                                                    "DownVennPathways.xlsx",
                                                    "SignificanceMatrix.csv")
  
  pathwaysListNames <- c("HallMark","Wiki")
  contrasts <- c("HF_Surgery","LF_Surgery","Surgery","Diet","Surgery_Diet")
  
  for (tissue in tissueOfInterest) {
    
    pathwayFilesToDeleteWithInTissueFullPaths <- paste0(dataPath,tissue,"/",pathwayFilesToDeleteWithInTissue)
    
    file.remove(pathwayFilesToDeleteWithInTissueFullPaths)
    
    for (pathwayList in pathwaysListNames) {
      
      pathwayFilesToDeleteWithInTissuePlotsFullPaths <- paste0(dataPath,tissue,"/",pathwayList,"_Plots/",pathwayList,"_",pathwayFilesToDeleteWithInTissuePlots[3:7])
      
      file.remove(pathwayFilesToDeleteWithInTissuePlotsFullPaths)
      
      for (contrast in contrasts) {
        
        pathwayFilesToDeleteWithInTissuePlotsFullPaths <- paste0(dataPath,tissue,"/",pathwayList,"_Plots/",pathwayList,"_",contrast,"_",pathwayFilesToDeleteWithInTissuePlots[1:2])
        file.remove(pathwayFilesToDeleteWithInTissuePlotsFullPaths)
        
      }
      
    }
    
  }
  
  
  
}

DoComprehensiveGeneAbundanceAnalysis <- function(tissueOfInterest){
  
  geneList <- c()
  dataDirs <- c()
  for(tissue in tissueOfInterest){
    
    dataPath <- paste0("data/LimmaResults/",tissue,"/VennGeneListMatrix.csv")
    
    dataDirs <- c(dataDirs,dataPath)
    
    genes <- read.csv(file = dataPath,row.names = 1)
    genes <- rownames(genes)
    geneList <- c(geneList,genes)
  }
  
  geneList <- unique(geneList)
  
  contrasts <- c("HF_Surgery","LF_Surgery","Surgery","Diet")
  
  comColNames <- c()
  for (tissue in tissueOfInterest) {
    
    comColNames <- c(comColNames,paste0(tissue,"_",contrasts))
    
  }
  
  comGeneDf <- as.data.frame(matrix(ncol = length(comColNames),nrow = length(geneList),data = 0))
  colnames(comGeneDf) <- comColNames
  rownames(comGeneDf) <- geneList
  
  
  
  for (i in 1:length(tissueOfInterest)) {
    
    geneListMatrix <- read.csv(file = dataDirs[i],row.names = 1)
    matrixColumns <- colnames(geneListMatrix)
    
    for (col in matrixColumns) {
      
      subGeneListMatrix <- geneListMatrix[geneListMatrix[,col]==1,]
      sigGenes <- rownames(subGeneListMatrix)
      comGeneDf[sigGenes,paste0(tissueOfInterest[i],"_",col)] <- 1
    }
    
  }
  
  #Update significance of genes by including up or down regulation
  
  comGeneDf$Gene <- rownames(comGeneDf)
  for (colName in comColNames) {
    
    if(colName == "Gene") next
    
    splitName <- strsplit(colName,"_")[[1]]
    if(length(splitName)==3){
      currentTissue <- splitName[1]
      currentContrast <- paste0(splitName[2],"_",splitName[3])
    }else if(length(splitName)==2){
      currentTissue <- splitName[1]
      currentContrast <- splitName[2]
    }
    
    limmaResults <- read.xlsx(paste0("data/LimmaResults/",currentTissue,"/",currentTissue,".xlsx"),sheet = currentContrast,rowNames = TRUE) 
    
    for (i in 1:nrow(comGeneDf)) {
      
      if(comGeneDf[i,colName]==1){
        logFC <- limmaResults$logFC[limmaResults$GeneSymbol==comGeneDf$Gene[i]]
        if(length(logFC)!=1){
          logFC <- limmaResults[limmaResults$GeneSymbol==comGeneDf$Gene[i],]
          logFC <- logFC[order(logFC$adj.P.Val),]
          logFC <- logFC$logFC[1]
          
        } 
        
        if(logFC < 0) comGeneDf[i,colName] <- -1
      }
      
    }
    
  }
  
  comGeneDf <- comGeneDf[,colnames(comGeneDf)!="Gene"]
  
  
  tissueOccurenceDf <- as.data.frame(matrix(ncol = length(tissueOfInterest),nrow = length(geneList),data = 0))
  rownames(tissueOccurenceDf) <- geneList
  colnames(tissueOccurenceDf) <- tissueOfInterest
  for (tissue in tissueOfInterest) {
    tissueColNames <- paste0(tissue,"_",contrasts)
    tissuecomGeneDf <- comGeneDf[,tissueColNames]
    tissuecomGeneDf <- tissuecomGeneDf %>%
      dplyr::mutate(TotalOccurence = rowSums(abs(.[,1:length(tissueColNames)])))
    tissueOccurenceDf[,tissue] <- tissuecomGeneDf$TotalOccurence
  }
  
  
  
  tissueOccurenceDf <- tissueOccurenceDf %>%
    mutate(TotalOccurence = rowSums(.[,1:length(tissueOfInterest)]))
  
  tissueOccurenceDf <- tissueOccurenceDf[order(tissueOccurenceDf$TotalOccurence,decreasing = TRUE),]
  
  
  write.csv(comGeneDf,"data/LimmaResults/ComprehensiveSignificantGeneInfo.csv",row.names = TRUE)
  write.csv(tissueOccurenceDf,"data/LimmaResults/TissueSignificantGeneInfo.csv",row.names = TRUE)
  
  
  comprehensiveData <- comGeneDf %>%
    mutate(TotalOccurence = rowSums(abs(.[,1:32])))
  
  comprehensiveData <- comprehensiveData[order(comprehensiveData$TotalOccurence,decreasing = TRUE),]
  
  comprehensiveData <- head(comprehensiveData,n=50)
  
  comprehensiveData <- comprehensiveData[,!colSums(abs(comprehensiveData))==0]
  
  comprehensiveData <- as.matrix(comprehensiveData)
  
  hmdata <- comprehensiveData[,-ncol(comprehensiveData)]
  
  
  
  
  
  hm <- Heatmap(hmdata,
                row_names_side = "left",
                column_names_side = "top",
                cluster_rows = FALSE,
                column_title = "Heatmap of significant genes(top 50) across all comparisons",
                column_title_gp = gpar(fontsize = 30, fontface = "bold"),
                row_names_gp = gpar(fontsize = 15),
                column_names_gp = gpar(fontsize = 15),
                column_title_side = "bottom",
                cluster_columns = FALSE,
                heatmap_width = unit(25,"cm"),
                heatmap_height = unit(40,"cm"),
                border = TRUE,
                show_heatmap_legend = FALSE,
                col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.rect(x, y, width, height, 
                            gp = gpar(col = "black", fill = NA))  # Draw cell borders
                }
  )
  
  lgd = Legend(
    labels = c("Up", "NotSig", "Down"),
    title = "Significance",
    legend_gp = gpar(fill = c("red", "white", "blue")),
    border = TRUE,
    direction = "horizontal", 
    title_gp = gpar(fontsize = 15),
    labels_gp = gpar(fontsize = 15)
  )
  
  
  jpeg(filename = "data/LimmaResults/ComprehensiveSignificantGeneHeatmap.jpeg",height = 1400,width = 1100, quality = 100)
  
  draw(hm,annotation_legend_list = list(lgd),
       annotation_legend_side = "left"
  )
  
  dev.off()
  
}


VisualizeMatches <- function(pathwayList, sheetNames, outFileName, topN = NULL) {
  
  # Prepare data frame for plotting
  matchData <- data.frame()
  
  for (j in seq_along(sheetNames)) {
    for (i in seq_along(pathwayList)) {
      # Skip self-comparison
      if (names(pathwayList[i]) == sheetNames[j]) next
      
      overlap <- intersect(pathwayList[[i]], pathwayList[[sheetNames[j]]])
      overlap <- na.omit(overlap)
      if (length(overlap) > 0) {
        matchData <- rbind(matchData,
                           data.frame(
                             Pathway = names(pathwayList[i]),
                             Reference = sheetNames[j],
                             Matches = length(overlap)
                           ))
      }
    }
  }
  
  # Optional: keep only top N pathways by total matches
  if(!is.null(topN)) {
    topPathways <- matchData %>%
      group_by(Pathway) %>%
      summarise(Total = sum(Matches)) %>%
      arrange(desc(Total)) %>%
      slice_head(n = topN)
    
    matchData <- matchData %>% filter(Pathway %in% topPathways$Pathway)
  }
  
  # Convert Matches to factor for discrete scale
  #matchData$MatchesFactor <- factor(matchData$Matches, levels = sort(unique(matchData$Matches)))
  
  # Plot heatmap with discrete scale and text labels
  p <- ggplot(matchData, aes(x = Reference, y = Pathway, fill = Matches)) +
    geom_tile(color = "black") +
    geom_text(aes(label = Matches), size = 5) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 16, hjust = 0.5)
    ) +
    labs(fill = "Number of Common genes", x = "Custom pathways", y = "Pathways in Database")+
    ggtitle("Common genes in custom and database pathways")
  
  # Save using ggsave
  dir.create(dirname(outFileName), recursive = TRUE, showWarnings = FALSE)
  ggsave(outFileName, plot = p, width = 12, height = 14, dpi = 300)
  
  return(matchData)
}
