

Compare_HFS_With_HFC <- function(cts,sampleInfo,geneInfo,tissueOfInterest){
  
  dirPath <- paste0("data/HFS_Vs_HFC")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  dirPath <- paste0("data/HFS_Vs_HFC/",tissueOfInterest)
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail == tissueOfInterest,]
  sampleInfo <- sampleInfo[sampleInfo$diet == "High fat",]
  
  cts <- cts[,colnames(cts) %in% sampleInfo$sampleid]
  
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
  
  
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  #Library size plot to understand the variability between samples
  libSizePlot <- ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=treatment)) +
    geom_bar(stat = "identity") +
    labs(y = "Library size (total number of mapped and quantified reads)",
         x = "Sample", fill = "Treatment") +
    coord_flip() 
  
  ggsave(paste0(dirPath,"/LibrarySizePlot.jpeg"),plot = libSizePlot,dpi = 300)
  
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
    ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes")) +
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
  densityPlot <- cowplot::plot_grid(beforeFiltering_plot, afterFiltering_plot)
  
  ggsave(paste0(dirPath,"/DensityPlot.jpeg"),plot = densityPlot,dpi = 600)
  
  
  #Filter the genes from the DGE object
  dge <- dge[flag,,keep.lib.sizes = FALSE]
  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + treatment, data = dge$samples)
  
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count=1)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    s_vs_c = (treatmentSurgery - treatmentControl)
  )
  
  contrastUsed <- "(treatmentSurgery - treatmentControl)"
  
  
  #Save design and contrast
  write.csv(as.data.frame(design),file = paste0("data/HFS_Vs_HFC/",tissueOfInterest,"/DesignMatrix.csv"))
  writeLines(contrastUsed, paste0("data/HFS_Vs_HFC/",tissueOfInterest,"/ContrastUsed.txt"))
  
  ###############################################Less weightage to LFC in gene ranking##############################
  
  fit <- lmFit(logCPM, design)%>%
    contrasts.fit(contrasts)
  
  
  fit <- eBayes(fit, trend=TRUE)
  
  
  #Extract all the differentially expressed genes
  allDEresults <- topTable(fit, 
                           coef = "s_vs_c", 
                           number = Inf, 
                           adjust.method = "fdr") %>%
    as.data.frame() 
  
  
  allDEresults$EnsemblID <-  rownames(allDEresults)
  ensemblIds <- allDEresults$EnsemblID
  
  geneNames <- geneInfo$gene_name[match(ensemblIds,geneInfo$ensembl_gene_id)]
  geneNames[is.na(geneNames)] <- ensemblIds[is.na(geneNames)]
  allDEresults$GeneSymbol <- geneNames
  
  
  rowNamesallDEresults <- paste0(allDEresults$EnsemblID,"_",allDEresults$GeneSymbol)
  rownames(allDEresults) <- rowNamesallDEresults
  
  allDEresults <- allDEresults[,c("EnsemblID","GeneSymbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
  
  
  
  write.csv(allDEresults,file = paste0("data/HFS_Vs_HFC/",tissueOfInterest,"/AllGenes.csv"))
  
  
  #Filter the differentially expressed genes based on FDR value < 0.05 and LFC absolute value > 1
  allDEresults <- allDEresults %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE, 
      TRUE ~ FALSE # If conditions in the line above are not met, gene is not DE. 
    ))
  
  
  
  sigDEresults <- allDEresults %>%
    dplyr::filter(isSignificant == TRUE)
  
  
  
  if(nrow(sigDEresults)!=0) write.csv(sigDEresults,file = paste0("data/HFS_Vs_HFC/",tissueOfInterest,"/SignificantGenes.csv"))
  
  sigdata <- filter(allDEresults,isSignificant == TRUE)
  
  top10 <- head(sigdata,n=10)
  
  
  #Volcano plot to visualize the results
  volcanoPlot <- allDEresults %>%
    ggplot(aes(x = logFC, 
               y = -log10(adj.P.Val),
               colour = isSignificant)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    geom_vline(xintercept = -1, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: HF_S Vs HF_C : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/VolcanoPlot.jpeg"),plot = volcanoPlot,height=12,width=10)
  
  
  #MA plot to visualize the results
  mAPlot <- allDEresults %>%
    ggplot(aes(x = AveExpr, 
               y = logFC,
               colour = isSignificant)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: HF_S Vs HF_C : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/MAPlot.jpeg"),plot = mAPlot,height=12,width=10)
  
  
  #Heatmap of top 10 genes
  if(nrow(sigDEresults)!=0){
    
    sigDEresults <- head(sigDEresults[order(sigDEresults$adj.P.Val),],n=50)
    
    sigGenes <- sigDEresults$EnsemblID
    
    sigGeneExpressionData <- cts[rownames(cts) %in% sigGenes,]
    
    rownames(sigGeneExpressionData) <- sigDEresults$GeneSymbol
    
    scaledSigExpData <- t(scale(t(sigGeneExpressionData)))
    
    jpeg(paste0(dirPath,"/TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          # row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          clustering_distance_columns = "euclidean"
    )
    
    print(topHeatmap)
    
    dev.off()
    
  }
  
  
}

Compare_LFS_With_LFC <- function(cts,sampleInfo,geneInfo,tissueOfInterest){
  
  dirPath <- paste0("data/LFS_Vs_LFC")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  dirPath <- paste0("data/LFS_Vs_LFC/",tissueOfInterest)
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail == tissueOfInterest,]
  sampleInfo <- sampleInfo[sampleInfo$diet == "Low fat",]
  
  cts <- cts[,colnames(cts) %in% sampleInfo$sampleid]
  
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
  
  
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  #Library size plot to understand the variability between samples
  libSizePlot <- ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=treatment)) +
    geom_bar(stat = "identity") +
    labs(y = "Library size (total number of mapped and quantified reads)",
         x = "Sample", fill = "Treatment") +
    coord_flip()
  
  ggsave(paste0(dirPath,"/LibrarySizePlot.jpeg"),plot = libSizePlot,dpi = 300)
  
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
    ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes")) +
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
  densityPlot <- cowplot::plot_grid(beforeFiltering_plot, afterFiltering_plot)
  
  ggsave(paste0(dirPath,"/DensityPlot.jpeg"),plot = densityPlot,dpi = 600)
  
  
  #Filter the genes from the DGE object
  dge <- dge[flag,,keep.lib.sizes = FALSE]
  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + treatment, data = dge$samples)
  
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count=1)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    s_vs_c = (treatmentSurgery - treatmentControl)
  )
  
  contrastUsed <- "(treatmentSurgery - treatmentControl)"
  
  
  #Save design and contrast
  write.csv(as.data.frame(design),file = paste0("data/LFS_Vs_LFC/",tissueOfInterest,"/DesignMatrix.csv"))
  writeLines(contrastUsed, paste0("data/LFS_Vs_LFC/",tissueOfInterest,"/ContrastUsed.txt"))
  
  
  ###############################################Less weightage to LFC in gene ranking##############################
  
  fit <- lmFit(logCPM, design)%>%
    contrasts.fit(contrasts)
  
  
  fit <- eBayes(fit, trend=TRUE)
  
  
  #Extract all the differentially expressed genes
  allDEresults <- topTable(fit, 
                           coef = "s_vs_c", 
                           number = Inf, 
                           adjust.method = "fdr") %>%
    as.data.frame()
  
  
  
  allDEresults$EnsemblID <-  rownames(allDEresults)
  ensemblIds <- allDEresults$EnsemblID
  
  geneNames <- geneInfo$gene_name[match(ensemblIds,geneInfo$ensembl_gene_id)]
  geneNames[is.na(geneNames)] <- ensemblIds[is.na(geneNames)]
  allDEresults$GeneSymbol <- geneNames
  
  
  rowNamesallDEresults <- paste0(allDEresults$EnsemblID,"_",allDEresults$GeneSymbol)
  rownames(allDEresults) <- rowNamesallDEresults
  
  allDEresults <- allDEresults[,c("EnsemblID","GeneSymbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
  
  
  
  
  write.csv(allDEresults,file = paste0("data/LFS_Vs_LFC/",tissueOfInterest,"/AllGenes.csv"))
  
  
  #Filter the differentially expressed genes based on FDR value < 0.05 and LFC absolute value > 1
  allDEresults <- allDEresults %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE, 
      TRUE ~ FALSE # If conditions in the line above are not met, gene is not DE. 
    ))
  
  
  
  sigDEresults <- allDEresults %>%
    dplyr::filter(isSignificant == TRUE)
  
  
  
  if(nrow(sigDEresults)!=0) write.csv(sigDEresults,file = paste0("data/LFS_Vs_LFC/",tissueOfInterest,"/SignificantGenes.csv"))
  
  sigdata <- filter(allDEresults,isSignificant == TRUE)
  
  top10 <- head(sigdata,n=10)
  
  
  #Volcano plot to visualize the results
  volcanoPlot <- allDEresults %>%
    ggplot(aes(x = logFC, 
               y = -log10(adj.P.Val),
               colour = isSignificant)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    geom_vline(xintercept = -1, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: LF_S Vs LF_C : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/VolcanoPlot.jpeg"),plot = volcanoPlot,height=12,width=10)
  
  
  #MA plot to visualize the results
  mAPlot <- allDEresults %>%
    ggplot(aes(x = AveExpr, 
               y = logFC,
               colour = isSignificant)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: LF_S Vs LF_C : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/MAPlot.jpeg"),plot = mAPlot,height=12,width=10)
  
  
  #Heatmap of top 10 genes
  if(nrow(sigDEresults)!=0){
    
    sigDEresults <- head(sigDEresults[order(sigDEresults$adj.P.Val),],n=50)
    
    sigGenes <- sigDEresults$EnsemblID
    
    sigGeneExpressionData <- cts[rownames(cts) %in% sigGenes,]
    
    rownames(sigGeneExpressionData) <- sigDEresults$GeneSymbol
    
    scaledSigExpData <- t(scale(t(sigGeneExpressionData)))
    
    jpeg(paste0(dirPath,"/TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          # row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          clustering_distance_columns = "euclidean"
    )
    
    print(topHeatmap)
    
    dev.off()
    
  }
  
}

Compare_HF_With_LF <- function(cts,sampleInfo,geneInfo,tissueOfInterest){
  
  dirPath <- paste0("data/HF_Vs_LF")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  dirPath <- paste0("data/HF_Vs_LF/",tissueOfInterest)
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail == tissueOfInterest,]
  
  cts <- cts[,colnames(cts) %in% sampleInfo$sampleid]
  
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
  
  
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  #Library size plot to understand the variability between samples
  libSizePlot <- ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=diet)) +
    geom_bar(stat = "identity") +
    labs(y = "Library size (total number of mapped and quantified reads)",
         x = "Sample", fill = "Diet") +
    coord_flip()
  
  ggsave(paste0(dirPath,"/LibrarySizePlot.jpeg"),plot = libSizePlot,dpi = 300)
  
  #Variability between samples
  max(dge$samples$lib.size)/min(dge$samples$lib.size)
  
  #Set a threshold for filtering genes
  flag <- (rowSums(cpm(dge, log = TRUE, prior.count = 1) > 1) >= smallestSampleNumber)
  
  
  #Generate a density plot before filtering
  beforeFiltering_plot <- dge %>% 
    cpm(log = TRUE, prior.count = 1) %>% 
    melt %>% 
    dplyr::filter(is.finite(value)) %>% 
    ggplot(aes(x = value, colour = Var2)) +
    geom_density() + 
    guides(colour = FALSE,scale = 'none') +
    ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes")) +
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
  densityPlot <- cowplot::plot_grid(beforeFiltering_plot, afterFiltering_plot)
  
  ggsave(paste0(dirPath,"/DensityPlot.jpeg"),plot = densityPlot,dpi = 600)
  
  
  #Filter the genes from the DGE object
  dge <- dge[flag,,keep.lib.sizes = FALSE]
  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + diet, data = dge$samples)
  colnames(design) <- c("dietHigh_fat", "dietLow_fat")
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count = 1)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    h_vs_l = (dietHigh_fat - dietLow_fat)
  )
  
  contrastUsed <- "(dietHigh_fat - dietLow_fat)"
  
  
  #Save design and contrast
  write.csv(as.data.frame(design),file = paste0("data/HF_Vs_LF/",tissueOfInterest,"/DesignMatrix.csv"))
  writeLines(contrastUsed, paste0("data/HF_Vs_LF/",tissueOfInterest,"/ContrastUsed.txt"))
  
  ###############################################Less weightage to LFC in gene ranking##############################
  
  fit <- lmFit(logCPM, design)%>%
    contrasts.fit(contrasts)
  
  
  fit <- eBayes(fit, trend=TRUE)
  
  
  #Extract all the differentially expressed genes
  allDEresults <- topTable(fit, 
                           coef = "h_vs_l", 
                           number = Inf, 
                           adjust.method = "fdr") %>%
    as.data.frame()
  
  
  allDEresults$EnsemblID <-  rownames(allDEresults)
  ensemblIds <- allDEresults$EnsemblID
  
  geneNames <- geneInfo$gene_name[match(ensemblIds,geneInfo$ensembl_gene_id)]
  geneNames[is.na(geneNames)] <- ensemblIds[is.na(geneNames)]
  allDEresults$GeneSymbol <- geneNames
  
  
  rowNamesallDEresults <- paste0(allDEresults$EnsemblID,"_",allDEresults$GeneSymbol)
  rownames(allDEresults) <- rowNamesallDEresults
  
  allDEresults <- allDEresults[,c("EnsemblID","GeneSymbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
  
  
  
  write.csv(allDEresults,file = paste0("data/HF_Vs_LF/",tissueOfInterest,"/AllGenes.csv"))
  
  
  #Filter the differentially expressed genes based on FDR value < 0.05 and LFC absolute value > 1
  allDEresults <- allDEresults %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE, 
      TRUE ~ FALSE # If conditions in the line above are not met, gene is not DE. 
    ))
  
  
  
  sigDEresults <- allDEresults %>%
    dplyr::filter(isSignificant == TRUE)
  
  
  
  if(nrow(sigDEresults)!=0) write.csv(sigDEresults,file = paste0("data/HF_Vs_LF/",tissueOfInterest,"/SignificantGenes.csv"))
  
  sigdata <- filter(allDEresults,isSignificant == TRUE)
  
  top10 <- head(sigdata,n=10)
  
  
  #Volcano plot to visualize the results
  volcanoPlot <- allDEresults %>%
    ggplot(aes(x = logFC, 
               y = -log10(adj.P.Val),
               colour = isSignificant)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    geom_vline(xintercept = -1, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: HF Vs LF : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/VolcanoPlot.jpeg"),plot = volcanoPlot,height=12,width=10)
  
  
  #MA plot to visualize the results
  mAPlot <- allDEresults %>%
    ggplot(aes(x = AveExpr, 
               y = logFC,
               colour = isSignificant)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: HF Vs LF : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/MAPlot.jpeg"),plot = mAPlot,height=12,width=10)
  
  
  #Heatmap of top 10 genes
  if(nrow(sigDEresults)!=0){
    
    sigDEresults <- head(sigDEresults[order(sigDEresults$adj.P.Val),],n=50)
    
    sigGenes <- sigDEresults$EnsemblID
    
    sigGeneExpressionData <- cts[rownames(cts) %in% sigGenes,]
    
    rownames(sigGeneExpressionData) <- sigDEresults$GeneSymbol
    
    scaledSigExpData <- t(scale(t(sigGeneExpressionData)))
    
    jpeg(paste0(dirPath,"/TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          # row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          clustering_distance_columns = "euclidean"
    )
    
    print(topHeatmap)
    
    dev.off()
    
  }
  
}

Compare_Surgery_With_Diet <- function(cts,sampleInfo,geneInfo,tissueOfInterest){
  
  dirPath <- paste0("data/Surgery_With_Diet")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  dirPath <- paste0("data/Surgery_With_Diet/",tissueOfInterest)
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail == tissueOfInterest,]
  
  cts <- cts[,colnames(cts) %in% sampleInfo$sampleid]
  
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
  
  
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  #Library size plot to understand the variability between samples
  libSizePlot <- ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=diet)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ treatment) +
    labs(y = "Library size (total number of mapped and quantified reads)",
         x = "Sample", fill = "Diet") +
    coord_flip()
  
  ggsave(paste0(dirPath,"/LibrarySizePlot.jpeg"),plot = libSizePlot,dpi = 300)
  
  #Variability between samples
  max(dge$samples$lib.size)/min(dge$samples$lib.size)
  
  #Set a threshold for filtering genes
  flag <- (rowSums(cpm(dge, log = TRUE, prior.count = 1) > 1) >= smallestSampleNumber)
  
  
  #Generate a density plot before filtering
  beforeFiltering_plot <- dge %>% 
    cpm(log = TRUE, prior.count = 1) %>% 
    melt %>% 
    dplyr::filter(is.finite(value)) %>% 
    ggplot(aes(x = value, colour = Var2)) +
    geom_density() + 
    guides(colour = FALSE,scale = 'none') +
    ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes")) +
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
  densityPlot <- cowplot::plot_grid(beforeFiltering_plot, afterFiltering_plot)
  
  ggsave(paste0(dirPath,"/DensityPlot.jpeg"),plot = densityPlot,dpi = 600)
  
  
  #Filter the genes from the DGE object
  dge <- dge[flag,,keep.lib.sizes = FALSE]
  
  #Create an interaction term in the sampleinfo
  dge$samples$diet <- factor(dge$samples$diet)
  dge$samples$treatment <- factor(dge$samples$treatment)
  dge$samples$diet_treatment <- interaction(dge$samples$diet,dge$samples$treatment)
  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + diet_treatment, data = dge$samples)
  colnames(design) <- c("Highfat_Control", "Lowfat_Control","Highfat_Surgery","Lowfat_Surgery")
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count=1)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    hfsvshf_vs_lfsvslf = ((Highfat_Surgery - Highfat_Control)/2) - ((Lowfat_Surgery - Lowfat_Control)/2)
  )
  
  contrastUsed <- "((Highfat_Surgery - Highfat_Control)/2) - ((Lowfat_Surgery - Lowfat_Control)/2)"
  
  
  #Save design and contrast
  write.csv(as.data.frame(design),file = paste0("data/Surgery_With_Diet/",tissueOfInterest,"/DesignMatrix.csv"))
  writeLines(contrastUsed, paste0("data/Surgery_With_Diet/",tissueOfInterest,"/ContrastUsed.txt"))
  
  ###############################################Less weightage to LFC in gene ranking##############################
  
  fit <- lmFit(logCPM, design)%>%
    contrasts.fit(contrasts)
  
  
  fit <- eBayes(fit, trend=TRUE)
  
  
  #Extract all the differentially expressed genes
  allDEresults <- topTable(fit, 
                           coef = "hfsvshf_vs_lfsvslf", 
                           number = Inf, 
                           adjust.method = "fdr") %>%
    as.data.frame()
  
  
  allDEresults$EnsemblID <-  rownames(allDEresults)
  ensemblIds <- allDEresults$EnsemblID
  
  geneNames <- geneInfo$gene_name[match(ensemblIds,geneInfo$ensembl_gene_id)]
  geneNames[is.na(geneNames)] <- ensemblIds[is.na(geneNames)]
  allDEresults$GeneSymbol <- geneNames
  
  
  rowNamesallDEresults <- paste0(allDEresults$EnsemblID,"_",allDEresults$GeneSymbol)
  rownames(allDEresults) <- rowNamesallDEresults
  
  allDEresults <- allDEresults[,c("EnsemblID","GeneSymbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
  
  
  write.csv(allDEresults,file = paste0("data/Surgery_With_Diet/",tissueOfInterest,"/AllGenes.csv"))
  
  
  #Filter the differentially expressed genes based on FDR value < 0.05 and LFC absolute value > 1
  allDEresults <- allDEresults %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE, 
      TRUE ~ FALSE # If conditions in the line above are not met, gene is not DE. 
    ))
  
  
  
  sigDEresults <- allDEresults %>%
    dplyr::filter(isSignificant == TRUE)
  
  
  
  if(nrow(sigDEresults)!=0) write.csv(sigDEresults,file = paste0("data/Surgery_With_Diet/",tissueOfInterest,"/SignificantGenes.csv"))
  
  sigdata <- filter(allDEresults,isSignificant == TRUE)
  
  top10 <- head(sigdata,n=10)
  
  
  #Volcano plot to visualize the results
  volcanoPlot <- allDEresults %>%
    ggplot(aes(x = logFC, 
               y = -log10(adj.P.Val),
               colour = isSignificant)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    geom_vline(xintercept = -1, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: Surgery With Diet effect : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/VolcanoPlot.jpeg"),plot = volcanoPlot,height=12,width=10)
  
  
  #MA plot to visualize the results
  mAPlot <- allDEresults %>%
    ggplot(aes(x = AveExpr, 
               y = logFC,
               colour = isSignificant)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: Surgery With Diet effect : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/MAPlot.jpeg"),plot = mAPlot,height=12,width=10)
  
  
  #Heatmap of top 10 genes
  if(nrow(sigDEresults)!=0){
    
    sigDEresults <- head(sigDEresults[order(sigDEresults$adj.P.Val),],n=50)
    
    sigGenes <- sigDEresults$EnsemblID
    
    sigGeneExpressionData <- cts[rownames(cts) %in% sigGenes,]
    
    rownames(sigGeneExpressionData) <- sigDEresults$GeneSymbol
    
    scaledSigExpData <- t(scale(t(sigGeneExpressionData)))
    
    jpeg(paste0(dirPath,"/TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          # row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          clustering_distance_columns = "euclidean"
    )
    
    print(topHeatmap)
    
    dev.off()
    
  }
  
}

Compare_Surgery_Without_Diet <- function(cts,sampleInfo,geneInfo,tissueOfInterest){
  
  dirPath <- paste0("data/Surgery_Without_Diet")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  dirPath <- paste0("data/Surgery_Without_Diet/",tissueOfInterest)
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail == tissueOfInterest,]
  
  cts <- cts[,colnames(cts) %in% sampleInfo$sampleid]
  
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
  
  
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  #Library size plot to understand the variability between samples
  libSizePlot <- ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=diet)) +
    geom_bar(stat = "identity") +
    facet_wrap(~ treatment) +
    labs(y = "Library size (total number of mapped and quantified reads)",
         x = "Sample", fill = "Diet") +
    coord_flip()
  
  ggsave(paste0(dirPath,"/LibrarySizePlot.jpeg"),plot = libSizePlot,dpi = 300)
  
  #Variability between samples
  max(dge$samples$lib.size)/min(dge$samples$lib.size)
  
  #Set a threshold for filtering genes
  flag <- (rowSums(cpm(dge, log = TRUE, prior.count = 1) > 1) >= smallestSampleNumber)
  
  
  #Generate a density plot before filtering
  beforeFiltering_plot <- dge %>% 
    cpm(log = TRUE, prior.count = 1) %>% 
    melt %>% 
    dplyr::filter(is.finite(value)) %>% 
    ggplot(aes(x = value, colour = Var2)) +
    geom_density() + 
    guides(colour = FALSE,scale = 'none') +
    ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes")) +
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
  densityPlot <- cowplot::plot_grid(beforeFiltering_plot, afterFiltering_plot)
  
  ggsave(paste0(dirPath,"/DensityPlot.jpeg"),plot = densityPlot,dpi = 600)
  
  
  #Filter the genes from the DGE object
  dge <- dge[flag,,keep.lib.sizes = FALSE]
  
  #Create an interaction term in the sampleinfo
  dge$samples$diet <- factor(dge$samples$diet)
  dge$samples$treatment <- factor(dge$samples$treatment)
  dge$samples$diet_treatment <- interaction(dge$samples$diet,dge$samples$treatment)
  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + diet_treatment, data = dge$samples)
  colnames(design) <- c("Highfat_Control", "Lowfat_Control","Highfat_Surgery","Lowfat_Surgery")
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count=1)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    hfslfs_vs_hfclfc = ((Highfat_Surgery + Lowfat_Surgery)/2) - ((Highfat_Control + Lowfat_Control)/2)
  )
  
  contrastUsed <- "((Highfat_Surgery + Lowfat_Surgery)/2) - ((Highfat_Control + Lowfat_Control)/2)"
  
  
  #Save design and contrast
  write.csv(as.data.frame(design),file = paste0("data/Surgery_Without_Diet/",tissueOfInterest,"/DesignMatrix.csv"))
  writeLines(contrastUsed, paste0("data/Surgery_Without_Diet/",tissueOfInterest,"/ContrastUsed.txt"))
  
  
  ###############################################Less weightage to LFC in gene ranking##############################
  
  fit <- lmFit(logCPM, design)%>%
    contrasts.fit(contrasts)
  
  
  fit <- eBayes(fit, trend=TRUE)
  
  
  #Extract all the differentially expressed genes
  allDEresults <- topTable(fit, 
                           coef = "hfslfs_vs_hfclfc", 
                           number = Inf, 
                           adjust.method = "fdr") %>%
    as.data.frame()

  allDEresults$EnsemblID <-  rownames(allDEresults)
  ensemblIds <- allDEresults$EnsemblID
  
  geneNames <- geneInfo$gene_name[match(ensemblIds,geneInfo$ensembl_gene_id)]
  geneNames[is.na(geneNames)] <- ensemblIds[is.na(geneNames)]
  allDEresults$GeneSymbol <- geneNames
  
  
  rowNamesallDEresults <- paste0(allDEresults$EnsemblID,"_",allDEresults$GeneSymbol)
  rownames(allDEresults) <- rowNamesallDEresults
  
  allDEresults <- allDEresults[,c("EnsemblID","GeneSymbol","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
  
  
  write.csv(allDEresults,file = paste0("data/Surgery_Without_Diet/",tissueOfInterest,"/AllGenes.csv"))
  
  
  #Filter the differentially expressed genes based on FDR value < 0.05 and LFC absolute value > 1
  allDEresults <- allDEresults %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE, 
      TRUE ~ FALSE # If conditions in the line above are not met, gene is not DE. 
    ))
  
  
  
  sigDEresults <- allDEresults %>%
    dplyr::filter(isSignificant == TRUE)
  
  
  
  if(nrow(sigDEresults)!= 0) write.csv(sigDEresults,file = paste0("data/Surgery_Without_Diet/",tissueOfInterest,"/SignificantGenes.csv"))
  
  
  sigdata <- filter(allDEresults,isSignificant == TRUE)
  
  top10 <- head(sigdata,n=10)
  
  
  #Volcano plot to visualize the results
  volcanoPlot <- allDEresults %>%
    ggplot(aes(x = logFC, 
               y = -log10(adj.P.Val),
               colour = isSignificant)) +
    geom_vline(xintercept = 1, linetype = "dotted") +
    geom_vline(xintercept = -1, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: Surgery Without Diet effect : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/VolcanoPlot.jpeg"),plot = volcanoPlot,height=12,width=10)
  
  
  #MA plot to visualize the results
  mAPlot <- allDEresults %>%
    ggplot(aes(x = AveExpr, 
               y = logFC,
               colour = isSignificant)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
    ggtitle(paste0("Differential gene expression results: Surgery Without Diet effect : ",tissueOfInterest," tissue")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  ggsave(paste0(dirPath,"/MAPlot.jpeg"),plot = mAPlot,height=12,width=10)
  
  
  #Heatmap of top 10 genes
  if(nrow(sigDEresults)!=0){
    
    sigDEresults <- head(sigDEresults[order(sigDEresults$adj.P.Val),],n=50)
    
    sigGenes <- sigDEresults$EnsemblID
    
    sigGeneExpressionData <- cts[rownames(cts) %in% sigGenes,]
    
    rownames(sigGeneExpressionData) <- sigDEresults$GeneSymbol
    
    scaledSigExpData <- t(scale(t(sigGeneExpressionData)))
    
    jpeg(paste0(dirPath,"/TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          # row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          clustering_distance_columns = "euclidean"
    )
    
    print(topHeatmap)
    
    dev.off()
    
  }
  
}

Compare_MultipleContrasts <- function(cts,sampleInfo,geneInfo,tissueOfInterest){
  
  dirPath <- paste0("data/LimmaResults")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  dirPath <- paste0("data/LimmaResults/",tissueOfInterest,"/")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  if(tissueOfInterest=="Small Bowel Mucosa"){
    
    sampleInfo <- sampleInfo[sampleInfo$tissuegeneral == tissueOfInterest,]
    sampleInfo <- sampleInfo[!sampleInfo$tissuedetail == "Ileum interposition",]
    
  }else{
    
    sampleInfo <- sampleInfo[sampleInfo$tissuedetail == tissueOfInterest,]
    
  }
  
  
  
  
  cts <- cts[,colnames(cts) %in% sampleInfo$sampleid]
  
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
  
  
  
  #Create a DGE object using the counts, sampleInfo and geneInfo
  dge <- DGEList(counts = cts,genes = geneInfo,samples = sampleInfo,remove.zeros = TRUE)
  
  #Normalization of the counts using TMM method
  dge <- calcNormFactors(dge,method = "TMM")
  
  #Library size plot to understand the variability between samples
  libSizePlot <- ggplot(data = dge$samples,aes(x = sampleid,y = lib.size,fill=treatment)) +
    facet_wrap(~diet)+
    geom_bar(stat = "identity") +
    labs(y = "Library size (total number of mapped and quantified reads)",
         x = "Sample", fill = "Treatment") +
    coord_flip() 
  
  ggsave(paste0(dirPath,"LibrarySizePlot.jpeg"),plot = libSizePlot,dpi = 300)
  
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
    ggtitle("A. Before filtering", subtitle = paste0(nrow(dge), " genes")) +
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
  densityPlot <- cowplot::plot_grid(beforeFiltering_plot, afterFiltering_plot)
  
  ggsave(paste0(dirPath,"DensityPlot.jpeg"),plot = densityPlot,dpi = 600)
  
  
  #Filter the genes from the DGE object
  dge <- dge[flag,,keep.lib.sizes = FALSE]
  
  #Create an interaction term in the sampleinfo
  dge$samples$diet <- factor(dge$samples$diet)
  dge$samples$treatment <- factor(dge$samples$treatment)
  dge$samples$diet_treatment <- interaction(dge$samples$diet,dge$samples$treatment)
  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + diet_treatment, data = dge$samples)
  colnames(design) <- c("Highfat_Control", "Lowfat_Control","Highfat_Surgery","Lowfat_Surgery")
  
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count=1)
  
  #Save logCPM counts
  SaveLogCPMCounts(logCPM,dirPath)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    HF_Surgery = (Highfat_Surgery - Highfat_Control),
    LF_Surgery = (Lowfat_Surgery - Lowfat_Control),
    Diet = (((Highfat_Surgery + Highfat_Control)/2) - ((Lowfat_Surgery + Lowfat_Control)/2)),
    Surgery_Diet = ((Highfat_Surgery - Highfat_Control)/2) - ((Lowfat_Surgery - Lowfat_Control)/2),
    Surgery = ((Highfat_Surgery + Lowfat_Surgery)/2) - ((Highfat_Control + Lowfat_Control)/2)
  )
  
  contrastsUsed <- c(
    "HF_Surgery = (Highfat_Surgery - Highfat_Control)",
                     "LF_Surgery = (Lowfat_Surgery - Lowfat_Control)",
                     "Diet = ((((Highfat_Surgery + Highfat_Control)/2)) - (((Lowfat_Surgery + Lowfat_Control)/2)) )",
                     "Surgery_Diet = ((Highfat_Surgery - Highfat_Control)/2) - ((Lowfat_Surgery - Lowfat_Control)/2)",
                     "Surgery = ((Highfat_Surgery + Lowfat_Surgery)/2) - ((Highfat_Control + Lowfat_Control)/2)")
  
  
  
  #Save design and contrast
  # write.csv(as.data.frame(design),file = paste0("data/LimmaResults/",tissueOfInterest,"_DesignMatrix.csv"))
  writeLines(contrastsUsed, paste0("data/LimmaResults/",tissueOfInterest,"/ContrastsUsed.txt"))
  
  ###############################################Less weightage to LFC in gene ranking##############################
  
  fit <- lmFit(logCPM, design)%>%
    contrasts.fit(contrasts)
  
  
  fit <- eBayes(fit, trend=TRUE)
  
  
  #Save all and significant genes
  for (contrast in colnames(contrasts)) {
    
    SaveAllNSigGenes(fit,contrast,geneInfo,dirPath)
    
  }
  
  
  #Plot true null fractions
  PlotTrueNullFractions(dirPath,colnames(contrasts))
  
  
  
  #Plot biological and technical correlations
  # PlotBioNTechCorr(fit,colnames(contrasts),dirPath)
  
  
  #Plot Venn Diagrams
  PlotVennDiagrams(fit,colnames(contrasts),geneInfo,tissueOfInterest,dirPath)
  
  
  #Plot Volcano plots
  PlotVolcano(colnames(contrasts),tissueOfInterest,dirPath)
  
  #Plot pValue histograms
  PlotPValueHistogram(colnames(contrasts),tissueOfInterest,dirPath)
  
  
  #Plot MA plots
  PlotMA(colnames(contrasts),tissueOfInterest,dirPath)

  #Plot heatmaps
  PlotHeatmaps(logCPM,colnames(contrasts),tissueOfInterest,dirPath)
  
  
}