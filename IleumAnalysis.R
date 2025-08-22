
#######################################################Ileum Interposition Vs Control#################################################################
#Read filtered counts and sample information
cts <- as.data.frame(read.csv("data/Filtered_counts.csv",row.names = 1))
sampleInfo <- as.data.frame(read.csv("data/Filtered_sampleInfo.csv",row.names = 1))
geneInfo <- as.data.frame(read.csv("data/geneInfoCorrected.csv",row.names = 1))


#Extracting Ileum Interposition, Ileum Control and Ileum Surgery samples
sampleInfo <- sampleInfo[sampleInfo$tissuedetail %in% c("Ileum","Ileum interposition"),]
sampleInfo <- sampleInfo[!(sampleInfo$tissuedetail == "Ileum" & sampleInfo$treatment == "Surgery"),]
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



dirPath <- paste0("data/PCAplots/IleumIVsC")
# Check if the directory exists to save plot
if (!dir.exists(dirPath)) {
  # Create the directory
  dir.create(dirPath)
}
  

plotWidth <- 12
plotHeight <- 10
titleSize <- 15
labelSize <- 5
tissueOfInterest <- "Ileum interposition Vs Control"
  
  
  pcaPlot <- ggbiplot::ggbiplot(pcaAnalysis,
                                groups = dge$samples$samplename,
                                # labels = dge$samples$sampleid,
                                var.axes = FALSE) +
    geom_text_repel(aes(label = dge$samples$sampleid), size = labelSize,max.overlaps = Inf) +
    ggtitle(paste0("PCA Plot of RNA-seq Data for ",tissueOfInterest)) +   
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
      title = paste0("Scree Plot of PCA for ",tissueOfInterest),
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
  
  dirPath <- paste0("data/LimmaResults")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  tissueOfInterest <- "Ileum_interposition_Control"
  dirPath <- paste0("data/LimmaResults/",tissueOfInterest,"/")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  #Read filtered counts and sample information
  cts <- as.data.frame(read.csv("data/Filtered_counts.csv",row.names = 1))
  sampleInfo <- as.data.frame(read.csv("data/Filtered_sampleInfo.csv",row.names = 1))
  geneInfo <- as.data.frame(read.csv("data/geneInfoCorrected.csv",row.names = 1))
  
  
  #Extracting Ileum Interposition, Ileum Control and Ileum Surgery samples
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail %in% c("Ileum","Ileum interposition"),]
  sampleInfo <- sampleInfo[!(sampleInfo$tissuedetail == "Ileum" & sampleInfo$treatment == "Surgery"),]
  sampleInfo <- sampleInfo[!(sampleInfo$sampleid == "LF_S_II_R1_B1"),]
  sampleInfo <- sampleInfo[order(sampleInfo$tissuedetail),]
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
    facet_wrap(~tissuedetail)+
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
  
  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + tissuedetail, data = dge$samples)
  colnames(design) <- c("Ileum_Control","Ileum_Interposition")
  
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count=1)
  
  #Save logCPM counts
  SaveLogCPMCounts(logCPM,dirPath)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    Ileum_Interposition_Control = Ileum_Interposition -Ileum_Control
  )
  
  contrastsUsed <- c(
    "Ileum_Interposition_Control = Ileum_Interposition -Ileum_Control")
  
  
  
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
  #PlotTrueNullFractions(dirPath,colnames(contrasts))
  
  
  
  #Plot biological and technical correlations
  # PlotBioNTechCorr(fit,colnames(contrasts),dirPath)
  
  
  #Plot Venn Diagrams
  #PlotVennDiagramsIleum(fit,colnames(contrasts),geneInfo,tissueOfInterest,dirPath)
  
  
  #Plot Volcano plots
  PlotVolcanoIleum(colnames(contrasts),tissueOfInterest,dirPath)
  
  #Plot pValue histograms
  PlotPValueHistogramIleum(colnames(contrasts),tissueOfInterest,dirPath)
  
  
  #Plot MA plots
  PlotMAIleum(colnames(contrasts),tissueOfInterest,dirPath)
  
  #Plot heatmaps #Check sample names in function before running
  PlotHeatmapsIleum(logCPM,colnames(contrasts),tissueOfInterest,dirPath)
  
  
  
  #######################################################Ileum Interposition Vs Surgery#################################################################
  #Read filtered counts and sample information
  cts <- as.data.frame(read.csv("data/Filtered_counts.csv",row.names = 1))
  sampleInfo <- as.data.frame(read.csv("data/Filtered_sampleInfo.csv",row.names = 1))
  geneInfo <- as.data.frame(read.csv("data/geneInfoCorrected.csv",row.names = 1))
  
  
  #Extracting Ileum Interposition, Ileum Control and Ileum Surgery samples
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail %in% c("Ileum","Ileum interposition"),]
  sampleInfo <- sampleInfo[!(sampleInfo$tissuedetail == "Ileum" & sampleInfo$treatment == "Control"),]
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
  
  
  
  dirPath <- paste0("data/PCAplots/IleumIVsS")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  plotWidth <- 12
  plotHeight <- 10
  titleSize <- 15
  labelSize <- 5
  tissueOfInterest <- "Ileum interposition Vs Surgery"
  
  
  pcaPlot <- ggbiplot::ggbiplot(pcaAnalysis,
                                groups = dge$samples$samplename,
                                # labels = dge$samples$sampleid,
                                var.axes = FALSE) +
    geom_text_repel(aes(label = dge$samples$sampleid), size = labelSize,max.overlaps = Inf) +
    ggtitle(paste0("PCA Plot of RNA-seq Data for ",tissueOfInterest)) +   
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
      title = paste0("Scree Plot of PCA for ",tissueOfInterest),
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
  
  
  
  dirPath <- paste0("data/LimmaResults")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  tissueOfInterest <- "Ileum_interposition_Surgery"
  dirPath <- paste0("data/LimmaResults/",tissueOfInterest,"/")
  # Check if the directory exists to save plot
  if (!dir.exists(dirPath)) {
    # Create the directory
    dir.create(dirPath)
  }
  
  
  #Read filtered counts and sample information
  cts <- as.data.frame(read.csv("data/Filtered_counts.csv",row.names = 1))
  sampleInfo <- as.data.frame(read.csv("data/Filtered_sampleInfo.csv",row.names = 1))
  geneInfo <- as.data.frame(read.csv("data/geneInfoCorrected.csv",row.names = 1))
  
  
  #Extracting Ileum Interposition, Ileum Control and Ileum Surgery samples
  sampleInfo <- sampleInfo[sampleInfo$tissuedetail %in% c("Ileum","Ileum interposition"),]
  sampleInfo <- sampleInfo[!(sampleInfo$tissuedetail == "Ileum" & sampleInfo$treatment == "Control"),]
  sampleInfo <- sampleInfo[order(sampleInfo$tissuedetail),]
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
    facet_wrap(~tissuedetail)+
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

  
  #Design a model for comparison of samples
  design <- model.matrix(~ 0 + tissuedetail, data = dge$samples)
  colnames(design) <- c("Ileum_Surgery","Ileum_Interposition")
  
  
  ############################LIMMA-TREND APPROACH###########################################################
  
  
  logCPM <- cpm(dge, log=TRUE, prior.count=1)
  
  #Save logCPM counts
  SaveLogCPMCounts(logCPM,dirPath)
  
  #Design a contrast for pairwise comparison of samples
  contrasts <- makeContrasts(
    levels = colnames(design),
    Ileum_Interposition_Surgery = Ileum_Interposition -Ileum_Surgery
  )
  
  contrastsUsed <- c(
    "Ileum_Interposition_Surgery = Ileum_Interposition -Ileum_Surgery")
  
  
  
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
  #PlotTrueNullFractions(dirPath,colnames(contrasts))
  
  
  
  #Plot biological and technical correlations
  # PlotBioNTechCorr(fit,colnames(contrasts),dirPath)
  
  
  #Plot Venn Diagrams
  #PlotVennDiagramsIleum(fit,colnames(contrasts),geneInfo,tissueOfInterest,dirPath)
  
  
  #Plot Volcano plots
  PlotVolcanoIleum(colnames(contrasts),tissueOfInterest,dirPath)
  
  #Plot pValue histograms
  PlotPValueHistogramIleum(colnames(contrasts),tissueOfInterest,dirPath)
  
  
  #Plot MA plots
  PlotMAIleum(colnames(contrasts),tissueOfInterest,dirPath)
  
  #Plot heatmaps
  PlotHeatmapsIleum(logCPM,colnames(contrasts),tissueOfInterest,dirPath)
  
  
  
  #######################################Pathway analysis for all contrasts#####################################################
  hallMarkPathways <- CreateMouseGenomeHallMarkEnsemblPathways()
  wikiPathways <- CreateMouseGenomeWikiEnsemblPathways()
  
  wb <- loadWorkbook("data/Bile\ Acid\ Genes_Sujoy.xlsx")
  sheetNames <- sheets(wb)
  
  for(i in 1:length(sheetNames)){
    
    v1 <- read.xlsx(wb, sheet = sheetNames[i])
    v1 <- v1$Gene
    v1 <- ConvertSymbolsToEnsembl(toupper(v1))
    wikiPathways[sheetNames[i]] <- list(v1)
    hallMarkPathways[sheetNames[i]] <- list(v1)
    
  }
  
  pathwaysList <- list(hallMarkPathways,wikiPathways)
  pathwaysListNames <- c("HallMark","Wiki")
  
  DoFGSEAForAllGenesIleum("Ileum_interposition_Surgery",pathways = pathwaysList,pathwayNames = pathwaysListNames)
  DoFGSEAForAllGenesIleum("Ileum_interposition_Control",pathways = pathwaysList,pathwayNames = pathwaysListNames)
  
  SaveAndPlotTopPathwaysIleum("Ileum_interposition_Surgery",pathways = pathwaysList,pathwayNames = pathwaysListNames)
  SaveAndPlotTopPathwaysIleum("Ileum_interposition_Control",pathways = pathwaysList,pathwayNames = pathwaysListNames)
  
  
  GenerateBarPlotsForPathwaysIleum("Ileum_interposition_Surgery", pathwaysListNames = pathwaysListNames, contrasts = "Ileum_interposition_Surgery")
  GenerateBarPlotsForPathwaysIleum("Ileum_interposition_Control", pathwaysListNames = pathwaysListNames, contrasts = "Ileum_interposition_Control")
  
  
  outFile <- "data/HallmarkMatches.txt"
  messages <- c()
  
  for (j in 1:length(sheetNames)) {
    for (i in 1:50) {
      overlap <- intersect(hallMarkPathways[[i]], hallMarkPathways[[sheetNames[j]]])
      if (length(overlap) > 0) {
        msg <- paste0(
          names(hallMarkPathways[i]), 
          " has ", 
          length(overlap),
          " matches with ",
          sheetNames[j]
        )
        messages <- c(messages, msg)
      }
    }
  }
  
  writeLines(messages, outFile)
  
  
  outFile <- "data/WikiMatches.txt"
  messages <- c()
  
  for (j in 1:length(sheetNames)) {
    for (i in 1:202) {
      overlap <- intersect(wikiPathways[[i]], wikiPathways[[sheetNames[j]]])
      if (length(overlap) > 0) {
        msg <- paste0(
          names(wikiPathways[i]), 
          " has ", 
          length(overlap),
          " matches with ",
          sheetNames[j]
        )
        messages <- c(messages, msg)
      }
    }
  }
  
  writeLines(messages, outFile)
  