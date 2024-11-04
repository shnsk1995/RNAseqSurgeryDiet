

ConvertEntrezToEnsembl <- function(entrezIds, mapping) {
  ensemblIDs <- mapping[match(entrezIds, mapping$entrezgene_id), "ensembl_gene_id"]
  return(ensemblIDs)
  
}

CreateMouseGenomeHallMarkEnsemblPathways <- function(){
  
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") 
  
  # gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")
  gmt.file <- "data/mh.all.v2024.1.Mm.entrez.gmt"
  pathways <- gmtPathways(gmt.file)
  
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "entrezgene_id",
    values = unique(unlist(pathways)),
    mart = ensembl
  )
  
  # Apply the function to all pathways
  pathways <- lapply(pathways, ConvertEntrezToEnsembl, mapping = gene_mapping)
  
  return(pathways)
  
}

CreateMouseGenomeWikiEnsemblPathways <- function(){
  
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl") 
  
  # gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")
  gmt.file <- "data/wikipathways-20240910-gmt-Mus_musculus.gmt"
  pathways <- gmtPathways(gmt.file)
  
  gene_mapping <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "entrezgene_id",
    values = unique(unlist(pathways)),
    mart = ensembl
  )
  
  # Apply the function to all pathways
  pathways <- lapply(pathways, ConvertEntrezToEnsembl, mapping = gene_mapping)
  
  pathwayNames <- names(pathways)
  
  for (i in 1:length(pathwayNames)) {
    
    pathwayNames[i] <- sub("%.*","",pathwayNames[i])
    
  }
  
  names(pathways) <- pathwayNames
  
  return(pathways)
  
}

ConvertToSymbolAll <- function(ensemblIds, geneInfo) {
  ensemblList <- unlist(strsplit(ensemblIds, ", "))
  
  
  
  symbols <- sapply(ensemblList, function(ensemblId) {
    symbol <- geneInfo$gene_name[geneInfo$ensembl_gene_id == ensemblId]
    
    
    if (length(symbol) == 0 || is.na(symbol)) {
      return(ensemblId)
    } else {
      return(symbol)
    }
  })
  
  
  return(paste(symbols, collapse = ","))
}

ConvertToSymbolMain <- function(ensemblIds, geneInfo) {
  ensemblList <- unlist(strsplit(ensemblIds, ",, "))
  
  
  
  symbols <- sapply(ensemblList, function(ensemblId) {
    symbol <- geneInfo$gene_name[geneInfo$ensembl_gene_id == ensemblId]
    
    
    if (length(symbol) == 0 || is.na(symbol)) {
      return(ensemblId)
    } else {
      return(symbol)
    }
  })
  
  
  return(paste(symbols, collapse = ","))
}

ConvertEnsemblToSymbol <- function(tissueOfInterest,geneInfo){
  
  for (tissue in tissueOfInterest) {
    
    dirPath <- paste0("data/FGSEAResults/",tissue,"/")
    
    wbAll <- loadWorkbook(paste0(dirPath,"AllPathways.xlsx"))
    
    wbMain <- loadWorkbook(paste0(dirPath,"MainPathways.xlsx"))
    
    allSheets <- names(wbAll)
    
    mainSheets <- names(wbMain)
    
    for (sheet in allSheets) {
      
      fgseaRes <- read.xlsx(wbAll, sheet = sheet)
      
      fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, ConvertToSymbolAll, geneInfo = geneInfo)
      
      writeData(wbAll, sheet = sheet, x = fgseaRes)
      
    }
    
    saveWorkbook(wbAll, paste0(dirPath,"AllPathways.xlsx"), overwrite = TRUE)
    
    for (sheet in mainSheets) {
      
      fgseaRes <- read.xlsx(wbMain, sheet = sheet)
      
      fgseaRes$leadingEdge <- sapply(fgseaRes$leadingEdge, ConvertToSymbolMain, geneInfo = geneInfo)
      
      writeData(wbMain, sheet = sheet, x = fgseaRes)
      
    }
    
    saveWorkbook(wbMain, paste0(dirPath,"MainPathways.xlsx"), overwrite = TRUE)
    
    
  }
  
  
  
  
}

DoFGSEAForAllGenes <- function(tissue,pathways=list(),pathwayNames=c()){
  
  
  
  if(!dir.exists("data/FGSEAResults/")){
    dir.create("data/FGSEAResults/")
  }
  
  readDirPath <- paste0("data/LimmaResults/",tissue,"/")
  writeDirPath <- paste0("data/FGSEAResults/",tissue,"/")
  
  if(!dir.exists(writeDirPath)){
    dir.create(writeDirPath)
  }
  
  contrasts <- readLines(paste0(readDirPath,"ContrastsUsed.txt"))
  
  for (i in 1:length(contrasts)) {
    
    contrasts[i] <- sub(" =.*","",contrasts[i])
    
  }
  
  for (i in 1:length(pathwayNames)) {
    
    for (contrast in contrasts) {
      
      allGenes <- read.xlsx(paste0(readDirPath,"AllGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
      
      allGenes <- allGenes %>%
        arrange(desc(logFC))
      
      generanks <- setNames(
        allGenes$logFC,
        allGenes$EnsemblID      
      )
      
      fgseaRes <- fgsea(pathways[[i]], generanks, minSize=15, maxSize=250)
      
      #Save FGSEA results
      #fwrite(fgseaRes, file=paste0(dirPath,pathwayNames[i],"_",contrast,"_FGSEAResults.txt"), sep="\t", sep2=c("", " ", ""))
      file <- paste0(writeDirPath,"/AllPathways.xlsx")
      
      if (file.exists(file)) {
        # Load the existing workbook
        wb <- loadWorkbook(file)
        
        sheetName <- paste0(pathwayNames[i],"_",contrast)
        
        # Add a new sheet (check if the sheet already exists first)
        if (!(sheetName %in% names(wb))) {
          addWorksheet(wb, sheetName)
        }
        
        # Optionally, write some data to the new sheet
        writeData(wb, sheetName,fgseaRes, startCol = 1, startRow = 1)
        
        # Save the workbook
        saveWorkbook(wb, file, overwrite = TRUE)
        
      } else {
        
        # Create a new workbook if the file doesn't exist
        wb <- createWorkbook()
        
        
        # Add a new sheet
        sheetName <- paste0(pathwayNames[i],"_",contrast)
        addWorksheet(wb, sheetName)
        
        # Optionally, write some data to the new sheet
        writeData(wb, sheetName,fgseaRes, startCol = 1, startRow = 1)
        
        # Save the new workbook
        saveWorkbook(wb, file, overwrite = TRUE)
      }
      
      
    }
    
    
  }
  
}

SaveAndPlotTopPathways <- function(tissue,pathways=list(),pathwayNames=c()){
  
  
  readDirPath <- paste0("data/LimmaResults/",tissue,"/")
  writeDirPath <- paste0("data/FGSEAResults/",tissue,"/")
  
  contrasts <- readLines(paste0(readDirPath,"ContrastsUsed.txt"))
  
  for (i in 1:length(contrasts)) {
    
    contrasts[i] <- sub(" =.*","",contrasts[i])
    
  }
  
  for (i in 1:length(pathwayNames)) {
    
    plotsPath <- paste0(writeDirPath,"/",pathwayNames[i],"_Plots/")
    
    if(!dir.exists(plotsPath)) dir.create(plotsPath)
    
    for (contrast in contrasts) {
      
      allGenes <- read.xlsx(paste0(readDirPath,"AllGenes.xlsx"), sheet = contrast, rowNames = TRUE, colNames = TRUE)
      
      allGenes <- allGenes %>%
        arrange(desc(logFC))
      
      generanks <- setNames(
        allGenes$logFC,
        allGenes$EnsemblID      
      )
      
      sheetName <- paste0(pathwayNames[i],"_",contrast)
      
      fgseaRes <- read.xlsx(paste0(writeDirPath,"AllPathways.xlsx"), sheet = sheetName, colNames = TRUE)
      fgseaRes <- as.data.table(fgseaRes)
      #fgseaRes <- fread(paste0(dirPath,pathwayNames[i],"_",contrast,"_FGSEAResults.txt"), sep = "\t")
      
      
      topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
      topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
      topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

      jpeg(filename = paste0(plotsPath,pathwayNames[i],"_",contrast,"_TopPathways.jpeg"),height = 1200 , width = 1000,quality = 100)


      gSEATable <- plotGseaTable(pathways[[i]][topPathways], generanks, fgseaRes,
                                 gseaParam=0.5)


      print(gSEATable)

      dev.off()


      collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
                                            pathways[[i]], generanks)

      mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
        order(-NES), pathway]

      jpeg(filename = paste0(plotsPath,pathwayNames[i],"_",contrast,"_MainPathways.jpeg"),height = 1200 , width = 1000,quality = 100)

      gSEATable <- plotGseaTable(pathways[[i]][mainPathways], generanks, fgseaRes,
                                 gseaParam = 0.5)

      print(gSEATable)

      dev.off()


      fgseaResMain <- fgseaRes[match(mainPathways, pathway)]
      
      fgseaResMain[, leadingEdge := strsplit(leadingEdge, split = " ")]

      # fgseaResMain[, leadingEdge := mapIdsList(
      #   x=org.Mm.eg.db,
      #   keys=leadingEdge,
      #   keytype="ENSEMBL",
      #   column="SYMBOL")]


      #fwrite(fgseaResMain, file=paste0(dirPath,pathwayNames[i],"_",contrast,"_FGSEAResMain.txt"), sep="\t", sep2=c("", " ", ""))
      
      #Save FGSEA results
      #fwrite(fgseaRes, file=paste0(dirPath,pathwayNames[i],"_",contrast,"_FGSEAResults.txt"), sep="\t", sep2=c("", " ", ""))
      file <- paste0(writeDirPath,"/MainPathways.xlsx")
      
      if (file.exists(file)) {
        # Load the existing workbook
        wb <- loadWorkbook(file)
        
        
        sheetName <- paste0(pathwayNames[i],"_",contrast)
        # Add a new sheet (check if the sheet already exists first)
        if (!(sheetName %in% names(wb))) {
          addWorksheet(wb, sheetName)
        }
        
        # Optionally, write some data to the new sheet
        writeData(wb, sheetName,fgseaResMain, startCol = 1, startRow = 1)
        
        # Save the workbook
        saveWorkbook(wb, file, overwrite = TRUE)
        
      } else {
        
        # Create a new workbook if the file doesn't exist
        wb <- createWorkbook()
        
        # Add a new sheet
        sheetName <- paste0(pathwayNames[i],"_",contrast)
        addWorksheet(wb, sheetName)
        
        # Optionally, write some data to the new sheet
        writeData(wb, sheetName,fgseaResMain, startCol = 1, startRow = 1)
        
        # Save the new workbook
        saveWorkbook(wb, file, overwrite = TRUE)
      }


      #Enrichment plots for top hits
      topHits <- fgseaRes[padj<0.05][order(padj), pathway]

      rl <- generanks
      setList <- pathways[[i]][topHits]
      gsea <- fgsea::fgsea(setList, rl, nperm=1000)

      #generate data for enrichment plots
      dataForEnrichPlots <- gseaCurve(rl,setList,gsea)


      #Generate pathway image
      enrichmentPlot <- ggplot2::ggplot() +
        geom_gsea(
          dataForEnrichPlots,
          linecolor = "purple",
          zeroline = TRUE,
          linesize = 1,
          ncol = 5,
        ) +
        ggtitle(paste0("Enrichment plot for top hits(padj < 0.05)")) +
        theme(plot.title = element_text(hjust = 0.5)) +
        theme_gsea(textsize = 10)
      
      
      plotHeight <- ceiling(length(topHits)/5) * 4
      if(plotHeight==4){
        plotWidth <- 4 * length(topHits)
      }else{
        plotWidth <- 20
      }

      ggsave(paste0(plotsPath,pathwayNames[i],"_",contrast,"_EnrichmentPlots.jpeg"),width = plotWidth ,height = plotHeight ,dpi=588)
      #ggsave(paste0(plotsPath,pathwayNames[i],"_",contrast,"_EnrichmentPlots.pdf"), width=12, height=12)
      
      
    }
    
    
  }
  
  
  
}

VennForPathways <- function(tissueOfInterest,pathways,pathwayNames){
  
  readPath <- "data/FGSEAResults/"
  
  
  
  for (p in 1:length(pathwayNames)) {
    
    
    for (tissue in tissueOfInterest) {
      
      writePath <- paste0(readPath,tissue,"/",pathwayNames[p],"_Plots/")
      
      contrasts <- c("HF_Surgery","LF_Surgery","Diet","Surgery")
      
      pathNames <- names(pathways[[p]])
      
      commonPathways <- as.data.frame(matrix(ncol = length(contrasts), nrow = length(pathNames)))
      rownames(commonPathways) <- pathNames
      colnames(commonPathways) <- contrasts
      
      for (contrast in contrasts) {
        
        fgseares <- read.xlsx(paste0(readPath,tissue,"/AllPathways.xlsx"),sheet = paste0(pathwayNames[p],"_",contrast))
        
        sigUpPathwayNames <- fgseares[fgseares$padj < 0.05 & fgseares$NES > 0,]$pathway
        sigDownPathwayNames <- fgseares[fgseares$padj < 0.05 & fgseares$NES < 0,]$pathway
        
        for (pathway in pathNames) {
          
          if(pathway %in% sigUpPathwayNames){
            commonPathways[pathway,contrast] <- 1
          }else if(pathway %in% sigDownPathwayNames){
            commonPathways[pathway,contrast] <- -1
          }else{
            commonPathways[pathway,contrast] <- 0
          }
          
        }
        
      }
      
      #Use function GenerateSignificanceMatrix()
      #write.csv(commonPathways,paste0(writePath,pathwayNames[p],"_SignificanceMatrix.csv"),row.names = TRUE)
      
      
      vennCounts <- vennCounts(commonPathways,include = "up")
      vennCounts <- as.data.frame(as.table(vennCounts))
      vennCounts <- pivot_wider(vennCounts,names_from = Var2, values_from = Freq)
      vennCounts <- vennCounts[,-1]
      
      #write.csv(vennCounts,paste0(writePath,pathwayNames[p],"_VennCountsUp.csv"),row.names = FALSE)
      
      vennCounts <- vennCounts[vennCounts$Counts!=0,]
      
      vennCounts <- vennCounts[!(vennCounts$HF_Surgery == 0 &
                                   vennCounts$LF_Surgery == 0 &
                                   vennCounts$Diet == 0 &
                                   vennCounts$Surgery == 0),]
      
      for (i in 1:nrow(vennCounts)) {
        
        selectedContrasts <- colnames(vennCounts)[which(vennCounts[i,(1:4)] == 1)]
        nonSelectedContrasts <- setdiff(colnames(commonPathways), selectedContrasts)
        
        filteredPathways <- commonPathways[
          rowSums(commonPathways[, selectedContrasts, drop = FALSE] == 1) == length(selectedContrasts) &  # All selected contrasts should be 1
            rowSums(commonPathways[, nonSelectedContrasts, drop = FALSE] == 1) == 0,  # No other columns should be 1
        ]
        
        mutualPathways <- rownames(filteredPathways)
        
        
        if(!file.exists(paste0(writePath,pathwayNames[p],"_UpVennPathways.xlsx"))){
          # Create a new workbook if the file doesn't exist
          wb <- createWorkbook()
          
          # Add a new sheet
          firstCharacters <- substr(selectedContrasts,1,1)
          sheetName <- paste(firstCharacters,collapse = "_")
          addWorksheet(wb, sheetName)
          
          # Optionally, write some data to the new sheet
          writeData(wb, sheetName,mutualPathways, startCol = 1, startRow = 1)
          
          # Save the new workbook
          saveWorkbook(wb, paste0(writePath,pathwayNames[p],"_UpVennPathways.xlsx"), overwrite = TRUE)
        }else{
          
          # Load the existing workbook
          wb <- loadWorkbook(paste0(writePath,pathwayNames[p],"_UpVennPathways.xlsx"))
          
          # Add a new sheet (check if the sheet already exists first)
          firstCharacters <- substr(selectedContrasts,1,1)
          sheetName <- paste(firstCharacters,collapse = "_")
          if (!(sheetName %in% names(wb))) {
            addWorksheet(wb, sheetName)
          }
          
          # Optionally, write some data to the new sheet
          writeData(wb, sheetName,mutualPathways, startCol = 1, startRow = 1)
          
          # Save the workbook
          saveWorkbook(wb, paste0(writePath,pathwayNames[p],"_UpVennPathways.xlsx"), overwrite = TRUE)
          
        }
        
      }
      
      vennCounts <- vennCounts(commonPathways,include = "down")
      vennCounts <- as.data.frame(as.table(vennCounts))
      vennCounts <- pivot_wider(vennCounts,names_from = Var2, values_from = Freq)
      vennCounts <- vennCounts[,-1]
      
      #write.csv(vennCounts,paste0(writePath,pathwayNames[p],"_VennCountsDown.csv"),row.names = FALSE)
      
      vennCounts <- vennCounts[vennCounts$Counts!=0,]
      
      vennCounts <- vennCounts[!(vennCounts$HF_Surgery == 0 &
                                   vennCounts$LF_Surgery == 0 &
                                   vennCounts$Diet == 0 &
                                   vennCounts$Surgery == 0),]
      
      for (i in 1:nrow(vennCounts)) {
        
        selectedContrasts <- colnames(vennCounts)[which(vennCounts[i,(1:4)] == 1)]
        nonSelectedContrasts <- setdiff(colnames(commonPathways), selectedContrasts)
        
        filteredPathways <- commonPathways[
          rowSums(commonPathways[, selectedContrasts, drop = FALSE] == -1) == length(selectedContrasts) &  
            rowSums(commonPathways[, nonSelectedContrasts, drop = FALSE] == -1) == 0, 
        ]
        
        mutualPathways <- rownames(filteredPathways)
        
        
        if(!file.exists(paste0(writePath,pathwayNames[p],"_DownVennPathways.xlsx"))){
          # Create a new workbook if the file doesn't exist
          wb <- createWorkbook()
          
          # Add a new sheet
          firstCharacters <- substr(selectedContrasts,1,1)
          sheetName <- paste(firstCharacters,collapse = "_")
          addWorksheet(wb, sheetName)
          
          # Optionally, write some data to the new sheet
          writeData(wb, sheetName,mutualPathways, startCol = 1, startRow = 1)
          
          # Save the new workbook
          saveWorkbook(wb, paste0(writePath,pathwayNames[p],"_DownVennPathways.xlsx"), overwrite = TRUE)
        }else{
          
          # Load the existing workbook
          wb <- loadWorkbook(paste0(writePath,pathwayNames[p],"_DownVennPathways.xlsx"))
          
          # Add a new sheet (check if the sheet already exists first)
          firstCharacters <- substr(selectedContrasts,1,1)
          sheetName <- paste(firstCharacters,collapse = "_")
          if (!(sheetName %in% names(wb))) {
            addWorksheet(wb, sheetName)
          }
          
          # Optionally, write some data to the new sheet
          writeData(wb, sheetName,mutualPathways, startCol = 1, startRow = 1)
          
          # Save the workbook
          saveWorkbook(wb, paste0(writePath,pathwayNames[p],"_DownVennPathways.xlsx"), overwrite = TRUE)
          
        }
        
      }
      
      baseColors <- c("red","blue","green")
      
      custom_palette <- colorRampPalette(baseColors)
      
      colors <- custom_palette(length(contrasts))
      
      
      jpeg(filename = paste0(writePath,pathwayNames[p],"_UpVenn.jpeg"),width = 1000,height = 800,quality = 100)
      
      
      vennDiagram(commonPathways,include = "up",
                  lwd = 1,
                  circle.col = colors
      )
      
      
      dev.off()
      
      
      
      jpeg(filename = paste0(writePath,pathwayNames[p],"_DownVenn.jpeg"),width = 1000,height = 800,quality = 100)
      
      
      vennDiagram(commonPathways,include = "down",
                  lwd = 1,
                  circle.col = colors
      )
      
      
      dev.off()
      
      
      
      
    }
    
    
    
    
  }
  
  
  
  
}



GenerateSignificanceMatrix <- function(tissueOfInterest,pathways,pathwayNames){
  
  readPath <- "data/FGSEAResults/"
  
  for (p in 1:length(pathwayNames)) {
    
    
    for (tissue in tissueOfInterest) {
      
      writePath <- paste0(readPath,tissue,"/",pathwayNames[p],"_Plots/")
      
      contrasts <- c("HF_Surgery","LF_Surgery","Diet","Surgery","Surgery_Diet")
      
      pathNames <- names(pathways[[p]])
      
      commonPathways <- as.data.frame(matrix(ncol = length(contrasts), nrow = length(pathNames)))
      rownames(commonPathways) <- pathNames
      colnames(commonPathways) <- contrasts
      
      for (contrast in contrasts) {
        
        fgseares <- read.xlsx(paste0(readPath,tissue,"/AllPathways.xlsx"),sheet = paste0(pathwayNames[p],"_",contrast))
        
        sigUpPathwayNames <- fgseares[fgseares$padj < 0.05 & fgseares$NES > 0,]$pathway
        sigDownPathwayNames <- fgseares[fgseares$padj < 0.05 & fgseares$NES < 0,]$pathway
        
        for (pathway in pathNames) {
          
          if(pathway %in% sigUpPathwayNames){
            commonPathways[pathway,contrast] <- 1
          }else if(pathway %in% sigDownPathwayNames){
            commonPathways[pathway,contrast] <- -1
          }else{
            commonPathways[pathway,contrast] <- 0
          }
          
        }
        
      }
      
      write.csv(commonPathways,paste0(writePath,pathwayNames[p],"_SignificanceMatrix.csv"),row.names = TRUE)
      
    }
    
  }
  
}


