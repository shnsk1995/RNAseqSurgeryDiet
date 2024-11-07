

PlotPathwayHeatmaps <- function(tissueOfInterest,pathways,pathwayNames){

  
  readPath <- "data/FGSEAResults/"
  
  if(!dir.exists("data/HigherLevelAnalysis")){
    dir.create("data/HigherLevelAnalysis")
  }
  
  if(!dir.exists("data/HigherLevelAnalysis/PathwayHeatmaps")){
    dir.create("data/HigherLevelAnalysis/PathwayHeatmaps")
  }
  
  writePath <- "data/HigherLevelAnalysis/PathwayHeatmaps/"
  
  for (i in 1:length(pathwayNames)) {
    
    if(!dir.exists(paste0("data/HigherLevelAnalysis/PathwayHeatmaps/",pathwayNames[i]))){
      dir.create(paste0("data/HigherLevelAnalysis/PathwayHeatmaps/",pathwayNames[i]))
    }
    
    writePath <- paste0("data/HigherLevelAnalysis/PathwayHeatmaps/",pathwayNames[i],"/")
    
    contrasts <- c("HF_Surgery","LF_Surgery","Diet","Surgery","Surgery_Diet")
    
    for (contrast in contrasts){
      
      pathNames <- names(pathways[[i]])
      
      heatmapData <- as.data.frame(matrix(ncol = length(tissueOfInterest),nrow = length(pathNames)))
      
      rownames(heatmapData) <- pathNames
      colnames(heatmapData) <- tissueOfInterest
      
      for (tissue in tissueOfInterest) {
        
        fgseares <- read.xlsx(paste0(readPath,tissue,"/AllPathways.xlsx"),sheet = paste0(pathwayNames[i],"_",contrast))
        
        sigUpPathwayNames <- fgseares[fgseares$padj < 0.05 & fgseares$NES > 0,]$pathway
        sigDownPathwayNames <- fgseares[fgseares$padj < 0.05 & fgseares$NES < 0,]$pathway
        
        for (pathway in pathNames) {
          
          if(pathway %in% sigUpPathwayNames){
            heatmapData[pathway,tissue] <- 1
          }else if(pathway %in% sigDownPathwayNames){
            heatmapData[pathway,tissue] <- -1
          }else{
            heatmapData[pathway,tissue] <- 0
          }
          
        }
        
        
      }
      
      heatmapData <- heatmapData %>%
        dplyr::mutate(UpCount = rowSums(across(all_of(tissueOfInterest)) == 1))
      
      heatmapData <- heatmapData %>%
        dplyr::mutate(DownCount = rowSums(across(all_of(tissueOfInterest)) == -1))
      
      heatmapData <- heatmapData %>%
        dplyr::mutate(NoSigCount = rowSums(across(all_of(tissueOfInterest)) == 0))
      
      heatmapDataDf <- heatmapData[!heatmapData$NoSigCount==8,-(9:11)]
      
      colFun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
      
      heatMap <- Heatmap(as.matrix(heatmapDataDf),
                         column_names_side = "top",
                         col = colFun,
                         border = TRUE,
                         row_names_gp = gpar(fontsize = 15),
                         column_names_gp = gpar(fontsize = 12),
                         heatmap_width = unit(28,"cm"),
                         heatmap_height = unit(40,"cm"),
                         show_heatmap_legend = FALSE,
                         rect_gp = gpar(col = "black", lwd = 1),
                         column_title =  paste0("Heatmap of significant pathways for ",contrast," contrast"),
                         column_title_gp = gpar(fontsize = 20, fontface = "bold")
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
      
      
      jpeg(paste0(writePath,pathwayNames[i],"_",contrast,"_PathwayHeatMap.jpeg"),height = 1400,width = 1600,quality = 100)
      draw(heatMap,annotation_legend_list = list(lgd),
           annotation_legend_side = "bottom"
      )
      dev.off()
      
      file <- paste0(writePath,pathwayNames[i],"_PathwaysSignificance.xlsx")
      
      if (file.exists(file)) {
        # Load the existing workbook
        wb <- loadWorkbook(file)
        
        
        sheetName <- contrast
        # Add a new sheet (check if the sheet already exists first)
        if (!(sheetName %in% names(wb))) {
          addWorksheet(wb, sheetName)
        }
        
        # Optionally, write some data to the new sheet
        writeData(wb, sheetName,heatmapData, startCol = 1, startRow = 1,rowNames = TRUE, colNames = TRUE)
        
        # Save the workbook
        saveWorkbook(wb, file, overwrite = TRUE)
        
      } else {
        
        # Create a new workbook if the file doesn't exist
        wb <- createWorkbook()
        
        # Add a new sheet
        sheetName <- contrast
        addWorksheet(wb, sheetName)
        
        # Optionally, write some data to the new sheet
        writeData(wb, sheetName,heatmapData, startCol = 1, startRow = 1,rowNames = TRUE, colNames = TRUE)
        
        # Save the new workbook
        saveWorkbook(wb, file, overwrite = TRUE)
      }
      
    }
    
    
    
  }
  
  
  
  
}