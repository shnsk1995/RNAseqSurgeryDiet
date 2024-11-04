source("Libraries.R")
source("Functions.R")
source("DGEFunctions.R")
source("FGSEAFunctions.R")
source("HigherLevelAnalysis.R")


comprehensiveData <- comprehensiveData %>%
  mutate(TotalOccurence = rowSums(.[,1:32]))

comprehensiveData <- comprehensiveData[order(comprehensiveData$TotalOccurence,decreasing = TRUE),]

comprehensiveData <- head(comprehensiveData,n=50)

comprehensiveData <- as.matrix(comprehensiveData)

hmdata <- comprehensiveData[,-33]

jpeg(filename = "data/test.jpeg",height = 1400,width = 1000, quality = 100)

Heatmap(hmdata,
        row_names_side = "left",
        column_names_side = "top",
        col = colorRamp2(c(min(hmdata),max(hmdata)),c("white","red")),
        heatmap_legend_param = list(title = "Significance"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        heatmap_width = unit(25,"cm"),
        heatmap_height = unit(40,"cm"),
        border = TRUE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x, y, width, height, 
                    gp = gpar(col = "black", fill = NA))  # Draw cell borders
        }
        )

lgd = Legend(
  labels = c("Sig", "NotSig"),
  title = "Significance",
  legend_gp = gpar(fill = c("red", "white")),
  border = TRUE,
  direction = "horizontal", 
  title_gp = gpar(fontsize = 14),
  labels_gp = gpar(fontsize = 12)
)

dev.off()
