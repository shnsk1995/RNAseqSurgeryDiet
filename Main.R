source("Libraries.R")
source("Functions.R")
source("DGEFunctions.R")
source("FGSEAFunctions.R")
source("HigherLevelAnalysis.R")


#Read data
cts <- as.data.frame(read.csv("data/counts.csv",row.names = 1))
sampleInfo <- as.data.frame(read.csv("data/SampleInfo.csv",row.names = 1))
geneInfo <- as.data.frame(read.csv("data/geneInfo.csv",row.names = 1))


#Specify the tissue of interest
tissueOfInterest <- unique(sampleInfo$tissuedetail)

#Specify any samples to be removed
samplesToRemove <- c("LF_C_I_R1_B1","LF_C_I_R4_B1","LF_S_II_R1_B1",
                     "HF_S_II_R1_B1","LF_S_J_R2_B1","LF_C_J_R1_B1",
                     "HF_C_J_R1_B2","LF_S_L_R3_B1","LF_S_SBr_R2_B1",
                     "LF_C_SBr_R3_B1","LF_S_SW_R2_B1","HF_S_SW_R1_B1")

samplesToRemove <- c()

#Do PCA for each tissue
for (tissue in tissueOfInterest) {
  
  DOPCA(cts = cts,sampleInfo = sampleInfo, geneInfo = geneInfo, tissueOfInterest = tissue, samplesToExclude = samplesToRemove)
  
}


#Save samples to be removed
samplesFiltered <- as.data.frame(samplesToRemove)
write.csv(samplesFiltered,file = "data/Filtered_samples.csv",row.names = FALSE)


#Do correlation plots for each sample group
sampleGroups <- unique(sampleInfo$samplename)

for (sampleGroup in sampleGroups) {
  
  DoGroupCorrelationPlots(cts,sampleInfo,geneInfo,sampleGroup)
  
}


#Do comprehensive correlation plot
DoComprehensiveCorrelationPlot(cts,sampleInfo,geneInfo)



#Do Cross Correlation plot as required
sampleGroup1 <- "HF_C_D"
sampleGroup2 <- "LF_S_D"
DoCrossCorrelationPlot(cts,sampleInfo, geneInfo, sampleGroup1,sampleGroup2)



#Do correlation plots across tissues and conditions
lowQualitySamples <- FilterLowQualitySamples(cts,sampleInfo,geneInfo)

#Save filtered samples and counts
ctsFiltered <- cts[,!colnames(cts) %in% lowQualitySamples]
sampleInfoFiltered <- sampleInfo[!sampleInfo$sampleid %in% lowQualitySamples,]
#Write the counts data to an excel file
write.csv(ctsFiltered,file = "data/Filtered_counts.csv",row.names = TRUE)
#Write the filtered sample info to an excel file
write.csv(sampleInfoFiltered,file  = "data/Filtered_sampleInfo.csv",row.names = TRUE)



#Read filtered counts and sample information
cts <- as.data.frame(read.csv("data/Filtered_counts.csv",row.names = 1))
sampleInfo <- as.data.frame(read.csv("data/Filtered_sampleInfo.csv",row.names = 1))
geneInfo <- as.data.frame(read.csv("data/geneInfoCorrected.csv",row.names = 1))


#Generate density and library size plots for all tissues
Generate_Density_N_Library_Plot_Grid(allCts = cts, allSampleInfo = sampleInfo, geneInfo = geneInfo)


#Specify the tissue of interest
tissueOfInterest <- unique(c(unique(sampleInfo$tissuedetail),unique(sampleInfo$tissuegeneral)))

tissueOfInterest <- tissueOfInterest[tissueOfInterest!="Ileum interposition"]

for (tissue in tissueOfInterest) {
  
    
    # Compare_HFS_With_HFC(cts,sampleInfo,geneInfo,tissue)
    # Compare_LFS_With_LFC(cts,sampleInfo,geneInfo,tissue)
    # Compare_HF_With_LF(cts,sampleInfo,geneInfo,tissue)
    # Compare_Surgery_With_Diet(cts,sampleInfo,geneInfo,tissue)
    # Compare_Surgery_Without_Diet(cts,sampleInfo,geneInfo,tissue)
    Compare_MultipleContrasts(cts,sampleInfo,geneInfo,tissue)

  
}


#Extract Summary of DE genes for each tissue and comparison
SummarizeDEGenesFDR(tissueOfInterest,0.05,1)
SummarizeDEGenesPValue(tissueOfInterest,0.001,1)


#Do FGSEA
hallMarkPathways <- CreateMouseGenomeHallMarkEnsemblPathways()
wikiPathways <- CreateMouseGenomeWikiEnsemblPathways()
pathwaysList <- list(hallMarkPathways,wikiPathways)
pathwaysListNames <- c("HallMark","Wiki")

for (tissue in tissueOfInterest) {
  
  DoFGSEAForAllGenes(tissue = tissue,pathways = pathwaysList,pathwayNames = pathwaysListNames)
  
}

#Extract Summary of pathways for each tissue and comparison
SummarizePathwayFDR(tissueOfInterest,0.05,pathwayNames = pathwaysListNames)
SummarizePathwayPValue(tissueOfInterest,0.001,pathwayNames = pathwaysListNames)

#Save and plot top pathways
for (tissue in tissueOfInterest) {
  
  SaveAndPlotTopPathways(tissue = tissue,pathways = pathwaysList,pathwayNames = pathwaysListNames)
  
}

#Generate box plots for top 5 genes in Surgery_Diet contrast
PlotBoxPlotsForTOPGenes(tissueOfInterest,c("Surgery_Diet"),5,"pvalue",1,sampleInfo)



#Higher level analysis for pathways
PlotPathwayHeatmaps(tissueOfInterest,pathwaysList,pathwaysListNames)


#Venn for pathways
VennForPathways(tissueOfInterest,pathwaysList,pathwaysListNames)


#Convert EnsemblIds to symbols in pathway result files
ConvertEnsemblToSymbol(tissueOfInterest,geneInfo)


#Generate heatmaps for pathway summary files
inputFiles <- list.files("data/FGSEAResults/",pattern = "*.csv",full.names = TRUE)
outputFiles <- sub("*.csv","Heatmap.jpeg",inputFiles)

GenerateHeatmapForSummaryValues(inputFiles = inputFiles,OutputFiles = outputFiles)



#Generate bar plots for pathways
contrasts <- c("HF_Surgery","LF_Surgery","Diet","Surgery","Surgery_Diet")
GenerateBarPlotsForPathways(tissueOfInterest,pathwaysListNames,contrasts)



#Generate significance matrix and heatmaps for pathways
GenerateSignificanceMatrix(tissueOfInterest,pathwaysList,pathwaysListNames)
GenerateSignificanceHeatmapForPathways(tissueOfInterest,pathwaysListNames)


#Convert VennGeneList to significanceMatrix
ConvertVennGeneListToSM(tissueOfInterest)


#Append relevant files to AllGenes
AppendRelevantFilesToAllGenes(tissueOfInterest)

#Append relevant files to AllPathways
AppendRelevantFilesToAllPathways(tissueOfInterest,pathwaysListNames)

#Add notes to Gene files
AddNotesToGeneFiles(tissueOfInterest = tissueOfInterest)

#Add notes to Pathway files
AddNotesToPathwayFiles(tissueOfInterest = tissueOfInterest)

#Correct Heatmaps - One time use - Not required if you are rerunning the entire analysis (11/01/2024)
#CorrectHeatmaps(tissueOfInterest)


#Delete files works in data/Share folder only
#Delete unnecessary gene files
DeleteUnnecessaryGeneFiles(tissueOfInterest = tissueOfInterest)

#Delete unnecessary pathway files
DeleteUnnecessaryPathwayFiles(tissueOfInterest = tissueOfInterest)


#Do Gene Abundance Analysis
resultsList <- DoComprehensiveGeneAbundanceAnalysis(tissueOfInterest = tissueOfInterest)
comprehensiveData <- resultsList[[1]]
tissueData <- resultsList[[2]]

#Append respective tissue names to all files
# AppendTissueNamesToAllFiles(tissueOfInterest)

