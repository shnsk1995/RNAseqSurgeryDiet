source("Libraries.R")

#Set directory path
projectPath <- "C:/Users/SuryadHN/PBRC BioInfo Core/PBRC_Labs/Kachmar/"


#Read sample batch Info
sampleBatchInfo <- read.csv("data/Cohort Samples by Time.csv")

#Extract all samples
allSamples <- as.data.frame(sampleBatchInfo[,-4])
colnames(allSamples) <- c("Spring2023","Spring2024","Spring20232")
allSamples <- allSamples[-(1:2),]
allSamples <- c(allSamples$Spring2023,allSamples$Spring2024,allSamples$Spring20232)
allSamples <- allSamples[allSamples != ""]
allSamplesFreq <- as.data.frame(table(allSamples))

#Extract sample time info
sampleInfoTime <- as.data.frame(sampleBatchInfo[,-4])
colnames(sampleInfoTime) <- c("Spring2023","Spring2024","Spring20232")
sampleInfoTime <- sampleInfoTime[-(1:2),]

#Extract sample diet info
sampleInfoDiet <- as.data.frame(sampleBatchInfo[,-4])
sampleInfoDiet <- sampleInfoDiet[-(1:2),]

#Extract samples of high diet
highFatDietSamples <- c(sampleInfoDiet$High.Fat.Cohort,sampleInfoDiet$High.Fat.Cohort.1)
highFatDietSamples <- highFatDietSamples[highFatDietSamples != ""]

#Extact samples of low diet
lowFatDietSamples <- c(sampleInfoDiet$LF.Cohorts)
lowFatDietSamples <- lowFatDietSamples[lowFatDietSamples != ""]


#Extract samples of batch 1
batch1Samples <- c(sampleInfoTime$`Spring2023`,sampleInfoTime$`Spring20232`)
batch1Samples <- batch1Samples[batch1Samples!=""]

#Extract samples of batch 2
batch2Samples <- c(sampleInfoTime$`Spring2024`)
batch2Samples <- batch2Samples[batch2Samples!=""]


#Function to extract batch number
GetSampleBatch <- function(title){
  
  if(title %in% batch1Samples) return(1)
  if(title %in% batch2Samples) return(2)
  return(NA)
  
}


#Function to extract diet
GetSampleDiet <- function(title){
  
  if(title %in% highFatDietSamples) return("High fat")
  if(title %in% lowFatDietSamples) return("Low fat")
  return(NA)
  
}


#Function to extract treatment info
GetTreatmentInfo <- function(title){
  
  if(length(grep("II",title))==1) return("Surgery")
  if(length(grep("C",title))==1) return("Control")
  return(NA)
  
  
}


#Function to extract detailed tissue info
GetDetailTissueInfo <- function(title){
  
  sampleData <- sampleInfo[sampleInfo$title==title,]
  
  if(sampleData$tissue == "Small Bowel Mucosa" && (length(grep("_I$",title))==1)) return("Ileum")
  if(sampleData$tissue == "Small Bowel Mucosa" && (length(grep("_II$",title))==1)) return("Ileum interposition")
  if(sampleData$tissue == "Small Bowel Mucosa" && length(grep("_D$",title))==1) return("Duodenum")
  if(sampleData$tissue == "Small Bowel Mucosa" && length(grep("_J$",title))==1) return("Jejunum")
  
  if(sampleData$tissue == "Adipose" && length(grep("_eWAT$",title))==1) return("Epididymal white adipose")
  if(sampleData$tissue == "Adipose" && length(grep("_scWAT$",title))==1) return("Subcutaneous white adipose")
  if(sampleData$tissue == "Adipose" && length(grep("_scBAT$",title))==1) return("Subcutaneous brown adipose")
  
  if(sampleData$tissue == "Hepatic") return("Liver")
  
  return(title)
  
}


#Function to construct a sample name
ConstructSampleName <- function(title){
  
  sampleData <- sampleInfo[sampleInfo$title==title,]
  
  #SampleDiet
  if(is.na(sampleData$diet)) {
    sampleDiet <- "NA"
  }else if(sampleData$diet=="High fat"){
    sampleDiet <- "HF"
  }else if(sampleData$diet=="Low fat"){
    sampleDiet <- "LF"
  }
  
  #Treatment
  if(is.na(sampleData$treatment)){
    sampleTreatment <- "N"
  }else if(sampleData$treatment=="Surgery"){
    sampleTreatment <- "S"
  }else if (sampleData$treatment=="Control"){
    sampleTreatment <- "C"
  }
  
  
  #Tissue
  if(sampleData$tissuedetail == "Ileum"){
    sampleTissue <- "I"
  }else if(sampleData$tissuedetail == "Ileum interposition"){
    sampleTissue <- "II"
  }else if(sampleData$tissuedetail == "Duodenum"){
    sampleTissue <- "D"
  }else if(sampleData$tissuedetail == "Jejunum") {
    sampleTissue <- "J"
  } else if(sampleData$tissuedetail == "Liver") {
    sampleTissue <- "L"
  }else if(sampleData$tissuedetail == "Epididymal white adipose"){
    sampleTissue <- "EA"
  }else if(sampleData$tissuedetail == "Subcutaneous white adipose"){
    sampleTissue <- "SW"
  }else if(sampleData$tissuedetail == "Subcutaneous brown adipose"){
    sampleTissue <- "SBr"
  }else{
    sampleTissue <- "N"
  }
  
  
  return(paste0(sampleDiet,"_",sampleTreatment,"_",sampleTissue))
  
}

#Function to construct a general sample name
ConstructSampleNameGeneral <- function(title){
  
  sampleData <- sampleInfo[sampleInfo$title==title,]
  
  #SampleDiet
  if(is.na(sampleData$diet)) {
    sampleDiet <- "NA"
  }else if(sampleData$diet=="High fat"){
    sampleDiet <- "HF"
  }else if(sampleData$diet=="Low fat"){
    sampleDiet <- "LF"
  }
  
  #Treatment
  if(is.na(sampleData$treatment)){
    sampleTreatment <- "N"
  }else if(sampleData$treatment=="Surgery"){
    sampleTreatment <- "S"
  }else if (sampleData$treatment=="Control"){
    sampleTreatment <- "C"
  }
  
  
  #Tissue
  if(sampleData$tissuegeneral == "Small Bowel Mucosa"){
    sampleTissue <- "SB"
  } else if(sampleData$tissuegeneral == "Liver") {
    sampleTissue <- "L"
  }else if(sampleData$tissuegeneral == "Epididymal white adipose"){
    sampleTissue <- "EA"
  }else if(sampleData$tissuegeneral == "Subcutaneous white adipose"){
    sampleTissue <- "SW"
  }else if(sampleData$tissuegeneral == "Subcutaneous brown adipose"){
    sampleTissue <- "SBr"
  }else{
    sampleTissue <- "N"
  }
  
  
  return(paste0(sampleDiet,"_",sampleTreatment,"_",sampleTissue))
  
}


#Function to construct a sample id
ConstructSampleId <- function(samplename,title){
  
  sampleData <- sampleInfo[sampleInfo$title==title,]
  
  #Batch
  if(is.na(sampleData$batch)){
    sampleBatch <- 0
  }else if(sampleData$batch==1){
    sampleBatch <- 1
  }else if (sampleData$batch==2){
    sampleBatch <- 2
  }
  
  #Replicate
  sampleReplicate <- sampleData$replicate
  
  return(paste0(samplename,"_R",sampleReplicate,"_B",sampleBatch))
  
}

#Function to construct a general sample id
ConstructSampleIdGeneral <- function(samplename,title){
  
  sampleData <- sampleInfo[sampleInfo$title==title,]
  
  #Batch
  if(is.na(sampleData$batch)){
    sampleBatch <- 0
  }else if(sampleData$batch==1){
    sampleBatch <- 1
  }else if (sampleData$batch==2){
    sampleBatch <- 2
  }
  
  #Replicate
  sampleReplicate <- sampleData$replicate.general
  
  return(paste0(samplename,"_R",sampleReplicate,"_B",sampleBatch))
  
}



#Read sample info
sampleInfo <- data.frame(readxl::read_xlsx("data/Kachmar_SAGES_geo.xlsx",sheet = "2. Metadata Template",range = "A30:P189"))

#Remove columns with all NAs
sampleInfo <- sampleInfo[,colSums(!is.na(sampleInfo))>0]

#Obtain the frequency of all the samples
samcts <- as.data.frame(table(sampleInfo$title))


#Check if there are any repetitions in samples
sum(samcts$Freq)
sum(allSamplesFreq$Freq)

#Check if all samples in time file are present in metadata
all(samcts$Var1 %in% allSamplesFreq$allSamples)

#Extract batch number, diet, treatment, tissuedetail using functions
sampleInfo <- sampleInfo %>%
  rowwise()%>%
  mutate(batch = GetSampleBatch(title)) %>%
  mutate(diet = GetSampleDiet(title))%>%
  mutate(treatment = GetTreatmentInfo(title))%>%
  mutate(tissuedetail = GetDetailTissueInfo(title))%>%
  ungroup()


#Construct a sample name based on diet, tissue and treatment
sampleInfo <- sampleInfo %>%
  rowwise()%>%
  mutate(samplename = ConstructSampleName(title))%>%
  ungroup()

#Generate a new column to save replicate information
sampleInfo <- sampleInfo %>%
  mutate(replicate = 0)

#Extract unique samples to count replicates
uniqueSamples <- as.data.frame(table(sampleInfo$samplename))
colnames(uniqueSamples) <- c("samplename","frequency")

#Order the sample info based on sample
sampleInfo <- sampleInfo[order(sampleInfo$sample),]

#For loop to assign replicate numbers
for(sample in uniqueSamples$samplename){
  ind <-  which(sampleInfo$samplename == sample)
  j=1
  for(i in ind){
    
    sampleInfo$replicate[i] <- j
    j <- j + 1
    
  }
}

#Construct sample id by adding batch number and replicate number to the sample name
sampleInfo <- sampleInfo %>%
  rowwise %>%
  mutate(sampleid = ConstructSampleId(samplename,title)) %>%
  ungroup()


#Extract general tissue info
sampleInfo <- sampleInfo %>%
  rowwise()%>%
  mutate(tissuegeneral = ifelse(tissuedetail == "Ileum" ||
                                  tissuedetail == "Duodenum" ||
                                  tissuedetail == "Jejunum" ||
                                  tissuedetail == "Ileum interposition","Small Bowel Mucosa",tissuedetail)) %>%
  ungroup()


#Construct general sample name
sampleInfo <- sampleInfo %>%
  rowwise()%>%
  mutate(samplenamegeneral = ConstructSampleNameGeneral(title)) %>%
  ungroup()


#Generate a new column to save general replicate number information
sampleInfo <- sampleInfo %>%
  mutate(replicate.general = 0)

#Extract unique samples to count replicates
uniqueSamples <- as.data.frame(table(sampleInfo$samplenamegeneral))
colnames(uniqueSamples) <- c("samplenamegeneral","frequency")


#Order the sample info based on sample
sampleInfo <- sampleInfo[order(sampleInfo$sample),]

#For loop to assign replicate numbers
for(sample in uniqueSamples$samplenamegeneral){
  ind <-  which(sampleInfo$samplenamegeneral == sample)
  j=1
  for(i in ind){
    
    sampleInfo$replicate.general[i] <- j
    j <- j + 1
    
  }
}


#Construct general sample id by adding batch number and general replicate number to the general sample name
sampleInfo <- sampleInfo %>%
  rowwise %>%
  mutate(sampleidgeneral = ConstructSampleIdGeneral(samplenamegeneral,title)) %>%
  ungroup()


#Convert sample information to a data frame
sampleInfo <- as.data.frame(sampleInfo)

#Assign row names to sample information
rownames(sampleInfo) <- sampleInfo$sample

#Order sample information based on sample id
sampleInfo <- sampleInfo[order(sampleInfo$sampleid),]

#Write the cleaned sample info to an excel file
write.csv(sampleInfo,file  = "data/SampleInfo.csv",row.names = TRUE)



#Read counts data
ctsInfo <- as.data.frame(data.table::fread("data/annotated_counts_Kachmar_080224.txt",header = TRUE))
rownames(ctsInfo) <- ctsInfo$ensembl_gene_id
samples <- paste0("smpl_",1:160)
cts <- ctsInfo %>% dplyr::select(all_of(samples))
cts <- cts[, !colnames(cts) %in% "smpl_103"]

#Set sampleid as sample name in counts data
samColNames <- colnames(cts)
samid <- c()
for (sample in samColNames) {
  
  
  samid <- c(samid,sampleInfo[sampleInfo$sample==sample,]$sampleid)
  
}

colnames(cts) <- samid


#Write the counts data to an excel file
write.csv(cts,file = "data/counts.csv",row.names = TRUE)


#Read Gene Info
geneInfo <- ctsInfo[,1:4]

#Write the gene information to an excel file
write.csv(geneInfo,file = "data/geneInfo.csv",row.names = TRUE)




