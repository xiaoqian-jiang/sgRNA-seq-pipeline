#This script is conduct the basic quality control of raw scRNA seq data.

#### Step1:Set the clean environment and Load the required libraries ####

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
set.seed(123456)
setwd("C:/Users/xqbus/Desktop/sg_rnaseq/")
getwd()

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr) 



#### Step2:Initial Seurat object #### 

#The Raw data were downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189357. 

project_dir = "GSE189357_RAW/"
files=list.files(project_dir,'^GSM')
#files

#Processing data, organizing the original files into barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz into their respective folders. 

samples=str_split(files,'_',simplify = T)[,1]
if(T){
  lapply(unique(samples),function(x){
    y=files[grepl(x,files)]
    subfolder=paste0(project_dir, paste(str_split(y[1],'_',simplify = T)[,1:2], collapse = "_"))
    dir.create(subfolder,recursive = T)
    file.rename(paste0(project_dir,y[1]),file.path(subfolder,"barcodes.tsv.gz"))
    file.rename(paste0(project_dir,y[2]),file.path(subfolder,"features.tsv.gz"))
    file.rename(paste0(project_dir,y[3]),file.path(subfolder,"matrix.mtx.gz"))
  })
}

dir_name <- list.files(project_dir)
#dir_name


#Initialize Seurat objects with minimal cells 10 and features 500  with raw counts

scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = 
                      paste(project_dir, dir_name[i], sep = ""))
  
  #Batch rename files to names recognized by the Read10X() function.
  name=str_split(dir_name[i],'_',simplify = T)[,2]
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = name,min.cells = 10, 
                                       min.features = 500)
}

#check the initial sample info
scRNAlist

#### Step3:Quality Control #### 

#Calculating mitochondrial, red blood cell, and ribosomal gene ratio 

for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]]#Obtain the i Seurat object in scRNAlist
  
  # Calculating mitochondrial ratio
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
  
  # Calculating red blood cell ratio
  sc[["HB_percent"]] <- PercentageFeatureSet(sc, pattern="^HB[ABDGQEMZ]") 
  
  # Calculate the ribosomal gene ratio
  sc[["RP_percent"]] <- PercentageFeatureSet(sc, pattern = "^Rp[sl]")  
  
  scRNAlist[[i]] <- sc
  rm(sc)
}

#Batch drawing of violin before quality control to check the data quality 
#Using function violin_plot in "all_functions.R"

output_directory <- "results/" 
featurename = c("nFeature_RNA", "nCount_RNA", "mt_percent","HB_percent", "RP_percent")

violin_plot(scRNAlist, output_directory, "violin_before_combined.pdf", featurename)


#Batch filter the cell according to the ratio of MT、HB, RP genes
#the set of nFeature_RNA should be set according to Violin plot
#The general default is that the mitochondrial content should be less than 10%, 
#the number of red blood cells should be less than 3%, 
#and the number of ribosomal should be less than 5%

scRNAlist <- lapply(X = scRNAlist, FUN = function(x){
  x <- subset(x, 
              subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & 
                mt_percent < 10 & 
                HB_percent < 3 & 
                RP_percent < 5 & 
                nCount_RNA < quantile(nCount_RNA,0.97) & nCount_RNA > 1000
  )})

#View(scRNAlist[[1]]@meta.data)
scRNAlist 

#Batch drawing of violin after quality control
violin_plot(scRNAlist, output_directory, "violin_after_combined.pdf", featurename)

#Save the file
save(scRNAlist,file = "./rdata/scRNAlist.Rdata")
#lnames = load("scRNAlist.Rdata")


#### Step4: Sketch A Subset of Cells #### 

#Due to computational limitations, I select a subset (‘sketch’) of 5,000 cells of each sample for downstream analysis.

sketchlist <- list()
for (i in 1:length(scRNAlist)) {
  obj <- scRNAlist[[i]]
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- SketchData(
    object = obj,
    ncells = 5000,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )

  count <- LayerData(obj, assay = "sketch", layer = "count")
  name = paste0('TD',i)
  sketchlist[[i]]  <- CreateSeuratObject(count, project = name)
}

#Batch drawing of violin after sketch to check the data
featurename = c("nFeature_RNA", "nCount_RNA")
violin_plot(sketchlist, output_directory, "violin_after_sketch.pdf", featurename)

#Save the file
save(sketchlist,file = "./rdata/sketchlist.Rdata")
#lnames = load("scRNAlist.Rdata")


