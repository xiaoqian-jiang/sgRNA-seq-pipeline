#### 1. #Batch drawing of violin for seurat object list ####

violin_plot <- function(object, directory, name, featurename){
  violin_before <- list()
  
  library(Seurat)
  
  pdf_filename <- paste0(directory, name)
  pdf(file = pdf_filename, width = 20, height = 10)
  
  for (i in 1:length(object)) {
    violin_before[[i]] <- VlnPlot(object[[i]],
                                  layer = "counts",
                                  features = featurename, 
                                  pt.size = 0.01, 
                                  ncol = 5)     
    print(violin_before[[i]])
  }
  dev.off()
  
}

