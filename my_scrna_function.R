

#2. The basic pacakges needed ####

library(devtools)

library(future)
library(patchwork)

library(ggplot2)
library(gplots)
library(ggsci)
library(ggpubr)
library(tidyr)
library(dplyr) 
library(tidyverse)
library(tinyarray)
library(pheatmap)
library(reshape2)
library(ggthemes)
library(RColorBrewer)
library(gridExtra)

library(Seurat)
library(MAST)
library(clustree)
library(cowplot)
library(data.table)

library(msigdbr)
library(GSEABase)

library(clusterProfiler)
library(enrichplot)
library(GOplot)
library(ComplexHeatmap)


#2. Set the clean environment ####

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls()) # remove all environment variables 
set.seed(123456)
plan("multicore", workers = 2) ### 6 compute cores
options(future.globals.maxSize = 8 * 1024^3) #4*4=16G maximum memory allocated
setwd("C:/Users/xqbus/Desktop/sgRNASeq_5/")
getwd()

#3. Get the 36 colors for NGS paper. ####

my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175'
)



#4. Quality Control of seurat for each sample ####

# Normalization, FindVariableFeatures, scale data, FindNeighbors, FindClusters for EACH SAMPLE
# from a initial seurat list-merge file, i.e., before quality control.
resolution_list = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)

seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
  # colon is seurat object, 
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  # all.genes <- rownames(colon)
  # colon <- ScaleData(colon, features = all.genes)
  
  #The following ScaleData is little time consuming, and the only purpose of scale data is to reduction,
  #and memory consuming, can be removed after run all reduction
  colon <- ScaleData(colon, vars.to.regress = c("mt_percent"))
  #scRNAmerge@array@RNA@layers[["scale.data"]] < NULL
  
  colon <- RunPCA(colon, features = VariableFeatures(object = colon),verbose=F)
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:50)
  return(colon)
}

# Quality control: function for basic seurat based qc and doubletfinder based doublet removal
make_seurat_object_and_doublet_removal <- function(data_directory, project_name, DoubletRate ){
  # DoubletRate: Assuming 5% doublet formation rate - tailor for your dataset, 8% for 20,000 cells in 10X data.
  colon.data <- Read10X(data.dir = data_directory)
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 5, min.features = 300)
  
  # Calculating mitochondrial ratio
  currentSample[["mt_percent"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
  # Calculating red blood cell ratio
  currentSample[["HB_percent"]] <- PercentageFeatureSet(currentSample, pattern="^HB[ABDGQEMZ]")   
  # Calculate the ribosomal gene ratio
  currentSample[["RP_percent"]] <- PercentageFeatureSet(currentSample, pattern = "^Rp[sl]")  
  
  #qc plot-pre filtering
  #!!!!!!  change the pathway according to your own environment
  setwd("C:/Users/xqbus/Desktop/sgRNASeq_5/qc.plot/")
  pdf(paste0("./qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "mt_percent"), ncol = 5, pt.size = 0.001))
  dev.off()
  
  # filter everything to 500 unique genes/cell
  currentSample <- subset(currentSample, subset =  nFeature_RNA > 500 & nCount_RNA > 1000 & 
                            mt_percent < 10  &  HB_percent < 3 & RP_percent < 5 )
  
  # Normalize and make UMAP
  currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
  
  # Run doublet finder:Assuming 8% doublet formation rate - tailor for your dataset
  pc.num <- 1:30
  DoubletRate = DoubletRate 
  
  # Find pK
  sweep.res <- paramSweep(currentSample, PCs = pc.num, sct = F) # we didn't run sctransform 
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  
  print(paste(pK_bcmvn, "pK_bcmvn!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", sep=" "))
  
  # # Calculate homotypic doublets probability and the expected doublets number
  # homotypic.prop <- modelHomotypic(currentSample$seurat_clusters)   
  # nExp_poi <- round(DoubletRate * ncol(currentSample))
  # nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  nExp_poi.adj <- round(0.08*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  
   
  seu_colon <- doubletFinder(currentSample, PCs = pc.num, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  #print(head(seu_colon@meta.data))
  
  # rename columns
  name = paste0("DF.classifications_0.25_", pK_bcmvn,"_")
  seu_colon$doublet.class <- seu_colon[[paste0(name,nExp_poi.adj )]]
  seu_colon[[paste0(name,nExp_poi.adj )]] <- NULL
  pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
  seu_colon$pANN <- seu_colon[[pann]]
  seu_colon[[pann]] <- NULL
  print(project_name)
  print(table(seu_colon$doublet.class))
  # plot pre and post doublet finder results
  pdf(paste0("./UMAP_pre_double_removal_", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
  dev.off()
  seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
  pdf(paste0("./UMAP_post_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", cols = c("#D51F26")))
  dev.off()
  
  # Remove extra stuff and return filtered Seurat object
  seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
  return(seu_colon)
}





#5. Add a shrunken coordinate axis to a DimPlot by using shrunk_dimplot(seuratobject) ####

##Create a function to generate axes data
create_shrunk_axes <- function(object) {
  # Extract UMAP coordinates
  umap_coords <- Embeddings(object = object, reduction = 'umap') %>%
    data.frame()
  
  # Calculate axis positions
  lower_bound <- floor(min(min(umap_coords$umap_1), min(umap_coords$umap_2))) - 2
  axis_length <- abs(0.3 * lower_bound) + lower_bound
  midpoint <- abs(0.3 * lower_bound) / 2 + lower_bound
  
  # Create axes lines data
  axes_data <- data.frame(
    x = c(lower_bound, lower_bound, lower_bound, axis_length),
    y = c(lower_bound, axis_length, lower_bound, lower_bound),
    group = c(1, 1, 2, 2),
    label = rep(c('UMAP_2', 'UMAP_1'), each = 2)
  )
  
  # Create axis labels data
  label_data <- data.frame(
    lab = c('UMAP_2', 'UMAP_1'),
    angle = c(90, 0),
    x = c(lower_bound - 3, midpoint),
    y = c(midpoint, lower_bound - 2.5)
  )
  
  # Return the data as a list
  return(list(axes = axes_data, label = label_data))
}

#create a left corner umap plot, can change group.by = "celltype" to others
shrunk_dimplot <- function (object) {
  result <- create_shrunk_axes(object)
  axes <- result$axes
  label <- result$label
  umap =DimPlot(object, reduction = "umap",cols = my36colors,pt.size = 0.8,
                group.by = "celltype",label = TRUE,label.box = TRUE) +
    NoAxes() + 
    theme(aspect.ratio = 1) +
    geom_line(data = axes,
              aes(x = x,y = y,group = group),
              arrow = arrow(length = unit(0.1, "inches"),
                            ends="last", type="closed")) +
    geom_text(data = label,
              aes(x = x,y = y,angle = angle,label = lab),fontface = 'italic')+
    theme(plot.title = element_blank())
  return(umap)
  
}

#5.Check result in different solutions and select proper solution combined with clustree figure ####

solutions_plot <- function(scRNAall, resolution_list, subname) {
  combined_plot <- NULL
  for (i in 1:length(resolution_list)) {
    name <- paste0("RNA_snn_res.", resolution_list[i])
    print(name) 
    p <- DimPlot(scRNAall, reduction = "umap", group.by = name,cols = my36colors)
    
    # Combine the plots in the loop using patchwork's '+' operator
    if (is.null(combined_plot)) {
      combined_plot <- p  
    } else {
      combined_plot <- combined_plot + p  
    }
  }
  combined_plot <- combined_plot + plot_layout(ncol = 3)
  path = paste0("results/",subname)
  filename <- paste0(path,"_resolution_diff.pdf")
  pdf(file = filename, width = 12, height = 10)
  print(combined_plot)
  dev.off()

}

#6. Get the standard geneset from msigdbr ####


# gmt files can be downloaded directly: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
get_gene_sets <- function(species, category) {
  # Retrieve the gene sets based on species and category
  all_gene_sets <- msigdbr(species = species, category = category)
  
  # Group gene symbols by their gene set names and remove duplicates
  gs <- split(all_gene_sets$gene_symbol, all_gene_sets$gs_name)
  gs <- lapply(gs, unique)
  
  # Convert the grouped gene sets to a GeneSetCollection object
  genesets <- GeneSetCollection(mapply(function(geneIds, setId) {
    GeneSet(geneIds,
            geneIdType = EntrezIdentifier(), # Create a GeneSet object using EntrezIdentifier as the gene ID type
            collectionType = KEGGCollection(setId), # Use KEGGCollection as the genoset type
            setName = setId) # Sets the gene set name to keggId
  }, gs, names(gs), SIMPLIFY = FALSE)) # gs is the list of gene ids, and names(gs) is the name of the gene set
  
  return(genesets)
}



#7.1 myGOChord is used for circle plot with input of enrichment results, changed from package GOplot::GOChord ####

myGOChord <- function (data, title, space, gene.order, gene.size, gene.space, 
                       nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, 
                       process.label, limit) 
{
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  if (missing(title)) 
    title <- ""
  if (missing(space)) 
    space = 0
  if (missing(gene.order)) 
    gene.order <- "none"
  if (missing(gene.size)) 
    gene.size <- 3
  if (missing(gene.space)) 
    gene.space <- 0.2
  if (missing(lfc.col)) 
    lfc.col <- c("brown1", "azure", "cornflowerblue")
  if (missing(lfc.min)) 
    lfc.min <- -3
  if (missing(lfc.max)) 
    lfc.max <- 3
  if (missing(border.size)) 
    border.size <- 0.5
  if (missing(process.label)) 
    process.label <- 11
  if (missing(limit)) 
    limit <- c(0, 0)
  if (gene.order == "logFC") 
    data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == "alphabetical") 
    data <- data[order(rownames(data)), ]
  if (sum(!is.na(match(colnames(data), "logFC"))) > 0) {
    if (nlfc == 1) {
      cdata <- GOplot:::check_chord(data[, 1:(Ncol - 1)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[match(x, 
                                                            rownames(data)), Ncol])
    }
    else {
      cdata <- GOplot:::check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[, 
                                                      (Ncol - nlfc + 1)])
    }
  }
  else {
    cdata <- GOplot:::check_chord(data, limit)
    lfc <- 0
  }
  if (missing(ribbon.col)) 
    colRib <- grDevices::rainbow(dim(cdata)[2])
  else colRib <- ribbon.col
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 
                                                      202 * nrib[b]))
  r1 <- 1
  r2 <- r1 + 0.1
  xmax <- c()
  x <- 0
  for (r in 1:length(nrib)) {
    perc <- nrib[r]/sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1) 
      x <- c(x, x[r] + pi * perc)
  }
  xp <- c()
  yp <- c()
  l <- 50
  for (s in 1:Ncol) {
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + 
                                                         xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), 
            r2 * sin(x[s]))
    yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + 
                                                         xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), 
            r2 * cos(x[s]))
  }
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), 
                                                    each = 4 + 2 * l))
  xp <- c()
  yp <- c()
  logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi/Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi/Nrow + space, length = Nrow)
  for (s in 1:Nrow) {
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1) {
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * 
                sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), 
              r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), 
              r2 * sin(x2[s]))
      yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * 
                cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), 
              r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), 
              r2 * cos(x2[s]))
    }
    else {
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc) {
        logs <- c(logs, data[s, (dim(data)[2] + 1 - 
                                   t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * 
                  sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), 
                tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 
                                                          1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * 
                  cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), 
                tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 
                                                          1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
      }
    }
  }
  if (lfc[1] != 0) {
    if (nlfc == 1) {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                      each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 
                                                                                       2 * l))
    }
    else {
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc * 
                                                             Nrow)), each = 4 + 2 * l), logFC = rep(logs, 
                                                                                                    each = 4 + 2 * l))
    }
  }
  else {
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), 
                                                    each = 4 + 2 * l))
  }
  aseq <- seq(0, 180, length = length(x2))
  angle <- c()
  for (o in aseq) if ((o + 270) <= 360) 
    angle <- c(angle, o + 270)
  else angle <- c(angle, o - 90)
  df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + 
                                                         xmax2/2), ygen = (r1 + gene.space) * cos(x2 + xmax2/2), 
                        labels = rownames(cdata), angle = angle)
  df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax/2), 
                        ypro = (r1 + 0.15) * cos(x + xmax/2), labels = colnames(cdata), 
                        stringsAsFactors = FALSE)
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c()
  y.end <- c()
  processID <- c()
  for (gs in 1:length(x2)) {
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 
                 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)) {
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID, -df_bezier$y.end), 
  ]
  x.start <- c()
  y.start <- c()
  for (rs in 1:length(x)) {
    val <- seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 
                 1)
    for (v in 1:(length(val) - 1)) {
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- GOplot:::bezier(df_bezier, colRib)
  if (length(df_genes$logFC) != 0) {
    tmp <- sapply(df_genes$logFC, function(x) ifelse(x > 
                                                       lfc.max, lfc.max, x))
    logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, 
                                            lfc.min, x))
    df_genes$logFC <- logFC
  }
  g <- ggplot() + geom_polygon(data = df_process, aes(x, y, 
                                                      group = id), fill = "gray70", inherit.aes = F, color = "black") + 
    geom_polygon(data = df_process, aes(x, y, group = id), 
                 fill = cols, inherit.aes = F, alpha = 0.6, color = "black") + 
    geom_point(aes(x = xpro, y = ypro, size = factor(labels, 
                                                     levels = labels), shape = NA), data = df_texp) + 
    guides(size = guide_legend("GO Terms", nrow = ncol(data) - 1, byrow = T, 
                               override.aes = list(shape = 22, fill = unique(cols), 
                                                   size = 5))) + theme(legend.text = element_text(size = process.label)) + 
    geom_text(aes(xgen, ygen, label = labels, angle = angle), 
              data = df_texg, size = gene.size) + geom_polygon(aes(x = lx, 
                                                                   y = ly, group = ID), data = df_path, fill = colRibb, 
                                                               color = "black", size = border.size, inherit.aes = F) + 
    labs(title = title) + GOplot:::theme_blank
  if (nlfc >= 1) {
    g + geom_polygon(data = df_genes, aes(x, y, group = id, 
                                          fill = logFC), inherit.aes = F, color = "black") + 
      scale_fill_gradient2("logFC", space = "Lab", low = lfc.col[3], 
                           mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", 
                                                                                       title.hjust = 0.5), breaks = c(min(df_genes$logFC), 
                                                                                                                      max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)), 
                                                                                                                                                       round(max(df_genes$logFC)))) + 
      theme(legend.position = "right", 
            text = element_text(size = 8, face = 'bold'),
            legend.background = element_rect(fill = "transparent"), 
            legend.box = "vertical", legend.direction = "vertical")
  }
  else {
    g + geom_polygon(data = df_genes, aes(x, y, group = id), 
                     fill = "gray50", inherit.aes = F, color = "black") + 
      theme(legend.position = "right", legend.background = element_rect(fill = "transparent"), 
            text = element_text(size = 8, face = 'bold'),
            legend.box = "vertical", legend.direction = "vertical")
  }
}



#7.2 Function to run the entire enriment analysis with 5 input parameters, output enriment results, including heatmap, bar, dot, circle plots ####
perform_enrichment_analysis <- function(species, category, subcategory, markers_group, top_n = 10) {
  # species = "Homo sapiens", category =  "H", subcategory = NULL, keep pathway number = top_n
  
  # get genesets from msigdbr
  gene_sets <- msigdbr(species = species, category = category, subcategory = subcategory)
  
  # loop for each group
  plot_list <- lapply(names(markers_group), function(comp) {
    # Obatain DEGs gene name list
    degs.list <- rownames(markers_group[[comp]])
    # Conduct enrichment analysis, can adjust the pvalueCutoff to get significant pathways
    res <- enricher(degs.list, TERM2GENE = gene_sets[, c(3, 4)], pvalueCutoff = 0.05)
    
    #1 Preparing data for GOplot :: myGOChord (my own function)
    res_cir <- head(res, top_n)
    # loop each pathway to get the top significant genes in each pathway
    significant_genes_by_pathway <- apply(res_cir, 1, function(pathway) {
      #split geneID to get gene name in pathway
      genes <- unlist(strsplit(as.character(pathway[11]), "/"))
      # Get the avg_log2FC of genes in pathway from markers_group, and keep the most significant 10 genes
      gene_expression <- markers_group[[comp]] %>%
        rownames_to_column('gene') %>%
        dplyr::filter(degs.list %in% genes) %>%
        arrange(desc(abs(avg_log2FC))) %>%
        dplyr::slice(1:10) %>%
        .$gene
    })
    g<- unique(unlist(c(significant_genes_by_pathway)))
    
    genedata1 = markers_group[[comp]][g, ] %>% arrange(desc(avg_log2FC))
    genedata2=data.frame(ID=rownames(genedata1),logFC=genedata1$avg_log2FC)
    circ <- circle_dat(res_cir %>%
                         dplyr::rename(adj_pval = p.adjust,
                                       Genes = geneID) %>%
                         dplyr::mutate(Category = category,
                                       Term = ID) %>%
                         dplyr::select(c("ID", "Term", "Genes", "adj_pval", "Category")) %>%
                         dplyr::mutate(Genes = stringr::str_replace_all(Genes, "/", ",")),
                       genedata2)
    chord<-chord_dat(circ, genedata2)
    
    # Circle plot using myGOChord 
    cir_plot <- myGOChord(chord, gene.order = 'log2FC', process.label = 5)
    
    #2 Heatmap plot using Heatmap {ComplexHeatmap}, actually show the gene log2FC is not much useful.
    df <- chord
    transformed_df <- apply(chord[, -ncol(chord)], 2, function(x) ifelse(x == 1, chord[,"logFC"], 0))
    transformed_df <- t(transformed_df)
    transformed_df <- round(transformed_df, 2)
    
    draw_heatmap <- function(data) {
      heatmap <- Heatmap(data,
                         name = "Expression Level",
                         column_title = "Gene",
                         row_title = "Pathway",
                         column_names_rot = 45,  # Rotate column names
                         column_names_gp = gpar(fontsize = 12),
                         row_names_gp = gpar(fontsize = 12),
                         heatmap_legend_param = list(
                           title = "logFC",
                           title_position = "topcenter",  # Title at the top of the legend
                           legend_direction = "vertical",  # Vertical legend
                           legend_width = unit(2, "cm"),
                           legend_height = unit(5, "cm"),
                           labels_fun = function(x) sprintf("%.2f", x),  # Format legend labels to 2 decimal places
                           at = seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = 5)
                         ),
                         color = colorRampPalette(c("blue", "white", "red"))(100)
      )
      draw(heatmap, heatmap_legend_side = "right")  # Draw the heatmap with the legend on the right side
    }
    heat_plot <- draw_heatmap(transformed_df)
    
    #3 plot dotplot using ggplot2 
    # The output data structure of enricher is similar to enrichGO, ggplot can re-plot the figure.
    k = data.frame(res)
    before <- as.numeric(sub("/\\d+$", "", k$GeneRatio))
    after <- as.numeric(sub("^\\d+/", "", k$GeneRatio))
    k$GeneRatio = before /after
    font.size =10
    dot_plot <- k %>% 
      # order according to p.ajust
      arrange(p.adjust) %>% 
      dplyr::slice(1:top_n) %>% 
      ggplot(aes(GeneRatio,forcats::fct_reorder(Description,Count)))+ 
      # plot points
      geom_point(aes(color=p.adjust, size = Count)) +
      # adjust color，guide_colorbar to adjust the figure direction
      scale_color_continuous(low="red", high="blue", guide=guide_colorbar(reverse=TRUE))+
      # buble size
      scale_size_continuous(range=c(3, 8))+
      # if use ylab(""), or left shows blank 
      labs(y=NULL) +
      # if no this setting, the top will reach ceiling of the figure 
      ggtitle("")+
      
      theme_bw() +
      theme(axis.text.x = element_text(colour = "black",
                                       size = font.size, vjust =1 ),
            axis.text.y = element_text(colour = "black",
                                       size = font.size, hjust =1 ),
            axis.title = element_text(margin=margin(10, 5, 0, 0),
                                      color = "black",size = font.size),
            axis.title.y = element_text(angle=90))
    
    # # Or can use default dotplot directly
    # dot_plot1 <- dotplot(res, showCategory = top_n) + ggtitle(paste("Enrichment Analysis for", comp))
    
    #4 plot dotplot using clusterprofile function barplot 
    bar_plot <- barplot(res, showCategory = top_n) + ggtitle(paste("Enrichment Analysis for", comp))
    
    return(list(result = res, 
                plot_dot = dot_plot, 
                plot_heat = heat_plot,
                plot_cir = cir_plot,
                plot_bar = bar_plot))
    
  })
  
  # Return all the results: data and plot
  return(plot_list)
}

