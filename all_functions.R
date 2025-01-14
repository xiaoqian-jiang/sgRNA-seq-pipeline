
#Created on  Sep 11 2024
#@author: Xiaoqian Jiang

# The all_functions.R file is a custom script designed to enhance efficiency and reproducibility in single-cell RNA sequencing (scRNA-seq) data analysis. 
# It encapsulates a collection of commonly used, modular functions tailored to streamline repetitive tasks encountered during scRNA-seq workflows. 
# By centralizing these reusable functions, the script minimizes redundancy and ensures consistency across analyses.

# load the library used in the following functions
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
library(cowplot)
library(data.table)

library(msigdbr)

library(enrichplot)
library(GOplot)
library(ComplexHeatmap)


#This script contains all the self-written functions in the scRNA seq analysis pipeline.

#### 1. Get the 36 colors for NGS paper. ####

my36colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FDB462",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7",
  "#BEBADA", "#FB8072", "#80B1D3", "#B3DE69", "#FCCDE5", "#FC8D62",
  "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C",
  "#FB9A99", "#B2DF8A", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#F4E242", "#B3B3B3", "#FFFF33"
)

#### 2. Batch drawing of violin for seurat object list ####

violin_plot <- function(object, directory, name, featurename){
  pdf_filename <- paste0(directory, name)
  pdf(file = pdf_filename, width = 20, height = 10)
  
  violin_before <- VlnPlot(object,
                                  layer = "counts",
                                  features = featurename, 
                                  pt.size = 0.01, 
                                  ncol = 5)     
  print(violin_before)
  dev.off()
}


#### 3. Batch drawing of dimplot for seurat object list ####

dim_plot <- function(object, directory, name, group_by, color){
  
  pdf_filename <- paste0(directory,name)
  pdf(file = pdf_filename, width = 20, height = 10)
  
  if (color){
    colors <- color
  }
  else
    colors <- my36colors

  dim_before <- DimPlot(object, reduction = "umap", 
                               group.by = group_by, 
                               cols = colors )
   
  print(dim_before)
  dev.off()
}


#### 4. Function: seurat_standard_normalize_and_scale  ####
# Normalization, FindVariableFeatures, scale data, FindNeighbors, FindClusters for each sample
# from a initial seurat list-merge file, i.e., before quality control.

seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
  
  # colon: the single seurat object
  # cluster: can be assigned or not
  # cluster_resolution: the resolution_list must be pre-set by author
  
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  
  #The following ScaleData is little time and memory consuming. 
  # The only purpose of scale data is to reduction, which can be removed after run all reduction
  if ("mt_percent" %in% names(colon@meta.data)){
    colon <- ScaleData(colon, vars.to.regress = c("mt_percent"))
  }
  else
    colon <- ScaleData(colon)
  
  
  colon <- RunPCA(colon, features = VariableFeatures(object = colon),verbose=F)
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:50)
  return(colon)
}


#### 5. Add a shrunken coordinate axis to a DimPlot by using shrunk_dimplot(seurat object) ####

#This function is changed based on other's script on website.
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

# Create a left corner umap plot, can change group.by = "celltype" to others
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

#### 6. Find marker genes in each sample between two groups (x,y) in a seurat object ####
library(Seurat)

My_findmarkers <- function(object, logfoldchange, pvalue, pctvalue, only.pos = FALSE,x,y){
  
  degs <- FindMarkers(object, group.by = "group",
                      
                      logfc.threshold = logfoldchange,  # log fold change threshold
                      only.pos = FALSE,        # keep the positive and negative expression genes
                      min.pct = pctvalue,      # Only keep genes detected in a minimum fraction of cells  
                      
                      ident.1 = x, ident.2 = y)
  # Add regulation direction
  degs$regulation <- ifelse(degs$avg_log2FC > 0, "Up", "Down")
  
  # Filter markers based on adjusted p-value
  significant_degs <- subset(degs, p_val_adj < pvalue)
  return(significant_degs)
  
}

#### 7. Get the standard geneset from msigdbr ####


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


#### 8, Comprehensive Analysis and Visualization for Enrichment Results #### 

# Complex processes are encapsulated within a section to simplify the user interface, offering a one-stop solution for comprehensive analysis and visualization of enrichment results. 
# Notice that: this section (the part #8. from file "all_functions.R") includes two key functions: myGOChord and perform_enrichment_analysis, which should be added firstly.  

# Function "myGOChord" is used for circle plot with input of enrichment results, changed from package GOplot::GOChord.

# Function "perform_enrichment_analysis" to run the entire enrichment analysis with 5 input parameters: 
# species: human or mouse 
# category: any category from msigdb, such as H, C1, C2 
# subcategory: subcategory from category, such as BP in C5 
# markers_group: Data (list) for comparison, such MIA_VS_IAC
# top_n : The top significant pathway number to show in the plot
# which are capable of generating four types of visual representations: dot plots, heat maps, circular plots, and bar plots. 

#8.1 myGOChord is used for circle plot with input of enrichment results, changed from package GOplot::GOChord ####

library(GOplot)

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

#8.2 Function to run the entire enriment analysis with 5 input parameters, output enriment results, including heatmap, bar, dot, circle plots ####
perform_enrichment_analysis <- function(species, category, subcategory, markers_group, top_n = 10) {
  # species = "Homo sapiens", category =  "H", subcategory = NULL, keep pathway number = top_n
  
  # get genesets from msigdbr
  gene_sets <- msigdbr(species = species, category = category, subcategory = subcategory)
  # # Modify the names of the pathways in your genesets list by removing the "first_" prefix
  # for (i in seq_along(gene_sets)) {
  #   # Using sub() to remove everything up to and including the first underscore
  #   gene_sets[[i]]@setName <- sub("^[^_]*_", "", gene_sets[[i]]@setName)
  # }
  
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


#### 9. Survival analysis and visualization ####

library(survival)
library(survminer)

# Function "my_survival_analysis" to run the entire survival analysis with 4 input parameters: 
# score: matrix of gsva scores or gene expression for each cluster of each patient
# meta: meta data of patient from clinical database
# type_input : using 'GSVA':gsva_score, or 'LIST':interested gene_list, should be consistent with parameter score
# category: assigning to conduct OS(overall survival) or PFS (Progression-Free Survival) analysis  

my_survival_analysis <- function(score, meta, type_input = 'S', category = 'OS', model = 'KM' ){
  
  splots <- list()
  g = 1

  
  
  # Initialize file name and y-axis title
  file_suffix <- if (category == 'PFS') 'PFS_plot.pdf' else 'OS_plot.pdf'
  Y_title <- if (category == 'PFS') "Progression-Free Survival" else "Overall Survival"
  
  
  # Iterate over each DE gene list
  for (i in rownames(score)) {
    v <- as.numeric(score[i,])
    
    # Classify subgroups based on the type of input: gsva_score use 0 and gene_list use median
    threshold <- if (type_input == 'GENELIST') median(v) else 0
    sub_group <- ifelse(v < threshold, "low", "high")
    
    # Update metadata with subgroup classification
    meta$sub_group <- sub_group
    
    if (model == 'COX'){
      # Fit survival curves, Cox Proportional Hazards models
      #fit_cox <- coxph(surv_object ~ covariate1 + covariate2, data = data)
    
      cox_model  <- coxph(Surv(meta[[time_var]], meta[[event_var]]) ~ sub_group, data = meta)
      fit <- survfit(cox_model, data = meta)
    }
    else{
       # Fit survival curves, Kaplan‒Meier comparative survival analyses
      fit <- survfit(Surv(meta[[time_var]], meta[[event_var]]) ~ sub_group, data = meta)
      
      # The following code is to solve the problem that function "ggsurvplot" cannot use time_var as passing parameter in my function "my_survival_analysis"
      # fit <- eval(parse(text = paste0("survfit(formula = Surv(as.numeric(meta[['", time_var, "']]), 
      #                     as.numeric(meta[['", event_var, "']])) ~ sub_group, data = meta)")))
      
    }
    
    # Create survival plot
    survp <- ggsurvplot(fit, data = meta,
                        surv.median.line = "hv", # Add medians survival
                        pval = TRUE,             # Add p-value and intervals 
                        conf.int = TRUE,        # Add the 95% confidence band
                        risk.table = TRUE,      # Add risk table
                        
                        tables.height = 0.2,
                        tables.theme = theme_cleantable(),
                        palette = "jco",
                        ggtheme = theme_bw(),
                        xlab = "Time(months)",
                        ylab = Y_title,
                        title = i)
    
    print(survp)
    splots[[g]] <- survp
    g <- g + 1
  }

  x_dims = ceiling(sqrt(length(splots)))
  dims = x_dims *5
  
  all_plot <- arrange_ggsurvplots(splots,
                                  print = F,
                                  ncol = x_dims, nrow = x_dims,
                                  risk.table.height = 0.3,
                                  surv.plot.height = 0.7)

  Filename = paste0(type_input,'_',model, '_', file_suffix)
  ggsave(all_plot, filename = Filename, path = "results/", width = dims,height = dims)
  all_plot <- NULL
  
}


