
# Single-cell RNA-seq Analysis Pipeline

## Overview
This repository provides a comprehensive pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data. The workflow integrates state-of-the-art tools in R for robust data preprocessing, clustering, and downstream analyses. The primary aim is to gain hands-on experience while deriving insights from a case study on the evolution of lung adenocarcinoma (LUAD) using data from Zhu et al. (2022).

## Project Objective
The main objective of this repository is to provide a hands-on exploration of scRNA-seq and spatial transcriptomics (ST) workflows. It covers:
- **Data preprocessing** (Quality control, Multi-sample integration)
- **Cell clustering and annotation**
- **Downstream analyses**: Differential expression, pathway enrichment, survival analysis, cell-cell communication, trajectory inference.

This repository serves as a resource for researchers interested in replicating or customizing workflows for similar datasets.

## Data Source
The datasets used in this project are from Zhu et al. (2022):
> "Delineating the dynamic evolution from preneoplasia to invasive lung adenocarcinoma by integrating single-cell RNA sequencing and spatial transcriptomics."  
DOI: [10.1038/s12276-022-00896-9](https://doi.org/10.1038/s12276-022-00896-9)

### Dataset Overview: scRNA-seq Data
- **Sequencing Platform**: Illumina NovaSeq 6000
- **Number of Cells**: 115,246
- **LUAD Stages**:
  - Adenocarcinoma in situ (AIS): 37,143 cells (32.2%)
  - Minimally invasive adenocarcinoma (MIA): 30,909 cells (26.8%)
  - Invasive adenocarcinoma (IAC): 47,194 cells (41.0%)

## Analysis Workflow
The pipeline is divided into ten modular steps, each focusing on critical components of scRNA-seq analysis:

### 1. Quality Control (`1.QC.R`)
Performs quality control by filtering cells based on:
- Unique gene counts
- Mitochondrial, red blood cell, and ribosomal gene content
- Doublet removal using the **DoubletFinder** algorithm.  
This ensures the dataset consists of high-quality singlet cells.

### 2. Multi-Sample Data Integration (`2.Integration.R`)
Integrates data from 9 scRNA-seq samples and addresses batch effects using the **Harmony** algorithm. This step aligns shared cell populations across samples while preserving biological signals.

### 3. Cell Type Annotation (`3.Celltype_Annotation.R`)
Cell type annotation is performed using:
- Differentially expressed genes
- Known marker genes
- Tools like **CellMarker2**, **SingleR**, and **Azimuth** for automated annotation.

### 4. Differential Gene Expression Analysis (`4.DEGs.R`)
Identifies and visualizes differentially expressed (DE) genes using:
- Volcano plots
- Heatmaps
- Other graphical formats to highlight significant gene expression changes across conditions or clusters.

### 5. Copy Number Variation Analysis (`5.Epithelia_InferCNV.R`)
Focuses on identifying epithelial subtypes, their transformation into cancerous states, and their role in lung cancer development using **inferCNV**.

### 6. Pathway and Functional Enrichment Analysis (`6.Pathway.R`)
Uses GSVA, AUC scores, and ORA to uncover the biological significance behind gene expression patterns.

### 7. Survival Analysis (`7.Survival.R`)
Integrates **Kaplan-Meier estimators** and **Cox Proportional Hazards models** using bulk RNA-seq data from TCGA-LUAD. Survival analysis is conducted using the **survival** R package.

### 8. Trajectory Inference (`8.Trajectory.R`)
Uses **Monocle 2** to infer cell trajectory and visualize the progression of cells through different states or conditions.

### 9. Cell-Cell Communication Analysis (`9.Communication.R`)
Explores cell-cell communication networks using **CellChat** and **NicheNet** for inference, analysis, and visualization.

### 10. Transcription Factor Analysis (`10.Transcription_Factor.R`)
SCENIC (Single-Cell Regulatory Network Inference and Clustering) computational framework was used here to reconstruct gene regulatory networks and identify cell states from epithelial cancer cells.
Note:  Although the script is written in R, the core part of this pipeline runs within PySCENIC, which requires setting up a Conda environment beforehand.

---

## Introduction to `all_functions.R`
The `all_functions.R` script is a custom library designed to enhance the efficiency and reproducibility of scRNA-seq data analysis. It encapsulates reusable, modular functions for common tasks, such as:
- Data preprocessing
- Clustering
- Visualization
- Other downstream analyses.

### How to Use
1. **Source the File**:  
   Include the script in your analysis by adding:
   ```R
   source("path_to_file/all_functions.R")
