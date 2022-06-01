# Canine spGEX Analysis.R
# Jacqueline Larouche
# Spring 2022

# Contains the code used for spGEX data analysis of canine VML tissues

# Workflow: 
# 1. Initialize variables
# 2. Load ST Data and Create List of Seurat Objects, Perform Seurat Label Transfer
# 3a. Generate seurat object with all ST datasets
# 3b. Plot (UMAPs, Overlays)
# 4. CellChat gene modules
# 5. Marker co-localization

#### 0. Initialize Environment ####
library(Seurat)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(Matrix)
library(data.table)
library(RColorBrewer)
library(tibble)
library(cowplot)
library(reshape)
library(amap)
library(hdf5r)
library(patchwork)
library(tidyverse)
library(ggpubr)

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/TA_VML/SpatialTranscriptomics/Visium/Data/spGEX")
#### 1. Initialize Variables ####
sample_dict <- list()
sample_dict[["Sample_5513-JL-S2-A_CGCGCACT-AGAATACA"]] = "D7_Defect"
sample_dict[["Sample_5513-JL-S2-B_CCTGTCAG-GTTACGGG"]] = "D7_Transition"
sample_dict[["Sample_5513-JL-S2-C_GTCCTTCG-CTGTGCAT"]] = "D14_Defect"
sample_dict[["Sample_5513-JL-S2-D_AATGTATC-TAAGCTCA"]] = "D14_Transition"
samples <- names(sample_dict)

st_list <- readRDS("Data/Canis_ST_list.RDS")
st_all <- readRDS(file = "Data/canis_merged_414T_415T_425T_426T.RDS")

#### 2. Load ST Data and Create List of Seurat Objects ####
st_list <- sapply(samples, USE.NAMES = TRUE, function(geo) {
  print(sample_dict[[geo]])
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/5513-JL/%s", geo)), 
                           assay = "Spatial", filter.matrix = TRUE, 
                           slice = sprintf("%s", sample_dict[[geo]]), 
                           filename = "filtered_feature_bc_matrix.h5")
  # remove spots with <1 gene detected
  st_se <- subset(st_se, subset = nFeature_Spatial > 0)
  st_se$orig.ident <- sample_dict[[geo]]
  st_se$timepoint <- strsplit(sample_dict[[geo]], '_')[[1]][1]
  st_se$location <- strsplit(sample_dict[[geo]], '_')[[1]][2]
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  st_se <- st_se %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(verbose = FALSE) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  return(list(st_se))
})
saveRDS(st_list, file = "Data/Canis_ST_list.RDS")

#### 3a. Generate seurat object with all ST datasets ####
st_1 <- st_list[[samples[1]]]; st_2 <- st_list[[samples[2]]] 
st_3 <- st_list[[samples[3]]]; st_4 <- st_list[[samples[4]]]

st_all <- merge(st_1, c(st_2, st_3, st_4))
table(st_all$timepoint, st_all$location)

st_all <- SCTransform(st_all, assay = "Spatial", verbose = FALSE)
st_all <- st_all %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)

saveRDS(st_all, file = "Data/canis_merged_414T_415T_425T_426T.RDS")

#### 3b. Plot (UMAPs, Overlays, Heatmap) ####
pdf(file="Plots/Finalized/canis_all_umap_dimplot.pdf", onefile = TRUE)
DimPlot(st_all, reduction = "umap", label = FALSE, pt.size = 3) + NoLegend()
dev.off()

pdf(file="Plots/canis_all_umap_qc_filtered.pdf", onefile = TRUE)
FeaturePlot(st_all, reduction = "umap", features = 'nCount_Spatial', pt.size = 1)
FeaturePlot(st_all, reduction = "umap", features = 'nFeature_Spatial', pt.size = 1)
DimHeatmap(st_all, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pdf("Plots/Finalized/Canine_umis_highres.pdf", width = 15, height = 5)
SpatialFeaturePlot(st_all, features = "nCount_Spatial", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/Canine_features_highres.pdf", width = 15, height = 5)
SpatialFeaturePlot(st_all, features = "nFeature_Spatial", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/canine_adgre1_overlay.pdf", onefile = TRUE, width = 15, height = 5)
SpatialFeaturePlot(st_all, features = "ADGRE1", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/canine_aspn_overlay.pdf", onefile = TRUE, width = 15, height = 5)
SpatialFeaturePlot(st_all, features = "ASPN", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/canine_myog_overlay.pdf", onefile = TRUE, width = 15, height = 5)
SpatialFeaturePlot(st_all, features = "MYOG", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/canine_tgfb_overlays.pdf", onefile = TRUE, width = 15, height = 5)
SpatialFeaturePlot(st_all, features = c("TGFB1"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
SpatialFeaturePlot(st_all, features = c("TGFB2"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
SpatialFeaturePlot(st_all, features = c("TGFB3"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
SpatialFeaturePlot(st_all, features = c("TGFBR1"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
SpatialFeaturePlot(st_all, features = c("TGFBR2"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0, max.cutoff = 2)
dev.off()

DefaultAssay(st_all) <- "SCT"
Idents(st_all) <- 'orig.ident'
degs <- FindAllMarkers(st_all, only.pos = TRUE)
top10 <- degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf("Plots/Finalized/canine_deg_heatmap_tissue.pdf", width = 8, height = 6)
DoHeatmap(st_all, features = top10$gene)
dev.off()

#### 4. CellChat Gene Modules ####
DefaultAssay(st_all) <- "SCT"

fap_r <- c('CD44', 'CD47', 'NCL', 'SDC2', 'LRP1', 'SDC1', 'SDC4', 'FGFR1', 'ITGA5', 'ITGB1', 'CDH11', 'ITGA11', 'ITGA4', 'ITGB7', 'AXL', 'EGFR', 'DAG1', 'LRP1', 'NCAM1', 'PDGFRB')
fap_l <- c('COL1A1', 'COL1A2', 'FN1', 'THBS4', 'PTN', 'COL6A1', 'COL6A2', 'THBS2', 'APP', 'MIF', 'MDK', 'HSPG2', 'TNC', 'LAMB1', 'COMP', 'THBS1', 'ANGPTL1', 'THBS3', 'ANGPTL2', 'COL6A3', 'LAMA4', 'LAMC1', 'SPP1', 'TNN', 'ANGPTL4', 'PDGFA', 'NCAM1', 'HBEGF', 'GAS6')

mac_l <- c('THBS1', 'APP', 'MIF', 'CCL6', 'FN1', 'CCL9', 'LGALS9', 'ITGA4', 'ITGB1', 'TNF', 'ITGB7', 'NAMPT', 'CCL3', 'CCL2', 'SELPG', 'SEMA4A')
mac_r <- c('CD74', 'CCR2', 'CCR1', 'CD44', 'SDC4', 'PTPRC', 'IGHM', 'CXCR4', 'SELL', 'PLXNB2', 'CD47', 'TNFRSF1A', 'TNFRSF1B', 'PIRB', 'SDC3', 'LRP1', 'NCL')

musc_l <- c('JAM3', 'COL4A2', 'LAMA2', 'CDH15', 'COL4A1')
musc_r <- c('ITGA7', 'CDH15', 'JAM3', 'VCAM1', 'IGF1R')

st_all <- AddModuleScore(object = st_all, features = list(fap_l), ctrl = 5, 
                        name = 'MP_Ligands')
st_all <- AddModuleScore(object = st_all, features = list(fap_r), ctrl = 5, 
                        name = 'MP_Receptors')

st_all <- AddModuleScore(object = st_all, features = list(mac_l), ctrl = 5, 
                        name = 'Mø_Ligands')
st_all <- AddModuleScore(object = st_all, features = list(mac_r), ctrl = 5, 
                        name = 'Mø_Receptors')

st_all <- AddModuleScore(object = st_all, features = list(musc_l), ctrl = 5, 
                        name = 'MuSC_Ligands')
st_all <- AddModuleScore(object = st_all, features = list(musc_r), ctrl = 5, 
                        name = 'MuSC_Receptors')


pdf("Plots/Finalized/Canine_Overlays_LR_Scores_vertical.pdf", height = 20, width = 12)
SpatialFeaturePlot(st_all, features = c('Mø_Ligands1', 'Mø_Receptors1', 
                                       'MP_Ligands1', 'MP_Receptors1',
                                       'MuSC_Ligands1', 'MuSC_Receptors1'), 
                   min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 2.5, alpha = c(0,1),
                   ncol = 2)
dev.off()

