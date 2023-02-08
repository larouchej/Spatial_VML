# VML spGEX Analysis.R
# Jacqueline Larouche
# Spring 2022

# Contains the code used for spGEX data analysis of mouse at both 7 and 14dpi time points (Figure 4 and Supp Fig. 3)

# Workflow: 
# 0.  Initialize Environment
# 1.  Initialize variables
# 2.  Create D14 scRNA-Seq reference
# 3a. Load ST Data and Create List of Seurat Objects, Perform Seurat Label Transfer
# 3b. Plot D14 Spatial
# 4a. Generate seurat object with all ST datasets
# 4b. UMAP Plots D7 & D14
# 4c. GOTerm Analysis D14
# 4d. Volcano Plots D7 vs D14

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
library(corrplot)
library(amap)
library(hdf5r)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tidyverse)
library(readxl)
library(ggpubr)
library(EnhancedVolcano)

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/TA_VML/SpatialTranscriptomics/Visium/Data/spGEX")

#### 1. Initialize Variables ####
sample_dict <- list()
sample_dict[["PBS_599L"]] = "3765-JL"
sample_dict[["ITD1_599R"]] = "3765-JL"
sample_dict[["ITD1_600L"]] = "3765-JL"
sample_dict[["PBS_600R"]] = "3765-JL"
sample_dict[["PBS_1197L"]] = "5347-JL"
sample_dict[["PBS_1203L"]] = "5347-JL"
sample_dict[["ITD1_1203R"]] = "5347-JL"
sample_dict[["ITD1_1202L"]] = "5347-JL"
sample_dict[["Sample_6834-JL-S1-A-GEX_GCGGGTAA-CTTAGTGC"]] = "D14_Mouse_1"
sample_dict[["Sample_6834-JL-S1-B-GEX_CCTATCCT-GTTAGTAT"]] = "D14_Mouse_2"
samples <- names(sample_dict)

st_list <- readRDS("Data/Mus_ST_list_label_transfer_1201.RDS")
st_all <- readRDS("Data/merged_D7_D14.RDS")
refd7 <- readRDS('Data/vml_cca_d7.rds')
refd14 <- readRDS('Data/vml_cca_d14.rds')

#### 2. Create D14 scRNA-Seq Reference ####
ref <- readRDS("Data/vml_cca_dec11.rds")
Idents(ref) <- "condition"
refd14 <- subset(ref, ident = c("D14.3mm"))
refd14 <- SCTransform(refd14, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.1)
DimPlot(refd14)
saveRDS(refd14, 'Data/vml_cca_d14.rds')

#### 3a. Load ST Data and Create List of Seurat Objects, Perform Seurat Label Transfer ####
st_list_1 <- sapply(samples[1:8], USE.NAMES = TRUE, function(geo) {
  print(geo)
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/%s/%s", sample_dict[[geo]], geo)), 
                           assay = "Spatial", filter.matrix = TRUE, slice = sprintf("%s", geo), 
                           filename = "filtered_feature_bc_matrix.h5")
  # remove spots with <1 gene detected
  st_se <- subset(st_se, subset = nFeature_Spatial > 0)
  st_se$orig.ident <- geo
  st_se$timepoint <- 'D7'
  st_se$species <- 'Mus'
  st_se$Zone <- read.csv(paste0(getwd(), sprintf("/10xOutputData/%s/%s/Zones.csv", sample_dict[[geo]], geo)), row.names = 1)
  
  comp <- paste0(st_se$timepoint, "_", st_se$Zone)
  m <- data.frame("comp" = comp)
  rownames(m) <- rownames(st_se@meta.data)
  st_se <- AddMetaData(object = st_se, metadata = m)
  
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  st_se <- st_se %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(verbose = FALSE) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  # Seurat Label Transfer
  anchors <- FindTransferAnchors(reference = refd7, query = st_se, 
                                 normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = refd7$Celltype, prediction.assay = TRUE,
                                    weight.reduction = st_se[["pca"]], dims = 1:30,
                                    k.weight = 45)
  st_se[["predictions"]] <- predictions.assay
  return(list(st_se))
})

st_list_2 <- sapply(samples[9:10], USE.NAMES = TRUE, function(geo) {
  print(sample_dict[[geo]])
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/6834-JL/%s", geo)), 
                           assay = "Spatial", filter.matrix = TRUE, 
                           slice = sprintf("%s", sample_dict[[geo]]), 
                           filename = "filtered_feature_bc_matrix.h5")
  
  # Convert coordinates of the spatial image to integers from characters.
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["tissue"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["tissue"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["row"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["row"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["col"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["col"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagerow"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagerow"]])
  st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagecol"]] <- as.integer(st_se@images[[sprintf("%s", sample_dict[[geo]])]]@coordinates[["imagecol"]])
  
  # remove spots with <1 gene detected
  st_se <- subset(st_se, subset = nFeature_Spatial > 0)
  st_se$orig.ident <- sample_dict[[geo]]
  st_se$timepoint <- 'D14'
  st_se$species <- 'Mus'
  st_se$Zone <- read.csv(paste0(getwd(), sprintf("/10xOutputData/6834-JL/%s/Zones.csv", geo)), row.names = 1)
  
  comp <- paste0(st_se$timepoint, "_", st_se$Zone)
  m <- data.frame("comp" = comp)
  rownames(m) <- rownames(st_se@meta.data)
  st_se <- AddMetaData(object = st_se, metadata = m)
  
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  st_se <- st_se %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(verbose = FALSE) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  # Seurat Label Transfer
  anchors <- FindTransferAnchors(reference = refd14, query = st_se, 
                                 normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = refd14$Celltype, prediction.assay = TRUE,
                                    weight.reduction = st_se[["pca"]], dims = 1:30,
                                    k.weight = 15)#was 45
  st_se[["predictions"]] <- predictions.assay
  return(list(st_se))
})

st_list <- c(st_list_1, st_list_2)
saveRDS(st_list, file = "Data/Mus_ST_list_label_transfer_1201.RDS")

#### 3b. D14 Spatial Plots (Fig 4B-D, SF 3A-C) ####
st_se <- st_list[["Sample_6834-JL-S1-A-GEX_GCGGGTAA-CTTAGTGC"]]

DefaultAssay(st_se) <- "SCT"

pdf("Plots/Finalized/D14_clusters_highres.pdf")
SpatialDimPlot(st_se, group.by = "seurat_clusters", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/D14_umis_highres.pdf")
SpatialFeaturePlot(st_se, features = "nCount_Spatial", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/D14_zones_highres.pdf")
SpatialDimPlot(st_se, group.by = "Zone", pt.size.factor = 1, crop = FALSE, alpha = 1, stroke = 1)
dev.off()

pdf("Plots/D14_gene_overlays.pdf", width = 15, height = 10)
SpatialFeaturePlot(st_se, features = c("Myh1", "Aspn", "Cd68", "Pax7",
                                       "Myh4", "Col1a1", "Adgre1", "Myod1"), 
                   pt.size.factor = 1, crop = FALSE, ncol = 4)
dev.off()

pdf("Plots/D14_gene_overlays_SAMs.pdf", width = 15, height = 5)
SpatialFeaturePlot(st_se, features = c("Cd9", 'Trem2', 'Spp1', 'Fabp5'), 
                   pt.size.factor = 1, crop = FALSE, ncol = 4)
dev.off()

DefaultAssay(st_se) <- "predictions"
cell_types_all <- c('Monocyte', 'Dendritic', 'Macrophage', 'Neutrophil',  'B-Cell', 'T/NK-Cell', 'Satellite-Cell', 'Myonuclei', 'Endothelial-Stem', 'Endothelial-Vein', 'Endothelial-Capillary', 'Endothelial-Artery', 'Pericyte',  'Smooth-Muscle', 'Lymph-Vessel', 'FAP-Stem', 'FAP-Pro-Remodeling', 'FAP-Matrix', 'FAP-Adipogenic', 'Tenocyte', 'Neural-Progenitor', 'Neural',  'Erythroblast')


pdf("Plots/Finalized/D14_tenocyte_highres.pdf")
SpatialFeaturePlot(st_se, features = "Tenocyte", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1))
dev.off()

pdf("Plots/Finalized/D14_musc_highres.pdf")
SpatialFeaturePlot(st_se, features = "Satellite-Cell", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1))
dev.off()

pdf("Plots/Finalized/D14_myonuclei_highres.pdf")
SpatialFeaturePlot(st_se, features = "Myonuclei", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1))
dev.off()

pdf("Plots/Finalized/D14_remaining_celltypes_v2.pdf", width = 20, height = 5)
SpatialFeaturePlot(st_se, features = c("Endothelial-Artery", "Smooth-Muscle", "FAP-Stem",
                                       "FAP-Pro-Remodeling", "FAP-Adipogenic", "Myonuclei"), 
                   pt.size.factor = 1, crop = FALSE, ncol = 5)
dev.off()

#### 4a. Generate seurat object with all ST datasets ####
st_1 <- st_list[[samples[1]]]; st_2 <- st_list[[samples[4]]]; 
st_3 <- st_list[[samples[5]]]; st_4 <- st_list[[samples[6]]];
st_5 <- st_list[[samples[9]]]; st_6 <- st_list[[samples[10]]];

st_1$round <- "1"; st_2$round <- "1"
st_3$round <- "2"; st_4$round <- "2";
st_5$round <- "3"; st_6$round <- "3";

st_all <- merge(st_1, c(st_2, st_3, st_4, st_5, st_6))
table(st_all$orig.ident, st_all$timepoint)
table(st_all$orig.ident, st_all$Zone)


Idents(st_all) <- 'Zone'
st_filt <- subset(st_all, ident = c('Defect', 'IntactMuscle', 'Transition'))
table(st_filt$orig.ident, st_filt$Zone)

st_all <- SCTransform(st_filt, assay = "Spatial", verbose = FALSE)
st_all <- st_all %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)

saveRDS(st_all, file = "Data/merged_D7_D14.RDS")

#### 4b. UMAP Plots D7 & D14 ####

pdf(file="Plots/Finalized/mouse_d7_d14_umap_dimplot.pdf", width = 8, height = 3)
p1 <- dittoDimPlot(st_all, 'orig.ident', color.panel = c(hue_pal()(6)))
p2 <- dittoDimPlot(st_all, 'Zone', color.panel = c(hue_pal()(3)))
p1 + p2 + plot_layout(ncol = 2)
dev.off()

pdf(file="Plots/mouse_d7_d14_umap_qc_filtered.pdf", onefile = TRUE)
FeaturePlot(st_all, reduction = "umap", features = 'nCount_Spatial', pt.size = 1)
FeaturePlot(st_all, reduction = "umap", features = 'nFeature_Spatial', pt.size = 1)
dev.off()

DefaultAssay(st_all) <- 'SCT'
Idents(st_all) <- 'Zone'
genes <- !grepl(pattern = "^Rp[l|s]|mt|Gm", x = rownames(st_all))
genes <- rownames(st_all)[genes]
my_levels <- c('Defect', 'Transition', 'Intact')
Idents(st_all) <- factor(Idents(st_all), levels= my_levels)

zone.degs <- FindAllMarkers(st_all, only.pos = TRUE, features = genes)
top10 <- zone.degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

comp.degs <- FindAllMarkers(st_all, only.pos = TRUE, features = genes)
top10.comp <- comp.degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

st_def_trans <- subset(st_all, idents = c('Intact'), invert = TRUE)
Idents(st_def_trans) <- 'orig.ident'
st_def_trans <- subset(st_def_trans, idents = c('D14_Mouse_2'), invert = TRUE)
Idents(st_def_trans) <- 'comp'
my_levels <- c('D7_Defect', 'D14_Defect',
               'D7_Transition', 'D14_Transition')
Idents(st_def_trans) <- factor(Idents(st_def_trans), levels= my_levels)


pdf("Plots/Finalized/D14_Zone_Markergenes_Heatmap.pdf", width = 6, height = 8)
DoHeatmap(st_def_trans, features = c('S100a8', 'S100a9', 'Ctss', 'Ccl7', 'Il1b',
                                     top10$gene[1:10],  
                                     top10$gene[11:20],
                                     top10.comp$gene[31:40]), 
          slot = 'scale.data') +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(text = element_text(color = 'black'))
dev.off()


#### 4c. GOTerm Analysis D14 ####

# PATHWAY ANALYSIS ON ZONES
DefaultAssay(st_def_trans) <- 'SCT'
all.genes <- rownames(st_def_trans)
st_def_trans <- ScaleData(st_def_trans, features = all.genes)
# remove mitochondrial and ribosomal genes
genes <- !grepl(pattern = "^Rp[l|s]|mt", x = rownames(st_def_trans))
genes <- rownames(st_def_trans)[genes]
Idents(st_def_trans) <- 'timepoint'
de_markers <- FindAllMarkers(st_def_trans, only.pos = TRUE, assay = 'SCT', features = genes)
de_markers$entrez <- mapIds(org.Mm.eg.db, rownames(de_markers),'ENTREZID', 'SYMBOL')
top10 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## GOTerm analysis
go_concat <- NULL
for (i in c('D7', 'D14')){
  go_temp <- enrichGO(de_markers$entrez[de_markers$cluster == i],
                      OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, qvalueCutoff = 1)
  go_concat <- c(go_concat, go_temp)
}

pdf(file="Plots/Finalized/Timepoint_enrichGO_D14_Defect_Transition.pdf", onefile = TRUE, width = 6, height = 3)
#dotplot(go_concat[[2]], showCategory=15)
dotplot(go_concat[[2]], showCategory=c('muscle structure development',
                                       'generation of precursor metabolites and energy',
                                       'muscle cell differentiation',
                                       'aerobic respiration',
                                       'oxidative phosphorylation'))
dev.off()
pdf(file="Plots/Finalized/Timepoint_enrichGO_D7_Defect_Transition.pdf", onefile = TRUE, width = 6, height = 3)
#dotplot(go_concat[[1]], showCategory=20)
dotplot(go_concat[[1]], showCategory=c('response to stress',
                                       'immune response',
                                       'cellular response to chemical stimulus',
                                       'cell migration',
                                       'immune system process',
                                       'cell death'))

dev.off()

DefaultAssay(st_trans) <- 'SCT'
all.genes <- rownames(st_trans)
st_trans <- ScaleData(st_trans, features = all.genes)
# remove mitochondrial and ribosomal genes
genes <- !grepl(pattern = "^Rp[l|s]|mt", x = rownames(st_trans))
genes <- rownames(st_trans)[genes]
Idents(st_trans) <- 'timepoint'
de_markers <- FindAllMarkers(st_trans, only.pos = TRUE, assay = 'SCT', features = genes)
de_markers$entrez <- mapIds(org.Mm.eg.db, rownames(de_markers),'ENTREZID', 'SYMBOL')
top10 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

## GOTerm analysis
go_concat <- NULL
for (i in c('D7', 'D14')){
  go_temp <- enrichGO(de_markers$entrez[de_markers$cluster == i],
                      OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, qvalueCutoff = 1)
  go_concat <- c(go_concat, go_temp)
}

pdf(file="Plots/Finalized/Timepoint_Transition_enrichGO_D14.pdf", onefile = TRUE, width = 6, height = 3)
dotplot(go_concat[[2]], showCategory=15)
dotplot(go_concat[[2]], showCategory=c('muscle structure development',
                                       'generation of precursor metabolites and energy',
                                       'muscle cell differentiation',
                                       'aerobic respiration',
                                       'oxidative phosphorylation'))
dev.off()
pdf(file="Plots/Finalized/Timepoint_Transition_enrichGO_D7.pdf", onefile = TRUE, width = 6, height = 3)
dotplot(go_concat[[1]], showCategory=15)
dotplot(go_concat[[1]], showCategory=c('response to stress',
                                       'cellular response to chemical stimulus',
                                       'cell migration',
                                       'immune system process',
                                       'cell death'))

dev.off()

#### 4d. Volcano Plots D7 vs D14 ####
Idents(st_all) <- 'Zone'
st_defect <- subset(st_all, idents = 'Defect')
Idents(st_defect) <- 'timepoint'
deg.defect <- FindMarkers(st_defect,ident.1 = 'D14',logfc.threshold = 0)
deg.defect$gene <- rownames(deg.defect)

goi_defect <- c('S100a8', 'S100a9', 'Cd68', 'Ctss', 'Myh1', 'Ccl8', 'Spp1', 'Ttn', 'Tpm1', 'Acta1', 'Trem2', 'Myh3', 'Acta2', 'Smad3')

vol_d <- EnhancedVolcano(deg.defect,
                         lab = rownames(deg.defect),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = '',
                         subtitle = '',
                         xlab = bquote(~Log[2]~ ("fold change")),
                         ylab = bquote(~-Log[10]("p-adj")),
                         axisLabSize = 25,
                         caption = NULL,
                         selectLab = goi_defect,
                         pCutoff = 0.05,
                         FCcutoff = 0.0585,
                         pointSize = 1.0,
                         labSize = 8.0,
                         boxedLabels = FALSE,
                         colAlpha = 0.5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         widthConnectors = 0.2,
                         colConnectors = 'black',
                         col = c('#999999', '#009E73', '#56B4E9', '#E69F00'))
vol_d

pdf("Plots/Finalized/volcano_defect_d14_v_d7.pdf", width = 7, height = 7)
vol_d
dev.off()

st_trans <- subset(st_all, idents = 'Transition')
Idents(st_trans) <- 'timepoint'
deg.trans <- FindMarkers(st_trans,ident.1 = 'D14',logfc.threshold = 0)
deg.trans$gene <- rownames(deg.trans)

goi_trans <- c('Ccl8', 'Myh3', 'Spp1', 'Trem2', 'Thbs2', 'Ms4a7', 'S100a8', 'C1qa', 'C1qb', 'Myh4', 'Krt18', 'Myh1', 'Tnnt3', 'Ltbp4')

vol_t <- EnhancedVolcano(deg.trans,
                         lab = rownames(deg.trans),
                         x = 'avg_log2FC',
                         y = 'p_val_adj',
                         title = '',
                         subtitle = '',
                         xlab = bquote(~Log[2]~ ("fold change")),
                         ylab = bquote(~-Log[10]("p-adj")),
                         axisLabSize = 25,
                         caption = NULL,
                         #labvjust = -4,
                         selectLab = goi_trans,
                         #ylim = c(0,5),
                         #xlim = c(-.5, .5),
                         pCutoff = 0.05,
                         FCcutoff = 0.0585,
                         pointSize = 1.0,
                         labSize = 8.0,
                         boxedLabels = FALSE,
                         colAlpha = 0.5,
                         legendPosition = 'none',
                         drawConnectors = TRUE,
                         widthConnectors = 0.2,
                         colConnectors = 'black',
                         col = c('#999999', '#009E73', '#56B4E9', '#E69F00'))
vol_t

pdf("Plots/Finalized/volcano_trans_d14_v_d7.pdf", width = 7, height = 7)
vol_t
dev.off()
