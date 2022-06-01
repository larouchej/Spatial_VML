# VML spGEX Analysis.R
# Jacqueline Larouche
# Spring 2022

# Contains the code used for spGEX data analysis of mouse

# Workflow: 
# 0. Initialize Environment
# 1. Initialize variables
# 2. Create scRNA-Seq reference
# 3. Load ST Data and Create List of Seurat Objects, Perform Seurat Label Transfer
# 3b. Plot PBS
# 3c. Plot ITD1
# 4a. Generate seurat object with all ST datasets
# 4b. Plot (UMAPs)
# 5a. PBS treated tissue plots/analysis
# 5b. ITD1 treated tissue plots/analysis
# 6. Pathway differences by treatment
# 7. CellChat
# 8. CellChat gene modules
# 9. Volcano plots of zone DEGs

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
library(SeuratDisk)
library(ggpubr)

setwd("~/Documents/UMichigan/Aguilar_Lab/ResearchProjects/TA_VML/SpatialTranscriptomics/Visium/Data/spGEX")

#### 1. Initialize Variables ####
cell_types_all <- c('Monocyte', 'Dendritic', 'Macrophage', 'Neutrophil',  'B-Cell', 'T/NK-Cell', 'Satellite-Cell', 'Myonuclei', 'Endothelial-Stem', 'Endothelial-Vein', 'Endothelial-Capillary', 'Endothelial-Artery', 'Pericyte',  'Smooth-Muscle', 'Lymph-Vessel', 'FAP-Stem', 'FAP-Pro-Remodeling', 'FAP-Matrix', 'FAP-Adipogenic', 'Tenocyte', 'Neural-Progenitor', 'Neural',  'Erythroblast')


sample_dict <- list()
sample_dict[["PBS_599L"]] = "3765-JL"
sample_dict[["ITD1_599R"]] = "3765-JL"
sample_dict[["ITD1_600L"]] = "3765-JL"
sample_dict[["PBS_600R"]] = "3765-JL"
sample_dict[["PBS_1197L"]] = "5347-JL"
sample_dict[["PBS_1203L"]] = "5347-JL"
sample_dict[["ITD1_1203R"]] = "5347-JL"
sample_dict[["ITD1_1202L"]] = "5347-JL"
samples <- names(sample_dict)

st_list <- readRDS("Data/Mus_ST_list_label_transfer.RDS")
st_all <- readRDS("Data/merged_filtered_599L_599R_600L_600R_1197L_1203L_1203R_1202L.RDS")
ref_d7 <- readRDS(file = "Data/vml_cca_d7.rds")

#### 2. Create scRNA-Seq Reference ####
ref <- readRDS("Data/vml_cca_dec11.rds")
Idents(ref) <- "condition"
ref_d7 <- subset(ref, ident = c("D7.3mm"))
ref_d7 <- SCTransform(ref_d7, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

cell_types_all <- c('Monocyte', 'Dendritic', 'Macrophage', 'Neutrophil',  'B_Cell', 'T/NK_Cell', 'Satellite_Cell', 'Myonuclei', 'Endothelial_Stem', 'Endothelial_Vein', 'Endothelial_Capillary', 'Endothelial_Artery', 'Pericyte',  'Smooth_Muscle', 'Lymph_Vessel', 'FAP_Stem', 'FAP_Pro-Remodeling', 'FAP_Matrix', 'FAP_Adipogenic', 'Tenocyte', 'Neural_Progenitor', 'Neural',  'Erythroblast')
cell_types_reduced <- c('Monocyte', 'Dendritic', 'Macrophage', 'Neutrophil',  'Lymphocyte', 'Lymphocyte', 'Satellite-Cell', 'Myonuclei', 'Endothelial', 'Endothelial', 'Endothelial', 'Endothelial', 'Pericyte',  'Pericyte', 'Endothelial', 'Mesenchymal', 'Mesenchymal', 'Mesenchymal', 'Mesenchymal', 'Mesenchymal', 'Neural', 'Neural',  'Erythroblast')

ref_d7$CelltypeReduced <- plyr::mapvalues(x = ref_d7$Celltype, from = cell_types_all, to = cell_types_reduced)

cellDefectReduced <- paste0(ref_d7@meta.data$CelltypeReduced, "_", ref_d7@meta.data$defect)
m <- data.frame("Celltype_Defect_Reduced" = cellDefectReduced)
rownames(m) <- rownames(ref_d7@meta.data)
ref_d7 <- AddMetaData(object = ref_d7, metadata = m)

saveRDS(ref_d7, file = "Data/vml_cca_d7.rds")
#### 3a. Load ST Data and Create List of Seurat Objects, Perform Seurat Label Transfer ####
st_list <- sapply(samples, USE.NAMES = TRUE, function(geo) {
  print(geo)
  st_se <- Load10X_Spatial(paste0(getwd(), sprintf("/10xOutputData/%s/%s", sample_dict[[geo]], geo)), 
                           assay = "Spatial", filter.matrix = TRUE, slice = sprintf("%s", geo), 
                           filename = "filtered_feature_bc_matrix.h5")
  # remove spots with <1 gene detected
  st_se <- subset(st_se, subset = nFeature_Spatial > 0)
  st_se$orig.ident <- geo
  st_se$Zone <- read.csv(paste0(getwd(), sprintf("/10xOutputData/%s/%s/Zones.csv", sample_dict[[geo]], geo)), row.names = 1)
  st_se <- SCTransform(st_se, assay = "Spatial", verbose = FALSE)
  st_se <- st_se %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
    FindNeighbors(reduction = "pca", dims = 1:30) %>% 
    FindClusters(verbose = FALSE) %>% 
    RunUMAP(reduction = "pca", dims = 1:30)
  # Seurat Label Transfer
  anchors <- FindTransferAnchors(reference = ref_d7, query = st_se, 
                                 normalization.method = "SCT")
  predictions.assay <- TransferData(anchorset = anchors, 
                                    refdata = ref_d7$Celltype_Defect_Reduced, prediction.assay = TRUE,
                                    weight.reduction = st_se[["pca"]], dims = 1:30,
                                    k.weight = 45)
  st_se[["predictions"]] <- predictions.assay
  return(list(st_se))
})
saveRDS(st_list, file = "Data/Mus_ST_list_label_transfer.RDS")

#### 3b. Plot PBS ####
st_se <- st_list[["PBS_1197L"]]

DefaultAssay(st_se) <- "SCT"

pdf("Plots/Finalized/PBS_1197L_umis_highres.pdf")
SpatialFeaturePlot(st_se, features = "nCount_Spatial", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/PBS_1197L_zones_highres.pdf")
SpatialDimPlot(st_se, group.by = "Zone", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/PBS_1197L_gene_overlays.pdf", onefile = TRUE)
SpatialFeaturePlot(st_se, features = "Myh1", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Aspn", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Myog", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Cd68", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/PBS_1197L_gene_overlays2.pdf", onefile = TRUE)
SpatialFeaturePlot(st_se, features = "Myh4", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Col1a1", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Myod1", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Ptprc", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Myh3", pt.size.factor = 1, crop = FALSE)
SpatialFeaturePlot(st_se, features = "Col1a1", pt.size.factor = 1, crop = FALSE)
dev.off()

DefaultAssay(st_se) <- "predictions"
pdf("Plots/Finalized/PBS_1197L_macrophages_highres.pdf")
SpatialFeaturePlot(st_se, features = "Macrophage", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1), min.cutoff = 0, max.cutoff = 1)
dev.off()

pdf("Plots/Finalized/PBS_1197L_tenocyte_highres.pdf")
SpatialFeaturePlot(st_se, features = "Tenocyte", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1), min.cutoff = 0, max.cutoff = 1)
dev.off()

pdf("Plots/Finalized/PBS_1197L_musc_highres.pdf")
SpatialFeaturePlot(st_se, features = "Satellite-Cell", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1), min.cutoff = 0, max.cutoff = 0.4)
dev.off()

pdf("Plots/Finalized/PBS_1197L_myonuclei_highres.pdf")
SpatialFeaturePlot(st_se, features = "Myonuclei", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1), min.cutoff = 0, max.cutoff = 0.4)
dev.off()

#### 3c. Plot ITD1 ####
st_se <- st_list[["ITD1_1203R"]]
DefaultAssay(st_se) <- "SCT"

pdf("Plots/Finalized/ITD1_1203R_umis_highres.pdf")
SpatialFeaturePlot(st_se, features = "nCount_Spatial", pt.size.factor = 1, crop = FALSE)
dev.off()

pdf("Plots/Finalized/ITD1_1203R_zones_highres.pdf")
SpatialDimPlot(st_se, group.by = "Zone", pt.size.factor = 1, crop = FALSE)
dev.off()

DefaultAssay(st_se) <- "predictions"
pdf("Plots/Finalized/ITD1_1203R_macrophages_highres.pdf")
SpatialFeaturePlot(st_se, features = "Macrophage", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1), min.cutoff = 0, max.cutoff = 1)
dev.off()

pdf("Plots/Finalized/ITD1_1203R_tenocyte_highres.pdf")
SpatialFeaturePlot(st_se, features = "Tenocyte", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1), min.cutoff = 0, max.cutoff = 1)
dev.off()

pdf("Plots/Finalized/ITD1_1203R_musc_highres.pdf")
SpatialFeaturePlot(st_se, features = "Satellite-Cell", pt.size.factor = 1,
                   crop = FALSE, alpha = c(0.1,1), min.cutoff = 0, max.cutoff = 0.4)
dev.off()
#### 4a. Generate seurat object with all ST datasets ####
st_1 <- st_list[[samples[1]]]; st_2 <- st_list[[samples[2]]]; st_3 <- st_list[[samples[3]]]
st_4 <- st_list[[samples[4]]]; st_5 <- st_list[[samples[5]]]; st_6 <- st_list[[samples[6]]]
st_7 <- st_list[[samples[7]]]; st_8 <- st_list[[samples[8]]]

st_1$treatment <- "PBS"; st_2$treatment <- "ITD1"; st_3$treatment <- "ITD1"
st_4$treatment <- "PBS"; st_5$treatment <- "PBS"; st_6$treatment <- "PBS"
st_7$treatment <- "ITD1"; st_8$treatment <- "ITD1"

st_1$mouseID <- "599"; st_2$mouseID <- "599"; st_3$mouseID <- "600"
st_4$mouseID <- "600"; st_5$mouseID <- "1197"; st_6$mouseID <- "1203"
st_7$mouseID <- "1203"; st_8$mouseID <- "1202"

st_1$round <- "1"; st_2$round <- "1"; st_3$round <- "1"; st_4$round <- "1"
st_5$round <- "2"; st_6$round <- "2"; st_7$round <- "2"; st_8$round <- "2"

st_all <- merge(st_1, c(st_2, st_3, st_4, st_5, st_6, st_7, st_8))
table(st_all$orig.ident, st_all$Zone)

Idents(st_all) <- 'Zone'
st_filt <- subset(st_all, ident = c('Defect', 'IntactMuscle', 'Transition'))
table(st_filt$orig.ident, st_filt$Zone)

st_all <- SCTransform(st_all, assay = "Spatial", verbose = FALSE)
st_all <- st_all %>% RunPCA(assay = "SCT", verbose = FALSE) %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30)

saveRDS(st_all, file = "Data/merged_filtered_599L_599R_600L_600R_1197L_1203L_1203R_1202L.RDS")

#### 4b. Plot (UMAPs) ####
pdf(file="Plots/all_umap_dimplot.pdf", onefile = TRUE)
DimPlot(st_all, reduction = "umap", label = TRUE, pt.size = 3) + NoLegend()
DimPlot(st_all, reduction = "umap", group.by = 'orig.ident', label = FALSE, pt.size = 3) + NoLegend()
DimPlot(st_all, reduction = "umap", group.by = 'orig.ident', label = FALSE, pt.size = 3)
DimPlot(st_all, reduction = "umap", group.by = 'Zone', label = FALSE, pt.size = 3)
DimPlot(st_all, reduction = "umap", group.by = 'treatment', label = FALSE, pt.size = 3)
SpatialDimPlot(st_all, label = TRUE, group.by = 'seurat_clusters', label.size = 3, pt.size.factor = 1, ncol=2) & NoLegend()
dev.off()

pdf(file="Plots/all_umap_qc_filtered.pdf", onefile = TRUE)
FeaturePlot(st_all, reduction = "umap", features = 'nCount_Spatial', pt.size = 1)
FeaturePlot(st_all, reduction = "umap", features = 'nFeature_Spatial', pt.size = 1)
DimHeatmap(st_all, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

DefaultAssay(st_all) <- 'SCT'
Idents(st_all) <- 'Celltype'
cell.degs <- FindAllMarkers(st_all, only.pos = TRUE)
top10 <- cell.degs %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("Plots/Finalized/Celltype_Markergenes_Dotplot.pdf", width = 8, height = 4)
DotPlot(st_all, features = c('Aspn', 'Postn', 'Acta2', 'Tnc', 'Myod1', 'Myog', 'Myf5', 'Myh4', 'Myh1', 'Ctss', 'Ptprc', 'Cd68', 'Adgre1'), group.by = 'Celltype', idents = c('Macrophage', 'Tenocyte', 'Satellite-Cell', 'Myonuclei')) + RotatedAxis()
dev.off()

pdf("Plots/mouse_051822_overlays.pdf", onefile = TRUE, width = 20, height = 5)
SpatialFeaturePlot(st_all, features = c("Axin2"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0)
SpatialFeaturePlot(st_all, features = c("Hic1"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0)
SpatialFeaturePlot(st_all, features = c("Ly6a"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0)
SpatialFeaturePlot(st_all, features = c("Thy1"), pt.size.factor = 1, crop = FALSE, min.cutoff = 0)
dev.off()

#### 5a. PBS treated tissue plots/analysis ####
Idents(st_all) <- 'treatment'
st_pbs <- subset(st_all, idents = 'PBS')
st_pbs <- ScaleData(st_pbs, assay = 'predictions')

zone_cell_avgs_pbs <- AverageExpression(st_pbs, assays = 'predictions', group.by = c('Zone', 'orig.ident'), slot = 'data') %>% as.data.frame()
zone_cell_avgs_pbs_melted <- zone_cell_avgs_pbs %>% melt()
zone_cell_avgs_pbs_melted$Celltype <- rep(rownames(zone_cell_avgs_pbs),12)
zone_cell_avgs_pbs_melted$Zone <- c(rep('Defect', 96), rep('IntactMuscle', 96), rep('Transition', 96))
zone_cell_avgs_pbs_melted$Tissue <- rep(c(rep('1197L', 24), rep('1203L', 24), rep('599L', 24), rep('600R', 24)), 3)
write.csv(zone_cell_avgs_pbs_melted, file = 'Data/PBS_Zone_Celltype_Averages.csv')

png("Plots/Finalized/PBS_Macrophages_byZone.png", res=200, units="in", width=3, height=6)
p <- zone_cell_avgs_pbs_melted %>% 
  filter(Celltype == 'Macrophage') %>%
  ggplot(aes(x=Zone, y=value, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.8)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45)
p
dev.off()

png("Plots/Finalized/PBS_MDC_byZone.png", res=200, units="in", width=3, height=6)
p <- zone_cell_avgs_pbs_melted %>% 
  filter(Celltype == 'Tenocyte') %>%
  ggplot(aes(x=Zone, y=value, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.8)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none")
p
dev.off()

png("Plots/Finalized/PBS_MuSC_byZone.png", res=200, units="in", width=3, height=6)
p <- zone_cell_avgs_pbs_melted %>% 
  filter(Celltype == 'Satellite-Cell') %>%
  ggplot(aes(x=Zone, y=value, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none")
p
dev.off()

# PATHWAY ANALYSIS ON ZONES
Idents(st_pbs) <- 'Zone'
DefaultAssay(st_pbs) <- 'SCT'
all.genes <- rownames(st_pbs)
st_pbs <- ScaleData(st_pbs, features = all.genes)
# remove mitochondrial and ribosomal genes
genes <- !grepl(pattern = "^Rp[l|s]|mt", x = rownames(st_pbs))
genes <- rownames(st_pbs)[genes]
de_markers <- FindAllMarkers(st_pbs, only.pos = TRUE, assay = 'SCT', features = genes)
de_markers$entrez <- mapIds(org.Mm.eg.db, rownames(de_markers),'ENTREZID', 'SYMBOL')
top10 <- de_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf(file="Plots/PBS_Zone_Markergenes.pdf", width = 11, height = 8)
DoHeatmap(st_pbs, features = unique(top10$gene), slot = 'scale.data') +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(text = element_text(color = 'black'))
dev.off()

## GOTerm analysis
go_concat <- NULL
for (i in c('Defect', 'IntactMuscle', 'Transition')){
  go_temp <- enrichGO(de_markers$entrez[de_markers$cluster == i],
                      OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, qvalueCutoff = 1)
  go_concat <- c(go_concat, go_temp)
}

pdf(file="Plots/Finalized/PBS_enrichGO_transition_zone.pdf", onefile = TRUE, width = 6, height = 3)
#dotplot(go_concat[[3]], showCategory=15)
dotplot(go_concat[[3]], showCategory=c('muscle structure development',
                                       'striated muscle tissue development',
                                       'muscle cell differentiation',
                                       'striated muscle contraction',
                                       'muscle adaptation'))
dev.off()
pdf(file="Plots/Finalized/PBS_enrichGO_defect_zone.pdf", onefile = TRUE, width = 6, height = 3)
#dotplot(go_concat[[1]], showCategory=20)
dotplot(go_concat[[1]], showCategory=c('response to stress',
                                       'cellular response to chemical stimulus',
                                       'positive regulation of metabolic process',
                                       'immune system process',
                                       'regulation of signal transduction'))

dev.off()

pdf(file="Plots/Finalized/PBS_enrichGO_intact_zone.pdf", onefile = TRUE, width = 6, height = 3)
#dotplot(go_concat[[2]], showCategory=20)
dotplot(go_concat[[2]], showCategory=c('generation of precursor metabolites and energy', 
                                       'energy derivation by oxidation of organic compounds',
                                       'muscle structure development',
                                       'oxidative phosphorylation',
                                       'ATP metabolic process'))
dev.off()




#### 5b. ITD1 treated tissue plots/analysis ####
Idents(st_all) <- 'treatment'
st_itd1 <- subset(st_all, idents = 'ITD1')
st_itd1 <- ScaleData(st_itd1, assay = 'predictions')

zone_cell_avgs_itd1 <- AverageExpression(st_itd1, assays = 'predictions', group.by = c('Zone', 'orig.ident'), slot = 'data') %>% as.data.frame() #%>% melt()
zone_cell_avgs_itd1_melted <- zone_cell_avgs_itd1 %>% melt()
zone_cell_avgs_itd1_melted$Celltype <- rep(rownames(zone_cell_avgs_itd1),12)
zone_cell_avgs_itd1_melted$Zone <- c(rep('Defect', 96), rep('IntactMuscle', 96), rep('Transition', 96))
zone_cell_avgs_itd1_melted$Tissue <- rep(c(rep('1202L', 24), rep('1203R', 24), rep('599R', 24), rep('600L', 24)), 3)
write.csv(zone_cell_avgs_itd1_melted, file = 'Data/ITD1_Zone_Celltype_Averages.csv')

png("Plots/Finalized/ITD1_Macrophages_byZone.png", res=200, units="in", width=6, height=6)
p <- zone_cell_avgs_itd1_melted %>% 
  filter(Celltype == 'Macrophage') %>%
  ggplot(aes(x=Zone, y=value, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, aes(color = Tissue)) +
  labs(x="", y="Average Macrophage Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.4)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none")
p
dev.off()

png("Plots/Finalized/ITD1_MDC_byZone.png", res=200, units="in", width=6, height=6)
p <- zone_cell_avgs_itd1_melted %>% 
  filter(Celltype == 'Tenocyte') %>%
  ggplot(aes(x=Zone, y=value, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, aes(color = Tissue)) +
  labs(x="", y="Average Myofibroblast Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none")
p
dev.off()

png("Plots/Finalized/ITD1_MuSC_byZone.png", res=200, units="in", width=6, height=6)
p <- zone_cell_avgs_itd1_melted %>% 
  filter(Celltype == 'Satellite-Cell') %>%
  ggplot(aes(x=Zone, y=value, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, aes(color = Tissue)) +
  labs(x="", y="Average MuSC Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="right")
p
dev.off()

# Compare PBS vs ITD1
zone_cell_avgs_both <- rbind(zone_cell_avgs_itd1_melted, zone_cell_avgs_pbs_melted)
zone_cell_avgs_both$Treatment <- c(rep('ITD1', 288), rep('PBS', 288))

png("Plots/Finalized/Comparison_MuSC_byZone.png", res=200, units="in", width=6, height=6)
pD <- zone_cell_avgs_both %>% 
  filter(Zone == 'Defect') %>% filter(Celltype == 'Satellite-Cell') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average MuSC Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  ggtitle('Defect')

pT <- zone_cell_avgs_both %>% 
  filter(Zone == 'Transition') %>% filter(Celltype == 'Satellite-Cell') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average MuSC Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- zone_cell_avgs_both %>% 
  filter(Zone == 'IntactMuscle') %>% filter(Celltype == 'Satellite-Cell') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average MuSC Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

png("Plots/Finalized/Comparison_MDC_byZone.png", res=200, units="in", width=6, height=6)
pD <- zone_cell_avgs_both %>% 
  filter(Zone == 'Defect') %>% filter(Celltype == 'Tenocyte') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Myofibroblast Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Defect')

pT <- zone_cell_avgs_both %>% 
  filter(Zone == 'Transition') %>% filter(Celltype == 'Tenocyte') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Myofibroblast Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- zone_cell_avgs_both %>% 
  filter(Zone == 'IntactMuscle') %>% filter(Celltype == 'Tenocyte') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Myofibroblast Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

png("Plots/Finalized/Comparison_Macrophage_byZone.png", res=200, units="in", width=6, height=6)
pD <- zone_cell_avgs_both %>% 
  filter(Zone == 'Defect') %>% filter(Celltype == 'Macrophage') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Macrophage Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.8)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Defect')

pT <- zone_cell_avgs_both %>% 
  filter(Zone == 'Transition') %>% filter(Celltype == 'Macrophage') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Macrophage Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.8)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- zone_cell_avgs_both %>% 
  filter(Zone == 'IntactMuscle') %>% filter(Celltype == 'Macrophage') %>%
  filter(Tissue != '599R') %>%
  ggplot(aes(x=Treatment, y=value, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Treatment)) +
  geom_point(aes(fill=Treatment),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Average Macrophage Prediction Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,0.8)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

#### 6. Pathway differences by treatment ####
# Defect Zone
Idents(st_all) <- 'Zone'
st_defect <- subset(st_all, idents = 'Defect')
Idents(st_defect) <- 'treatment'
DefaultAssay(st_defect) <- 'SCT'
st_defect <- ScaleData(st_defect, features = rownames(st_defect))
genes <- !grepl(pattern = "^Rp[l|s]|mt", x = rownames(st_defect))
genes <- rownames(st_defect)[genes]
de_markers_defect <- FindAllMarkers(st_defect, only.pos = T, assay = 'SCT', features = genes)

mac_genes <- degs.mac$gene
de_markers_defect <- FindAllMarkers(st_defect, only.pos = T, assay = 'SCT', features = degs.mac$gene, logfc.threshold = 0, min.pct = 0.01, return.thresh = 0)

de_markers_defect$entrez <- mapIds(org.Mm.eg.db, rownames(de_markers_defect),'ENTREZID', 'SYMBOL')
top20 <- de_markers_defect %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf(file="Plots/Defect_Treatment_Markergenes.pdf", width = 11, height = 8)
DoHeatmap(st_defect, features = unique(top20$gene), slot = 'scale.data') +
  theme(text = element_text(color = 'black'))
dev.off()

## GOTerm analysis defect zone
go_defect <- enrichGO(de_markers_defect$entrez[de_markers_defect$cluster == 'ITD1'],
                      OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, qvalueCutoff = 1)

pdf(file="Plots/Finalized/Defect_Treatment_enrichGO.pdf", width = 6, height = 5)
dotplot(go_defect)
dev.off()

# Transition Zone
Idents(st_all) <- 'Zone'
st_trans <- subset(st_all, idents = 'Transition')
Idents(st_trans) <- 'treatment'
DefaultAssay(st_trans) <- 'SCT'
st_trans <- ScaleData(st_trans, features = rownames(st_trans))
genes <- !grepl(pattern = "^Rp[l|s]|mt", x = rownames(st_trans))
genes <- rownames(st_trans)[genes]
de_markers_trans <- FindAllMarkers(st_trans, only.pos = T, assay = 'SCT', features = genes)
de_markers_trans$entrez <- mapIds(org.Mm.eg.db, rownames(de_markers_trans),'ENTREZID', 'SYMBOL')
top10 <- de_markers_trans %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf(file="Plots/TransitionZone_Treatment_Markergenes.pdf", width = 11, height = 8)
DoHeatmap(st_trans, features = unique(top10$gene), slot = 'scale.data') +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme(text = element_text(color = 'black'))
dev.off()

## GOTerm analysis defect zone
go_trans <- enrichGO(de_markers_trans$entrez[de_markers_trans$cluster == 'ITD1'],
                      OrgDb = "org.Mm.eg.db", ont = "ALL", maxGSSize = 5000, qvalueCutoff = 1)

pdf(file="Plots/Finalized/TransitionZone_Treatment_enrichGO.pdf", width = 6, height = 5)
dotplot(go_trans)
dev.off()

#### 7. CellChat ####
library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)
setwd("/nas/homes/jlarouche/TA_VML/spGEX")

ref_d7 <- readRDS(file = "Data/vml_cca_d7.rds")

cellchat <- createCellChat(object = ref_d7, group.by = "Celltype")

CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
saveRDS(cellchat, file = "Data/cellchat.RDS")

pdf("Plots/Chord_LR_MuSC_MDC_Mac.pdf", height = 10, width = 10)
netVisual_chord_gene(cellchat, sources.use = c('Macrophage', 'Tenocyte', 'Satellite_Cell'), 
                     targets.use = c('Macrophage', 'Tenocyte', 'Satellite_Cell'), 
                     lab.cex = 0.5, legend.pos.y = 30,
                     small.gap = 1)
dev.off()

#### 8. CellChat Gene Modules ####
st_se <- st_list[["PBS_1197L"]]
DefaultAssay(st_se) <- "SCT"

fap_r <- c('Cd44', 'Cd47', 'Ncl', 'Sdc2', 'Lrp1', 'Sdc1', 'Sdc4', 'Fgfr1', 'Itga5', 'Itgb1', 'Cdh11', 'Itga11', 'Itga4', 'Itgb7', 'Axl', 'Egfr', 'Dag1', 'Lrp1', 'Fgfr1', 'Ncam1', 'Pdgfrb')
fap_l <- c('Col1a1', 'Col1a2', 'Fn1', 'Thbs4', 'Ptn', 'Col6a1', 'Col6a2', 'Thbs2', 'App', 'Mif', 'Mdk', 'Hspg2', 'Tnc', 'Lamb1', 'Comp', 'Thbs1', 'Angptl1', 'Thbs3', 'Angptl2', 'Col6a3', 'Lama4', 'Lamc1', 'Spp1', 'Tnn', 'Angptl4', 'Pdgfa', 'Ncam1', 'Hbegf', 'Gas6')

mac_l <- c('Thbs1', 'App', 'Mif', 'Ccl6', 'Fn1', 'Ccl9', 'Lgals9', 'Itga4', 'Itgb1', 'Tnf', 'Itgb7', 'Nampt', 'Ccl3', 'Ccl2', 'Selpg', 'Sema4a')
mac_r <- c('Cd74', 'Ccr2', 'Ccr1', 'Cd44', 'Sdc4', 'Ptprc', 'Ighm', 'Cxcr4', 'Sell', 'Plxnb2', 'Cd47', 'Tnfrsf1a', 'Tnfrsf1b', 'Pirb', 'Sdc3', 'Lrp1', 'Ncl')

musc_l <- c('Jam3', 'Col4a2', 'Lama2', 'Cdh15', 'Col4a1')
musc_r <- c('Itga7', 'Cdh15', 'Jam3', 'Vcam1', 'Igf1r')

st_se <- AddModuleScore(object = st_se, features = list(fap_l), ctrl = 5, 
                        name = 'FAP_Ligands')
st_se <- AddModuleScore(object = st_se, features = list(fap_r), ctrl = 5, 
                        name = 'FAP_Receptors')

st_se <- AddModuleScore(object = st_se, features = list(mac_l), ctrl = 5, 
                        name = 'Mø_Ligands')
st_se <- AddModuleScore(object = st_se, features = list(mac_r), ctrl = 5, 
                        name = 'Mø_Receptors')

st_se <- AddModuleScore(object = st_se, features = list(musc_l), ctrl = 5, 
                        name = 'MuSC_Ligands')
st_se <- AddModuleScore(object = st_se, features = list(musc_r), ctrl = 5, 
                        name = 'MuSC_Receptors')


pdf("Plots/Finalized/Overlays_LR_Scores_vertical.pdf", height = 10, width = 7.5)
SpatialFeaturePlot(st_se, features = c('Mø_Ligands1', 'Mø_Receptors1', 
                                       'FAP_Ligands1', 'FAP_Receptors1',
                                       'MuSC_Ligands1', 'MuSC_Receptors1'), 
                   min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 2.5, alpha = c(0,1),
                   ncol = 2)
dev.off()


DefaultAssay(st_all) <- "SCT"

st_all <- AddModuleScore(object = st_all, features = list(fap_l), ctrl = 5, 
                        name = 'FAP_Ligands')
st_all <- AddModuleScore(object = st_all, features = list(fap_r), ctrl = 5, 
                        name = 'FAP_Receptors')

st_all <- AddModuleScore(object = st_all, features = list(mac_l), ctrl = 5, 
                        name = 'Mø_Ligands')
st_all <- AddModuleScore(object = st_all, features = list(mac_r), ctrl = 5, 
                        name = 'Mø_Receptors')

st_all <- AddModuleScore(object = st_all, features = list(musc_l), ctrl = 5, 
                        name = 'MuSC_Ligands')
st_all <- AddModuleScore(object = st_all, features = list(musc_r), ctrl = 5, 
                        name = 'MuSC_Receptors')


Idents(st_all) <- 'Zone'
my_levels <- c('Defect', 'Transition', 'IntactMuscle')
Idents(st_all) <- factor(Idents(st_all), levels= my_levels)

pdf("Plots/Mouse_Violin_LRModules_byZone.pdf", height = 8, width = 6)
VlnPlot(st_all, features = c('Mø_Ligands1', 'Mø_Receptors1', 
                             'FAP_Ligands1', 'FAP_Receptors1',
                              'MuSC_Ligands1', 'MuSC_Receptors1'), 
        pt.size = 0, split.by = 'treatment', ncol = 2) + theme(legend.position = 'right')
dev.off()


st_all <- AddModuleScore(object = st_all, features = list(c(fap_l, fap_r)), ctrl = 5, 
                         name = 'FAP_LR')

st_all <- AddModuleScore(object = st_all, features = list(c(mac_l, mac_r)), ctrl = 5, 
                         name = 'Mø_LR')

st_all <- AddModuleScore(object = st_all, features = list(c(musc_l, musc_r)), ctrl = 5, 
                         name = 'MuSC_LR')


pdf("Plots/Mouse_Zone_UMAP.pdf", height = 5, width = 6)
DimPlot(st_all, group.by = 'Zone', pt.size = 2)
dev.off()

LR_byZone <- data.frame(matrix(ncol=9,nrow=24, dimnames=list(NULL, c('Tissue', 'Zone', 'Treatment', 'MuSC_L', 'MuSC_R', 'FAP_L', 'FAP_R', 'Mac_L', 'Mac_R'))))
i = 1
for (zone in c('Defect', 'Transition', 'IntactMuscle')){
  for (tissue in c('ITD1_1202L', 'ITD1_1203R',  'ITD1_599R',  'ITD1_600L',  'PBS_1197L',  'PBS_1203L',   'PBS_599L',   'PBS_600R')){
    LR_byZone$Tissue[i] <- tissue
    LR_byZone$Zone[i] <- zone
    LR_byZone$Treatment[i] <- strsplit(tissue, '_')[[1]][1]
    LR_byZone$MuSC_L[i] <- mean(st_all$MuSC_Ligands1[((st_all$orig.ident == tissue) & (st_all$Zone == zone))])
    LR_byZone$MuSC_R[i] <- mean(st_all$MuSC_Receptors1[((st_all$orig.ident == tissue) & (st_all$Zone == zone))])
    LR_byZone$Mac_L[i] <- mean(st_all$Mø_Ligands1[((st_all$orig.ident == tissue) & (st_all$Zone == zone))])
    LR_byZone$Mac_R[i] <- mean(st_all$Mø_Receptors1[((st_all$orig.ident == tissue) & (st_all$Zone == zone))])
    LR_byZone$FAP_L[i] <- mean(st_all$FAP_Ligands1[((st_all$orig.ident == tissue) & (st_all$Zone == zone))])
    LR_byZone$FAP_R[i] <- mean(st_all$FAP_Receptors1[((st_all$orig.ident == tissue) & (st_all$Zone == zone))])
    i <- i+1
  }
}

png("Plots/Finalized/PBS_Macrophage_Ligand_Score_byZone.png", res=200, units="in", width=3, height=6)
p <- LR_byZone %>%
  filter(Treatment == 'PBS') %>%
  ggplot(aes(x=Zone, y=Mac_L, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Mø Ligand Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.2,0.8)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45)
p
dev.off()

png("Plots/Finalized/PBS_Macrophage_Receptor_Score_byZone.png", res=200, units="in", width=3, height=6)
p <- LR_byZone %>%
  filter(Treatment == 'PBS') %>%
  ggplot(aes(x=Zone, y=Mac_R, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="Mø Receptor Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.2,0.8)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45)
p
dev.off()

png("Plots/Finalized/PBS_FAP_Ligand_Score_byZone.png", res=200, units="in", width=3, height=6)
p <- LR_byZone %>%
  filter(Treatment == 'PBS') %>%
  ggplot(aes(x=Zone, y=FAP_L, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="MDC Ligand Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,1)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45)
p
dev.off()

png("Plots/Finalized/PBS_FAP_Receptor_Score_byZone.png", res=200, units="in", width=3, height=6)
p <- LR_byZone %>%
  filter(Treatment == 'PBS') %>%
  ggplot(aes(x=Zone, y=FAP_R, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="MDC Receptor Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.4,0.4)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45)
p
dev.off()

png("Plots/Finalized/PBS_MuSC_Ligand_Score_byZone.png", res=200, units="in", width=3, height=6)
p <- LR_byZone %>%
  filter(Treatment == 'PBS') %>%
  ggplot(aes(x=Zone, y=MuSC_L, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="MuSC Ligand Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45)
p
dev.off()

png("Plots/Finalized/PBS_MuSC_Receptor_Score_byZone.png", res=200, units="in", width=3, height=6)
p <- LR_byZone %>%
  filter(Treatment == 'PBS') %>%
  ggplot(aes(x=Zone, y=MuSC_R, fill=Zone)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_point(aes(fill=Zone),size=4,shape=21, position = position_dodge(0.2)) +
  labs(x="", y="MuSC Receptor Score") +
  theme_classic() +
  scale_x_discrete(limits = c("Defect", "Transition", "IntactMuscle")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.4,0.4)) +
  stat_compare_means(comparisons = 
                       list(c("Defect", "Transition"), 
                            c("Transition", "IntactMuscle"),
                            c("Defect", "IntactMuscle")), 
                     bracket.size = 1, 
                     method = "t.test", paired = TRUE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.signif", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45)
p
dev.off()

png("Plots/MuSC_Ligands_byZone.png", res=200, units="in", width=6, height=6)
pD <- LR_byZone %>% 
  filter(Zone == 'Defect') %>%
  filter(! Tissue %in% c('ITD1_599R', 'PBS_1197L')) %>%
  ggplot(aes(x=Treatment, y=MuSC_L, fill = Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="MuSC Ligand Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  ggtitle('Defect')

pT <- LR_byZone %>% 
  filter(Zone == 'Transition') %>%
  ggplot(aes(x=Treatment, y=MuSC_L, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="MuSC Ligand Score (%)") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- LR_byZone %>% 
  filter(Zone == 'IntactMuscle') %>%
  ggplot(aes(x=Treatment, y=MuSC_L, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="MuSC Ligand Score (%)") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

png("Plots/MuSC_Receptors_byZone.png", res=200, units="in", width=6, height=6)
pD <- LR_byZone %>% 
  filter(Zone == 'Defect') %>%
  filter(! Tissue %in% c('ITD1_599R', 'PBS_1197L')) %>%
  ggplot(aes(x=Treatment, y=MuSC_R, fill = Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="MuSC Receptor Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  ggtitle('Defect')

pT <- LR_byZone %>% 
  filter(Zone == 'Transition') %>%
  ggplot(aes(x=Treatment, y=MuSC_R, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- LR_byZone %>% 
  filter(Zone == 'IntactMuscle') %>%
  ggplot(aes(x=Treatment, y=MuSC_R, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

png("Plots/Mac_Ligands_byZone.png", res=200, units="in", width=6, height=6)
pD <- LR_byZone %>% 
  filter(Zone == 'Defect') %>%
  filter(! Tissue %in% c('PBS_599L', 'ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=Mac_L, fill = Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="Mø Ligand Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  ggtitle('Defect')

pT <- LR_byZone %>% 
  filter(Zone == 'Transition') %>%
  filter(! Tissue %in% c('PBS_599L', 'ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=Mac_L, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- LR_byZone %>% 
  filter(Zone == 'IntactMuscle') %>%
  filter(! Tissue %in% c('PBS_599L', 'ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=Mac_L, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

png("Plots/Mac_Receptors_byZone.png", res=200, units="in", width=6, height=6)
pD <- LR_byZone %>% 
  filter(Zone == 'Defect') %>%
  filter(! Tissue %in% c('PBS_599L', 'ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=MuSC_R, fill = Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="Mø Receptor Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  ggtitle('Defect')

pT <- LR_byZone %>% 
  filter(Zone == 'Transition') %>%
  filter(! Tissue %in% c('PBS_599L', 'ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=MuSC_R, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- LR_byZone %>% 
  filter(Zone == 'IntactMuscle') %>%
  filter(! Tissue %in% c('PBS_599L', 'ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=MuSC_R, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()



png("Plots/FAP_Ligands_byZone.png", res=200, units="in", width=6, height=6)
pD <- LR_byZone %>% 
  filter(Zone == 'Defect') %>%
  filter(! Tissue %in% c('ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=FAP_L, fill = Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="FAP Ligand Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  ggtitle('Defect')

pT <- LR_byZone %>% 
  filter(Zone == 'Transition') %>%
  filter(! Tissue %in% c('ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=FAP_L, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- LR_byZone %>% 
  filter(Zone == 'IntactMuscle') %>%
  filter(! Tissue %in% c('ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=FAP_L, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.5,0.5)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

png("Plots/Mac_Receptors_byZone.png", res=200, units="in", width=6, height=6)
pD <- LR_byZone %>% 
  filter(Zone == 'Defect') %>%
  filter(! Tissue %in% c('ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=FAP_R, fill = Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  labs(x="", y="Mø Receptor Score") +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 20, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="none") +
  ggtitle('Defect')

pT <- LR_byZone %>% 
  filter(Zone == 'Transition') %>%
  filter(! Tissue %in% c('ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=FAP_R, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar") +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('Transition')

pI <- LR_byZone %>% 
  filter(Zone == 'IntactMuscle') %>%
  filter(! Tissue %in% c('ITD1_599R')) %>%
  ggplot(aes(x=Treatment, y=FAP_R, fill=Treatment)) + 
  stat_summary(fun.data="mean_se", geom="errorbar", width=0.5) +
  stat_summary(fun="mean", geom="crossbar", aes(fill=Zone)) +
  geom_jitter(width=0.1, size = 5, aes(color = Treatment)) +
  theme_classic() +
  scale_x_discrete(limits = c("PBS", "ITD1")) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(ylim=c(-0.25,0.25)) +
  stat_compare_means(comparisons = 
                       list(c("PBS", "ITD1")), 
                     bracket.size = 1, 
                     method = "t.test", paired = FALSE,
                     method.args = c(var.equal = FALSE), 
                     label = "p.format", size = 8) +
  theme(axis.title.y = element_text(size = 0, face = "plain", color = "black")) +
  theme(axis.text.y = element_text(size = 15, color = "black")) +
  theme(axis.text.x = element_text(size = 20, colour = "black")) +
  theme(legend.position="none") +
  rotate_x_text(angle = 45) +
  ggtitle('IntactMuscle')
pD + pT + pI
dev.off()

#### 9. Volcano Plots ####
library(EnhancedVolcano)
library(MAST)

Idents(st_all) <- 'Zone'
st_defect <- subset(st_all, idents = 'Defect')
st_trans <- subset(st_all, idents = 'Transition')

fap_r <- c('Cd44', 'Cd47', 'Ncl', 'Sdc2', 'Lrp1', 'Sdc1', 'Sdc4', 'Fgfr1', 'Itga5', 'Itgb1', 'Cdh11', 'Itga11', 'Itga4', 'Itgb7', 'Axl', 'Egfr', 'Dag1', 'Lrp1', 'Fgfr1', 'Ncam1', 'Pdgfrb')
fap_l <- c('Col1a1', 'Col1a2', 'Fn1', 'Thbs4', 'Ptn', 'Col6a1', 'Col6a2', 'Thbs2', 'App', 'Mif', 'Mdk', 'Hspg2', 'Tnc', 'Lamb1', 'Comp', 'Thbs1', 'Angptl1', 'Thbs3', 'Angptl2', 'Col6a3', 'Lama4', 'Lamc1', 'Spp1', 'Tnn', 'Angptl4', 'Pdgfa', 'Ncam1', 'Hbegf', 'Gas6')

mac_l <- c('Thbs1', 'App', 'Mif', 'Ccl6', 'Fn1', 'Ccl9', 'Lgals9', 'Itga4', 'Itgb1', 'Tnf', 'Itgb7', 'Nampt', 'Ccl3', 'Ccl2', 'Selpg', 'Sema4a')
mac_r <- c('Cd74', 'Ccr2', 'Ccr1', 'Cd44', 'Sdc4', 'Ptprc', 'Ighm', 'Cxcr4', 'Sell', 'Plxnb2', 'Cd47', 'Tnfrsf1a', 'Tnfrsf1b', 'Pirb', 'Sdc3', 'Lrp1', 'Ncl')

musc_l <- c('Jam3', 'Col4a2', 'Lama2', 'Cdh15', 'Col4a1')
musc_r <- c('Itga7', 'Cdh15', 'Jam3', 'Vcam1', 'Igf1r')

goi <- c(fap_r, fap_l, mac_l, mac_r, musc_l, musc_r)

Idents(st_trans) <- 'treatment'
deg.trans <- FindMarkers(st_trans,ident.1 = 'ITD1',test.use = 'MAST',logfc.threshold = 0)

Idents(st_defect) <- 'treatment'
deg.defect <- FindMarkers(st_defect,ident.1 = 'ITD1',test.use = 'MAST',logfc.threshold = 0)

goi_defect <- c('Ccl9', 'Sema4a', 'Gas6', 'Lgals9', 'Fn1', 'Tnfrsf1b', 'App', 'Sell', 'Ccl2', 'Ccl3', 'Ptprc', 'Ptn', 'Sdc3', 'Thbs4', 'Ncam1', 'Cdh15', 'Itga7', 'Dag1', 'Col1a1')

vol_d <- EnhancedVolcano(deg.defect,
                       lab = rownames(deg.defect),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = 'Defect Zone',
                       subtitle = '',
                       xlab = bquote(~Log[2]~ ("fold change")),
                       ylab = bquote(~-Log[10]("p-adj")),
                       axisLabSize = 18,
                       titleLabSize = 18,
                       caption = NULL,
                       #labvjust = -4,
                       selectLab = goi_defect,
                       ylim = c(0,5),
                       xlim = c(-.75, .75),
                       pCutoff = 0.05,
                       FCcutoff = 0.0585,
                       pointSize = 1.0,
                       labSize = 6.0,
                       boxedLabels = FALSE,
                       colAlpha = 0.5,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = TRUE,
                       widthConnectors = 0.2,
                       colConnectors = 'black',
                       col = c('#999999', '#009E73', '#56B4E9', '#E69F00'),
                       max.overlaps = 25)
vol_d

goi_trans <- c('Col6a3', 'Lrp1', 'Pirb', 'App', 'Itga11', 'Plxnb2', 'Gas6', 'Cd44', 'Lamb1', 'Tnfrsf1a','Myh1', 'Tpm1', 'Myh4', 'Myoz1', 'Ccl5', 'Nkg7')
vol_t <- EnhancedVolcano(deg.trans,
                       lab = rownames(deg.trans),
                       x = 'avg_log2FC',
                       y = 'p_val_adj',
                       title = 'Transition Zone',
                       subtitle = '',
                       xlab = bquote(~Log[2]~ ("fold change")),
                       ylab = bquote(~-Log[10]("p-adj")),
                       axisLabSize = 18,
                       titleLabSize = 18,
                       caption = NULL,
                       #labvjust = -4,
                       selectLab = goi_trans,
                       ylim = c(0,5),
                       xlim = c(-0.75, 0.75),
                       pCutoff = 0.05,
                       FCcutoff = 0.0585,
                       pointSize = 1.0,
                       labSize = 6.0,
                       boxedLabels = FALSE,
                       colAlpha = 0.5,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = TRUE,
                       widthConnectors = 0.2,
                       colConnectors = 'black',
                       col = c('#999999', '#009E73', '#56B4E9', '#E69F00'),
                       max.overlaps = 25)
vol_t

pdf("Plots/Finalized/volcano_defect_pbs_v_itd1_padjusted_051022.pdf", width = 10, height = 7)
vol_d
dev.off()

pdf("Plots/Finalized/volcano_transition_pbs_v_itd1_padjusted_051022.pdf", width = 10, height = 7)
vol_t
dev.off()
