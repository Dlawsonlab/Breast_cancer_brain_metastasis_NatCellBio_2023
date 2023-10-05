### Load in Libraries ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(Seurat)
library(DoubletFinder)
library(Seurat.utils)

#### Load in Data ####
BL6_D4.seurat <- readRDS("BL6_D4.rds")
BL6_D10.seurat <- readRDS("BL6_D10.rds")
RAG_D4.seurat <- readRDS("RAG_D4.rds")
RAG_D10.seurat <- readRDS("RAG_D10.rds")

###B6D4 Doublets Marked
BL6_D4.seurat <- NormalizeData(BL6_D4.seurat)
BL6_D4.seurat <- FindVariableFeatures(BL6_D4.seurat)
BL6_D4.seurat <- ScaleData(BL6_D4.seurat)
BL6_D4.seurat <- RunPCA(BL6_D4.seurat)
ElbowPlot(BL6_D4.seurat)
BL6_D4.seurat <- FindNeighbors(BL6_D4.seurat, dims = 1:18)
BL6_D4.seurat <- RunUMAP(BL6_D4.seurat, dims = 1:18)
BL6_D4.seurat <- subset(BL6_D4.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000
                        & percent.mt < 20)

sweep.res.B6D4 <- paramSweep_v3(BL6_D4.seurat, PCs = 1:18, sct = FALSE)
sweep.stats.B6D4 <- summarizeSweep(sweep.res.B6D4, GT = FALSE)

pdf("doubletsB6D4.pdf")
bcmvn.B6D4 <- find.pK(sweep.stats.B6D4)
dev.off()

# update expected doublet rate
nExp_poi.B6D4 <- round(0.061*nrow(BL6_D4.seurat@meta.data))
annotations.B6D4 <- BL6_D4.seurat@meta.data$ClusteringResults
homotypic.prop.B6D4 <- modelHomotypic(annotations.B6D4)

# update pK based on simulation
BL6_D4.seurat <- doubletFinder_v3(BL6_D4.seurat, PCs = 1:18, pN = 0.25, pK = 0.04, 
                                  nExp = nExp_poi.B6D4, reuse.pANN = FALSE, sct = FALSE)
singlet_or_doublet = colnames(BL6_D4.seurat@meta.data)[grepl("DF.classification", 
                                                             colnames(BL6_D4.seurat@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(BL6_D4.seurat, group.by = "orig.ident") + NoAxes(),
                   DimPlot(BL6_D4.seurat, group.by = singlet_or_doublet) + NoAxes())

VlnPlot(BL6_D4.seurat, features = "nFeature_RNA", group.by = singlet_or_doublet, pt.size = 0.1)
BL6_D4.seurat[["single_or_doublet"]] <- BL6_D4.seurat@meta.data$DF.classification

#update path
saveRDS(BL6_D4.seurat, "BL6_D4_doublets_marked.rds")

###B6D10 Doublets Marked
BL6_D10.seurat <- NormalizeData(BL6_D10.seurat)
BL6_D10.seurat <- FindVariableFeatures(BL6_D10.seurat)
BL6_D10.seurat <- ScaleData(BL6_D10.seurat)
BL6_D10.seurat <- RunPCA(BL6_D10.seurat)
ElbowPlot(BL6_D10.seurat)
BL6_D10.seurat <- FindNeighbors(BL6_D10.seurat, dims = 1:16)
BL6_D10.seurat <- RunUMAP(BL6_D10.seurat, dims = 1:16)
BL6_D10.seurat <- subset(BL6_D10.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 
                         & percent.mt < 20)

sweep.res.B6D10 <- paramSweep_v3(BL6_D10.seurat, PCs = 1:16, sct = FALSE)
sweep.stats.B6D10 <- summarizeSweep(sweep.res.B6D10, GT = FALSE)

pdf("doubletsB6D10.pdf")
bcmvn.B6D10 <- find.pK(sweep.stats.B6D10)
dev.off()

# update expected doublet rate
nExp_poi.B6D10 <- round(0.061*nrow(BL6_D10.seurat@meta.data))
annotations.B6D10 <- BL6_D10.seurat@meta.data$ClusteringResults
homotypic.prop.B6D10 <- modelHomotypic(annotations.B6D10)

# update pK based on simulation
BL6_D10.seurat <- doubletFinder_v3(BL6_D10.seurat, PCs = 1:16, pN = 0.25, pK = 0.13, 
                                   nExp = nExp_poi.B6D10, reuse.pANN = FALSE, sct = FALSE)
singlet_or_doublet = colnames(BL6_D10.seurat@meta.data)[grepl("DF.classification", 
                                                              colnames(BL6_D10.seurat@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(BL6_D10.seurat, group.by = "orig.ident") + NoAxes(),
                   DimPlot(BL6_D10.seurat, group.by = singlet_or_doublet) + NoAxes())

VlnPlot(BL6_D10.seurat, features = "nFeature_RNA", group.by = singlet_or_doublet, pt.size = 0.1)
BL6_D10.seurat[["single_or_doublet"]] <- BL6_D10.seurat@meta.data$DF.classification

#update path
saveRDS(BL6_D10.seurat, "BL6_D10_doublets_marked.rds")

###RAGD4 Doublets Marked
RAG_D4.seurat <- NormalizeData(RAG_D4.seurat)
RAG_D4.seurat <- FindVariableFeatures(RAG_D4.seurat)
RAG_D4.seurat <- ScaleData(RAG_D4.seurat)
RAG_D4.seurat <- RunPCA(RAG_D4.seurat)
ElbowPlot(RAG_D4.seurat)
RAG_D4.seurat <- FindNeighbors(RAG_D4.seurat, dims = 1:17)
RAG_D4.seurat <- RunUMAP(RAG_D4.seurat, dims = 1:17)
RAG_D4.seurat <- subset(RAG_D4.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 
                        & percent.mt < 20)

sweep.res.RAGD4 <- paramSweep_v3(RAG_D4.seurat, PCs = 1:17, sct = FALSE)
sweep.stats.RAGD4 <- summarizeSweep(sweep.res.RAGD4, GT = FALSE)

pdf("doubletsRAGD4.pdf")
bcmvn.RAGD4 <- find.pK(sweep.stats.RAGD4)
dev.off()

# update expected doublet rate
nExp_poi.RAGD4 <- round(0.054*nrow(RAG_D4.seurat@meta.data))
annotations.RAGD4 <- RAG_D4.seurat@meta.data$ClusteringResults
homotypic.prop.RAGD4 <- modelHomotypic(annotations.RAGD4)

# update pK based on simulation
RAG_D4.seurat <- doubletFinder_v3(RAG_D4.seurat, PCs = 1:17, pN = 0.25, pK = 0.02, 
                                  nExp = nExp_poi.RAGD4, reuse.pANN = FALSE, sct = FALSE)
singlet_or_doublet = colnames(RAG_D4.seurat@meta.data)[grepl("DF.classification", 
                                                             colnames(RAG_D4.seurat@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(RAG_D4.seurat, group.by = "orig.ident") + NoAxes(),
                   DimPlot(RAG_D4.seurat, group.by = singlet_or_doublet) + NoAxes())

VlnPlot(RAG_D4.seurat, features = "nFeature_RNA", group.by = singlet_or_doublet, pt.size = 0.1)
RAG_D4.seurat[["single_or_doublet"]] <- RAG_D4.seurat@meta.data$DF.classification

#update path
saveRDS(RAG_D4.seurat, "RAG_D4_doublets_marked.rds")

###RAGD10 Doublets Marked
RAG_D10.seurat <- NormalizeData(RAG_D10.seurat)
RAG_D10.seurat <- FindVariableFeatures(RAG_D10.seurat)
RAG_D10.seurat <- ScaleData(RAG_D10.seurat)
RAG_D10.seurat <- RunPCA(RAG_D10.seurat)
ElbowPlot(RAG_D10.seurat)
RAG_D10.seurat <- FindNeighbors(RAG_D10.seurat, dims = 1:16)
RAG_D10.seurat <- RunUMAP(RAG_D10.seurat, dims = 1:16)
RAG_D10.seurat <- subset(RAG_D10.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 
                         & percent.mt < 20)

sweep.res.RAGD10 <- paramSweep_v3(RAG_D10.seurat, PCs = 1:16, sct = FALSE)
sweep.stats.RAGD10 <- summarizeSweep(sweep.res.RAGD10, GT = FALSE)

pdf("doubletsRAGD10.pdf")
bcmvn.RAGD10 <- find.pK(sweep.stats.RAGD10)
dev.off()

# update expected doublet rate
nExp_poi.RAGD10 <- round(0.046*nrow(RAG_D10.seurat@meta.data))
annotations.RAGD10 <- RAG_D10.seurat@meta.data$ClusteringResults
homotypic.prop.RAGD10 <- modelHomotypic(annotations.RAGD10)

# update pK based on simulation
RAG_D10.seurat <- doubletFinder_v3(RAG_D10.seurat, PCs = 1:16, pN = 0.25, pK = 0.09, 
                                   nExp = nExp_poi.RAGD10, reuse.pANN = FALSE, sct = FALSE)
singlet_or_doublet = colnames(RAG_D10.seurat@meta.data)[grepl("DF.classification", 
                                                              colnames(RAG_D10.seurat@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(RAG_D10.seurat, group.by = "orig.ident") + NoAxes(),
                   DimPlot(RAG_D10.seurat, group.by = singlet_or_doublet) + NoAxes())

VlnPlot(RAG_D10.seurat, features = "nFeature_RNA", group.by = singlet_or_doublet, pt.size = 0.1)
RAG_D10.seurat[["single_or_doublet"]] <- RAG_D10.seurat@meta.data$DF.classification

#update path
saveRDS(RAG_D10.seurat, "RAG_D10_doublets_marked.rds")

### Merge Seurat objects into one ###
BL6vsRAG.seurat <- merge(BL6_D4.seurat, y = c(BL6_D10.seurat, RAG_D4.seurat, RAG_D10.seurat), 
                         add.cell.ids = c("BL6_D4", "BL6_D10", "RAG_D4", "RAG_D10"), project = "CD45_BL6_vs_RAG")
BL6vsRAG.seurat

###############################################
######################### start basic analysis
###############################################

BL6vsRAG.seurat <- NormalizeData(BL6vsRAG.seurat)
BL6vsRAG.seurat <- FindVariableFeatures(BL6vsRAG.seurat)
BL6vsRAG.seurat <- ScaleData(BL6vsRAG.seurat)
BL6vsRAG.seurat <- RunPCA(BL6vsRAG.seurat)
ElbowPlot(BL6vsRAG.seurat)
BL6vsRAG.seurat <- FindNeighbors(BL6vsRAG.seurat, dims = 1:18)
BL6vsRAG.seurat <- RunUMAP(BL6vsRAG.seurat, dims = 1:18)
BL6vsRAG.seurat <- FindClusters(BL6vsRAG.seurat, resolution = 0.2)

# add relevant metadata columns
BL6vsRAG.seurat@meta.data$mouse_strain <- sapply(
        strsplit(BL6vsRAG.seurat@meta.data$orig.ident, "_"), "[",1)
BL6vsRAG.seurat@meta.data$timepoint <- sapply(
        strsplit(BL6vsRAG.seurat@meta.data$orig.ident, "_"), "[",2)
BL6vsRAG.seurat@meta.data$mouse_id <- sapply(
        strsplit(BL6vsRAG.seurat@meta.data$orig.ident, "_"), "[",3)

# annotate top-level clusters
cluster_markers <- FindAllMarkers(BL6vsRAG.seurat, only.pos = T)

top10<- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(BL6vsRAG.seurat, features = top10$gene, group.by = "seurat_clusters") + NoLegend()

# UMAP colored by cluster
DimPlot(BL6vsRAG.seurat, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.4) & NoAxes()

# Cluster 9 Microglia + monocyte w/ overlapping Cldn5 -> Doublets
seurat.sub <- subset(BL6vsRAG.seurat, seurat_clusters != "9")

seurat.sub <- NormalizeData(seurat.sub)
seurat.sub <- FindVariableFeatures(seurat.sub)
seurat.sub <- ScaleData(seurat.sub)
seurat.sub <- RunPCA(seurat.sub)
ElbowPlot(seurat.sub)
seurat.sub <- FindNeighbors(seurat.sub, dims=1:15)
seurat.sub <- RunUMAP(seurat.sub, dims=1:15)
seurat.sub <- FindClusters(seurat.sub, res=0.2)

# annotate top-level clusters
cluster_markers_v2 <- FindAllMarkers(seurat.sub, only.pos = T)


top10_v2<- cluster_markers_v2 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat.sub, features = top10_v2$gene, group.by = "seurat_clusters") + NoLegend()

# UMAP colored by cluster
DimPlot(seurat.sub, group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.4) & NoAxes()

# annotate clusters based on markers
current.cluster.ids <- 0:11
new.cluster.ids <- c("Microglia","Monocytes","Microglia",
                     "T.cells", "NK.cells",
                     "Macrophages.DCs", "Monocytes","mDCs",
                     "Neutrophils","B.cells","pv.Macrophages","T.cells")
seurat.sub@meta.data$named_clusters <- plyr::mapvalues(x = seurat.sub@meta.data$seurat_clusters, 
                                                       from = current.cluster.ids, 
                                                       to = new.cluster.ids)

################### Full UMAP subfigures ################### 
############################################################

seurat.sub@meta.data$named_clusters <- factor(seurat.sub@meta.data$named_clusters,
                                              levels = c("Microglia","Monocytes",
                                                         "T.cells", "NK.cells",
                                                         "Macrophages.DCs","mDCs",
                                                         "Neutrophils","B.cells","pv.Macrophages"))


color_celltype_map <- c("B.cells" = "palevioletred3",
                        "Macrophages.DCs" = "burlywood2",
                        "mDCs" = "burlywood3",
                        "Microglia" = "#5D9D52FF",
                        "pv.Macrophages" = "#004949",
                        "Neutrophils" = "#B2DEF0FF",
                        "Monocytes" = "#b66dff",
                        "NK.cells" = "#DA0001FF", 
                        "T.cells" = "#F4988DFF")

# UMAP colored by cell type
DimPlot(seurat.sub, group.by = "named_clusters", pt.size = 0.4, label = TRUE, repel = TRUE, label.size = 11) + 
        scale_color_manual(values = color_celltype_map) & NoAxes() & NoLegend()

#Generate feature plots
FeaturePlot(seurat.sub, features = c("Cd3e", "Tmem119", "P2ry12", "Itgam", "Itgax", 
                                     "Ncr1", "Cd19", "Ly6c2", "Ly6g"))

# UMAP colored by cell type & split by strain/timepoint
seurat.sub$timepoint_strain <- paste0(seurat.sub$timepoint, "_",
                                      seurat.sub$mouse_strain)
seurat.sub@meta.data$timepoint_strain <- factor(seurat.sub@meta.data$timepoint_strain, 
                                                levels = c("D4_BL6","D10_BL6",
                                                           "D4_RAG","D10_RAG"))

# extract umap coordinates and add to metadata
umap_coords <- as.data.frame(seurat.sub@reductions$umap@cell.embeddings)
seurat.sub <- AddMetaData(seurat.sub, umap_coords)

ggplot(seurat.sub@meta.data,
       aes(x=UMAP_1, y=UMAP_2, color=named_clusters)) +  
        geom_point(size=0.1) + 
        facet_wrap(mouse_strain~timepoint, ncol=2) + 
        scale_color_manual(values = color_celltype_map) + 
        theme_classic() + NoAxes()

###Dotplot showing marker genes for labeled clusters
celltype_dot_genes <- c("Tmem119", "P2ry12", "Aif1", "Ctss", "C1qa", "C1qb", "Hexb",
                        "Itgam", "Ly6c2","Cd14", "Lgals3", "Tgfbi", "Trem1",
                        "Trbc2", "Cd3e", "Cd8a", "Ccl5", "Cd2", "Gzma", "Gzmb",
                        "Nkg7", "Ncr1",  "Klre1", "Klrk1", 
                        "Lyz2", "H2-Aa", "Cd74", "Ifitm1", "Itgax", 
                        "Ccr7", "Tmem123", "Serpinb6b", "Relb", "Socs2",
                        "Ly6g", "S100a8", "S100a9", "Retnlg", "Ngp",
                        "Cd19", "Cd20", "Iglc1", "Iglc2", "Cd79b",
                         "Cd163", "Mrc1", "Pf4", "Lyve1", "F13a1")

seurat.sub@meta.data$named_clusters <- factor(seurat.sub@meta.data$named_clusters,
                                                         levels = c("pv.Macrophages","B.cells",
                                                                    "Neutrophils","mDCs",
                                                                    "Macrophages.DCs", "NK.cells",
                                                                    "T.cells","Monocytes","Microglia"))
DotPlot(seurat.sub, features = celltype_dot_genes,
        group.by="named_clusters", scale = TRUE, cols = c("gray80","black")) + RotatedAxis()



# stacked barplot of all cell types
seurat.sub@meta.data$timepoint <- factor(seurat.sub@meta.data$timepoint, levels = c("D4","D10"))
ggplot(seurat.sub@meta.data, aes(x=orig.ident, fill=named_clusters)) + 
        geom_bar(position="fill", color="black") + 
        facet_wrap(timepoint~ mouse_strain, scales="free_x", ncol=4) + 
        scale_fill_manual(values = color_celltype_map) + theme_bw() + RotatedAxis()

# add a column to split myeloid and lymphoid cells
seurat.sub@meta.data$cell_type_category <- ifelse(
        seurat.sub@meta.data$named_clusters %in%
                c("Macrophages.DCs",
                  "mDCs","Microglia","pv.Macrophages",
                  "Neutrophils", "Monocytes"), "Myeloid", 
        ifelse(seurat.sub@meta.data$named_clusters %in%
                       c("B.cells","NK.cells","T.cells"), "Lymphoid", 
               "Neither"))

# look at plot splitting by lymphoid/myeloid info
ggplot(seurat.sub@meta.data,
       aes(x=orig.ident, fill=named_clusters)) + 
        geom_bar(position="fill", color="black") + 
        facet_wrap(cell_type_category~timepoint_strain, scales="free_x", ncol=4) + 
        scale_fill_manual(values = color_celltype_map) + theme_classic() + RotatedAxis()


################### Full UMAP cluster markers ################### 
#################################################################

# every cell type strain comparison
seurat.sub$celltype_strain <- paste0(seurat.sub$named_clusters, "_",
                                     seurat.sub$mouse_strain)
Idents(seurat.sub) <- "celltype_strain"
seurat.sub.celltype_strain_markers <- FindAllMarkers(seurat.sub, only.pos = T)

View(seurat.sub.celltype_strain_markers %>% group_by(cluster) %>% top_n(5, avg_log2FC))

# mouse timepoint
seurat.sub$timepoint_strain <- paste0(seurat.sub$timepoint, "_",
                                      seurat.sub$mouse_strain)
Idents(seurat.sub) <- "timepoint_strain"
seurat.sub.timepoint_strain_markers <- FindAllMarkers(seurat.sub, only.pos = T)
seurat.sub.timepoint_strain_markers_D10 <- FindAllMarkers(subset(seurat.sub, timepoint == "D10"),
                                                          only.pos = T)

###########################################################
############### Microglia subclustering ###################
###########################################################

seurat.sub.microglia <- subset(seurat.sub, named_clusters %in% c("Microglia"))

seurat.sub.microglia <- NormalizeData(seurat.sub.microglia)
seurat.sub.microglia <- FindVariableFeatures(seurat.sub.microglia)
seurat.sub.microglia <- ScaleData(seurat.sub.microglia)
seurat.sub.microglia <- RunPCA(seurat.sub.microglia)
ElbowPlot(seurat.sub.microglia)
seurat.sub.microglia <- FindNeighbors(seurat.sub.microglia, dims=1:14)
seurat.sub.microglia <- RunUMAP(seurat.sub.microglia, dims=1:14)
seurat.sub.microglia <- FindClusters(seurat.sub.microglia, res=0.3)

# identify markers to remove non-microglia cells from proliferating cluster
# cluster 4 is neutrophils and monocytes, cluster 5 is low quality
seurat.sub.microglia.cluster_markers <- FindAllMarkers(seurat.sub.microglia, only.pos = T)
DimPlot(seurat.sub.microglia, label = TRUE, repel = TRUE, pt.size = 0.4) & NoAxes()

# subcluster again
seurat.sub.microglia <- subset(seurat.sub.microglia, 
                               seurat_clusters %in% c("0","1","2","3"))

seurat.sub.microglia <- NormalizeData(seurat.sub.microglia)
seurat.sub.microglia <- FindVariableFeatures(seurat.sub.microglia)
seurat.sub.microglia <- ScaleData(seurat.sub.microglia)
seurat.sub.microglia <- RunPCA(seurat.sub.microglia)
ElbowPlot(seurat.sub.microglia)
seurat.sub.microglia <- FindNeighbors(seurat.sub.microglia, dims=1:13)
seurat.sub.microglia <- RunUMAP(seurat.sub.microglia, dims=1:13)
seurat.sub.microglia <- FindClusters(seurat.sub.microglia, res=0.8)
seurat.sub.microglia <- RunTSNE(seurat.sub.microglia, dims=1:13)

seurat.sub.microglia.cluster_markers <- FindAllMarkers(seurat.sub.microglia, only.pos = T)

DimPlot(seurat.sub.microglia, label = TRUE, repel = TRUE, pt.size = 0.4) & NoAxes()
FeaturePlot(seurat.sub.microglia, features = "Tmem119")

# subcluster again, removing Tmem119 low, Lgals3 very high cluster
# likely macrophages
seurat.sub.microglia <- subset(seurat.sub.microglia, 
                               seurat_clusters != "12")

seurat.sub.microglia <- NormalizeData(seurat.sub.microglia)
seurat.sub.microglia <- FindVariableFeatures(seurat.sub.microglia)
seurat.sub.microglia <- ScaleData(seurat.sub.microglia)
seurat.sub.microglia <- RunPCA(seurat.sub.microglia)
ElbowPlot(seurat.sub.microglia)
seurat.sub.microglia <- FindNeighbors(seurat.sub.microglia, dims=1:15)
seurat.sub.microglia <- RunTSNE(seurat.sub.microglia, dims=1:15)
seurat.sub.microglia <- FindClusters(seurat.sub.microglia, res=0.6)

seurat.sub.microglia.cluster_markers <- FindAllMarkers(seurat.sub.microglia, only.pos = T)

top10_microglia<- seurat.sub.microglia.cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat.sub.microglia, features = top10_microglia$gene, group.by = "seurat_clusters") + NoLegend()

DimPlot(seurat.sub.microglia, label = TRUE, repel = TRUE, pt.size = 0.4) & NoAxes()

# order samples by strain & timepoint
seurat.sub.microglia@meta.data$mouse_id <- factor(
        seurat.sub.microglia@meta.data$mouse_id, 
        levels = c("tag487","tag484","tag481","tag480","tag483","tag486",
                   "tag494","tag493","tag489", "tag492","tag490","tag495"))
VlnPlot(seurat.sub.microglia, "Cd74", group.by="mouse_id", split.by = "mouse_strain", pt.size = 0.01)
VlnPlot(seurat.sub.microglia, "Bst2", group.by="mouse_id", split.by = "mouse_strain", pt.size = 0.01)
DotPlot(seurat.sub.microglia, features=c("Ccl3","Ccl4","Bst2","Isg15","Cd74","H2-DMa"), 
        group.by="mouse_id", scale = TRUE, cols=c("gray80","black")) 

Idents(seurat.sub.microglia) <- "mouse_id"
main_marker_avg_expr <- AverageExpression(seurat.sub.microglia, 
                                          features = c("Ccl3","Ccl4","Cd63",
                                                       "Bst2","Isg15","Ifitm3",
                                                       "Cd74","H2-DMa","H2-Ab1"))

anno_col = unique(seurat.sub.microglia@meta.data[,c("timepoint","mouse_id","mouse_strain")])
rownames(anno_col) <- anno_col$mouse_id
anno_col <- anno_col[,-2]

library(RColorBrewer)
library(pheatmap)
annot_colors=list(
        timepoint = c(D4 = "#9CE55C", D10 = "#F88379"),
        mouse_strain = c(BL6 = "#90FFFF", RAG = "#CBA4D4")
)

pheatmap(main_marker_avg_expr$RNA, scale="row", cluster_cols = F,  cluster_rows = F,
         annotation_col = anno_col, 
         color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(100), 
         annotation_colors = annot_colors, fontsize = 18,
         gaps_col = c(6), gaps_row = c(3,6,9))


##Downsample to highlight small groups##
seurat.sub.microglia.small <- subset(seurat.sub.microglia, downsample = 1000)

####Recluster downsampled cells###
seurat.sub.microglia.small <- NormalizeData(seurat.sub.microglia.small)
seurat.sub.microglia.small <- FindVariableFeatures(seurat.sub.microglia.small)
seurat.sub.microglia.small <- ScaleData(seurat.sub.microglia.small)
seurat.sub.microglia.small <- RunPCA(seurat.sub.microglia.small)
ElbowPlot(seurat.sub.microglia.small)
seurat.sub.microglia.small <- FindNeighbors(seurat.sub.microglia.small, dims=1:17)
seurat.sub.microglia.small <- RunUMAP(seurat.sub.microglia.small, dims=1:17)
seurat.sub.microglia.small <- FindClusters(seurat.sub.microglia.small, res=1.8)
seurat.sub.microglia.small <- RunTSNE(seurat.sub.microglia.small, dims=1:17)


seurat.sub.microglia.small.cluster_markers <- FindAllMarkers(seurat.sub.microglia.small, only.pos = T)
DimPlot(seurat.sub.microglia.small, reduction = "umap", pt.size = 1.0, label = TRUE, repel = TRUE)

top10.small <- seurat.sub.microglia.small.cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat.sub.microglia.small, features = top10.small$gene) + NoLegend()

View(top10.small)
seurat.sub.microglia.small <- subset(seurat.sub.microglia.small, 
                               seurat_clusters != "21")

seurat.sub.microglia.small <- NormalizeData(seurat.sub.microglia.small)
seurat.sub.microglia.small <- FindVariableFeatures(seurat.sub.microglia.small)
seurat.sub.microglia.small <- ScaleData(seurat.sub.microglia.small)
seurat.sub.microglia.small <- RunPCA(seurat.sub.microglia.small)
ElbowPlot(seurat.sub.microglia.small)
seurat.sub.microglia.small <- FindNeighbors(seurat.sub.microglia.small, dims=1:17)
seurat.sub.microglia.small <- RunUMAP(seurat.sub.microglia.small, dims=1:17)
seurat.sub.microglia.small <- FindClusters(seurat.sub.microglia.small, res=1.8)
seurat.sub.microglia.small <- RunTSNE(seurat.sub.microglia.small, dims=1:17)

seurat.sub.microglia.small.cluster_markers <- FindAllMarkers(seurat.sub.microglia.small, only.pos = T)
DimPlot(seurat.sub.microglia.small, reduction = "umap", pt.size = 1.0, label = TRUE, repel = TRUE)

top10.small <- seurat.sub.microglia.small.cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat.sub.microglia.small, features = top10.small$gene) + NoLegend()

seurat.sub.microglia.small$timepoint_strain <- factor(x = seurat.sub.microglia.small$timepoint_strain, 
                                                      levels = c("D4_BL6", "D10_BL6", "D4_RAG", "D10_RAG"))

FeaturePlot(seurat.sub.microglia.small, features = c("Tmem119", "P2ry12", "Itgam"), ncol = 3, pt.size = 0.8)

FeaturePlot(seurat.sub.microglia.small, features = c("Bst2", "Isg15", "Ifitm3", "Cd74", 
                                                     "H2-DMa", "H2-Ab1"), order = TRUE, pt.size = 0.8)
FeaturePlot(seurat.sub.microglia.small, features = c("Cd40", "Cd70", "Ifngr1", 
                                                     "Tnfrsf1a", "Tnfrsf1b", "Tnfsf4", 
                                                     "Tnfsf9", "Tnfsf14", "Tnfsf18"), order = TRUE, pt.size = 0.8)

FeaturePlot(seurat.sub.microglia.small, features= "Bst2", order = TRUE, pt.size = 1.2, ncol = 2, split.by = "timepoint_strain")
FeaturePlot(seurat.sub.microglia.small, features= "Cd74", order = TRUE, pt.size = 1.2, ncol = 2, split.by = "timepoint_strain")
FeaturePlot(seurat.sub.microglia.small, features= "H2-DMa", order = TRUE, pt.size = 1.2, ncol = 2, split.by = "timepoint_strain")

VlnPlot(seurat.sub.microglia.small, "Cd74", group.by="timepoint_strain", split.by = "mouse_strain", pt.size = 0.4)
VlnPlot(seurat.sub.microglia.small, "Bst2", group.by="timepoint_strain", split.by = "mouse_strain", pt.size = 0.4)
VlnPlot(seurat.sub.microglia.small, "H2-Ab1", group.by="timepoint_strain", split.by = "mouse_strain", pt.size = 0.4)

# annotate downsampled microglia clusters based on markers
microglia.small.cluster.ids <- 0:19
new.microglia.small.cluster.ids <- c("Homeostatic", "Homeostatic", "Homeostatic", "Secretory",
                                     "APC", "Ifn.Responsive", "APC", "Homeostatic",
                                     "Secretory", "Secretory", "Homeostatic", "Secretory",
                                     "Homeostatic", "Homeostatic", "Cycling", "APC",
                                     "Secretory", "Ifn.Responsive", "Homeostatic", "Cycling")
seurat.sub.microglia.small@meta.data$microglia_named_clusters <- plyr::mapvalues(x = seurat.sub.microglia.small@meta.data$seurat_clusters, 
                                                                                 from = microglia.small.cluster.ids, 
                                                                                 to = new.microglia.small.cluster.ids)
################### Renamed UMAP ################### 
####################################################

seurat.sub.microglia.small@meta.data$microglia_named_clusters <- factor(seurat.sub.microglia.small@meta.data$microglia_named_clusters,
                                                                        levels = c("Homeostatic","Secretory","Ifn.Responsive", "APC",
                                                                                   "Cycling"))
seurat.sub.microglia.small@meta.data$timepoint_strain <- factor(seurat.sub.microglia.small@meta.data$timepoint_strain,
                                                                levels = c("D4_BL6","D4_RAG","D10_BL6",
                                                                           "D10_RAG"))

# stacked bar
color_celltype_map <- c("Homeostatic" = "#b1b3bc",
                        "Secretory" = "#5AE6AC",
                        "Ifn.Responsive" = "#63E5FF",
                        "APC" = "#F88379",
                        "Cycling" ="#DCD0FF")

ggplot(seurat.sub.microglia.small@meta.data, aes(x=orig.ident, fill=microglia_named_clusters)) + 
        geom_bar(position="fill", width = 0.95, color="black") + 
        facet_wrap(timepoint~ mouse_strain, scales="free_x", ncol=4) + 
        scale_fill_manual(values = color_celltype_map) + 
        theme_classic(base_size = 24) + RotatedAxis()


# UMAP colored by cell type
DimPlot(seurat.sub.microglia.small, group.by = "microglia_named_clusters", cols = color_celltype_map, 
        label = TRUE, pt.size = 1.0, label.size = 10, repel = TRUE) & NoAxes()
DimPlot(seurat.sub.microglia.small, group.by = "timepoint_strain", pt.size = 1.0, 
        label = FALSE, order = TRUE) & NoAxes()

###Dotplot showing marker genes for labeled clusters
micro_dot_genes <- c("Tmem119", "P2ry12", "Olfml3", "Sall1", "Cx3cr1",  
                     "Ccl3", "Ccl4", "Csf1", "Il1b", "Aif1", "Apoe", "Nfkbiz",
                     "Bst2", "Isg15", "Ly6e", "Ifitm3", "Fcgr1", "Ccl12", "Cd74",
                     "H2-Ab1", "H2-Aa", "H2-DMa", "Ccl5", "Mki67", "Top2a", 
                     "Stmn1", "Hist1h2ap", "Hist1h2ae")
seurat.sub.microglia.small@meta.data$microglia_named_clusters <- factor(seurat.sub.microglia.small@meta.data$microglia_named_clusters,
                                                                        levels = c("Cycling", "APC", "Ifn.Responsive", "Secretory", "Homeostatic"))
DotPlot(seurat.sub.microglia.small, features = micro_dot_genes,
        group.by="microglia_named_clusters", scale = TRUE, cols = c("gray80","black")) + RotatedAxis()

#########################################################
################ NK/T cell subclustering ################
#########################################################

seurat.sub.NKT <- subset(seurat.sub, named_clusters %in% c("NK.cells","T.cells","Cycling.Lymphoid"))

seurat.sub.NKT <- NormalizeData(seurat.sub.NKT)
seurat.sub.NKT <- FindVariableFeatures(seurat.sub.NKT)
seurat.sub.NKT <- ScaleData(seurat.sub.NKT)
seurat.sub.NKT <- RunPCA(seurat.sub.NKT)
ElbowPlot(seurat.sub.NKT)
seurat.sub.NKT <- FindNeighbors(seurat.sub.NKT, dims=1:16)
seurat.sub.NKT <- RunTSNE(seurat.sub.NKT, dims=1:16)
seurat.sub.NKT <- FindClusters(seurat.sub.NKT, res=0.8)
seurat.sub.NKT <- RunUMAP(seurat.sub.NKT, dims=1:16)

seurat.sub.NKT.cluster_markers <- FindAllMarkers(seurat.sub.NKT, only.pos = T)
View(seurat.sub.NKT.cluster_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC))


Idents(seurat.sub.NKT) <- "mouse_strain"
seurat.sub.NKT.strain_markers <- FindAllMarkers(seurat.sub.NKT, only.pos = T)

# cluster 9 is monocyte doublets (Lyz2, Ly6c2, Ccr2)
# cluster 12 is microglia doublets (P2ry12, Sparc, Cd81)
seurat.sub.NKT <- subset(seurat.sub.NKT, 
                         seurat_clusters != "9" & 
                                 seurat_clusters != "12")
seurat.sub.NKT <- NormalizeData(seurat.sub.NKT)
seurat.sub.NKT <- FindVariableFeatures(seurat.sub.NKT)
seurat.sub.NKT <- ScaleData(seurat.sub.NKT)
seurat.sub.NKT <- RunPCA(seurat.sub.NKT)
ElbowPlot(seurat.sub.NKT)
seurat.sub.NKT <- FindNeighbors(seurat.sub.NKT, dims=1:15)
seurat.sub.NKT <- RunUMAP(seurat.sub.NKT, dims=1:15)
seurat.sub.NKT <- FindClusters(seurat.sub.NKT, res=0.8)

seurat.sub.NKT.cluster_markers <- FindAllMarkers(seurat.sub.NKT, only.pos = T)
View(seurat.sub.NKT.cluster_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC))

DimPlot(seurat.sub.NKT , reduction = "umap", pt.size = 1.0, label = TRUE, repel = TRUE)

top10.NKT <- seurat.sub.NKT.cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat.sub.NKT, features = top10.NKT$gene) + NoLegend()

#########################################################
################ T cell subclustering ###################
#########################################################

seurat.sub.T <- subset(seurat.sub.NKT, 
                       seurat_clusters != "0" & 
                               seurat_clusters != "1" & 
                               seurat_clusters != "2" & 
                               seurat_clusters != "5" &
                               seurat_clusters != "12")

seurat.sub.T <- NormalizeData(seurat.sub.T)
seurat.sub.T <- FindVariableFeatures(seurat.sub.T)
seurat.sub.T <- ScaleData(seurat.sub.T)
seurat.sub.T <- RunPCA(seurat.sub.T)
ElbowPlot(seurat.sub.T)
seurat.sub.T <- FindNeighbors(seurat.sub.T, dims=1:15)
seurat.sub.T <- RunUMAP(seurat.sub.T, dims=1:15)
seurat.sub.T <- FindClusters(seurat.sub.T, res=0.6)

seurat.sub.T.cluster_markers <- FindAllMarkers(seurat.sub.T, only.pos = T)
View(seurat.sub.T.cluster_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC))

# remove NK cluster & remaining, misclassified RAG cells
seurat.sub.T <- subset(seurat.sub.T, 
                       seurat_clusters != "8" & 
                        seurat_clusters != "9" &
                       mouse_strain != "RAG")
seurat.sub.T <- NormalizeData(seurat.sub.T)
seurat.sub.T <- FindVariableFeatures(seurat.sub.T)
seurat.sub.T <- ScaleData(seurat.sub.T)
seurat.sub.T <- RunPCA(seurat.sub.T)
ElbowPlot(seurat.sub.T)
seurat.sub.T <- FindNeighbors(seurat.sub.T, dims=1:13)
seurat.sub.T <- RunUMAP(seurat.sub.T, dims=1:13)
seurat.sub.T <- FindClusters(seurat.sub.T, res=1.8)

seurat.sub.T.cluster_markers <- FindAllMarkers(seurat.sub.T, only.pos = T)
View(seurat.sub.T.cluster_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC))

top10.T <- seurat.sub.T.cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(seurat.sub.T, features = top10.T$gene) + NoLegend()

DimPlot(seurat.sub.T , reduction = "umap", pt.size = 1.0, label = TRUE, repel = TRUE)

seurat.sub.T <- FlipReductionCoordinates(seurat.sub.T)
DimPlot(seurat.sub.T , reduction = "umap", pt.size = 1.0, label = TRUE, repel = TRUE)

# annotate clusters based on markers
current.cluster.ids.T <- 0:14
new.cluster.ids.T <- c("Memory","Cycling", "Tregs", "Tregs",
                       "CD4", "Naive", "CD8", "CD8", "CD8.effector",
                       "CD8.effector", "CD8.effector","Naive","CD4.effector",
                       "Late.stage.effector","Gamma.delta")
seurat.sub.T@meta.data$granular_named_clusters <- plyr::mapvalues(x = seurat.sub.T@meta.data$seurat_clusters, 
                                                                  from = current.cluster.ids.T, 
                                                                  to = new.cluster.ids.T)

seurat.sub.T@meta.data$granular_named_clusters <- factor(seurat.sub.T@meta.data$granular_named_clusters,
                                                         levels = c("Naive","Memory", "CD4", "CD4.effector",
                                                                    "CD8", "CD8.effector",
                                                                    "Late.stage.effector",
                                                                    "Tregs", "Gamma.delta",
                                                                    "Cycling"))

color_celltype_map <- c("Naive" = "gray70", 
                        "CD8.effector" = "#FF61CC",
                        "CD8" = "#F8766D",
                        "CD4" = "#CAB2D6",
                        "CD4.effector" = "#6A3D9A",
                        "Late.stage.effector" = "#FF0090",
                        "Tregs" ="#00A9FF",
                        "Gamma.delta" = "forestgreen",
                        "Memory" = "#FFD1AD",
                        "Cycling" = "#00C19a")

# UMAP colored by cell type
DimPlot(seurat.sub.T, pt.size = 1.8, label.size = 10, group.by = "granular_named_clusters", label = TRUE, repel = TRUE) + 
        scale_color_manual(values = color_celltype_map) & NoAxes() & NoLegend()

color_timepoint_map <- c("D4_BL6" = "#85BB61", "D10_BL6" = "#E68613")

# UMAP colored by Timepoint
DimPlot(seurat.sub.T, pt.size = 1.8, group.by = "timepoint_strain", label = FALSE, order = c("D4_BL6", "D10_BL6")) + 
        scale_color_manual(values = color_timepoint_map)  & NoAxes()

# stacked bar
ggplot(seurat.sub.T@meta.data, aes(x=orig.ident, fill=granular_named_clusters)) + 
        geom_bar(position="fill", color="black") + 
        facet_wrap(timepoint~ mouse_strain, scales="free_x", ncol=4) + 
        scale_fill_manual(values = color_celltype_map) + 
        theme_classic(base_size = 24) + RotatedAxis()

###Dotplot showing marker genes for labeled clusters
T_dot_genes <- c("Sell", "Ccr7", "Cd127", "Il7r", "Il2rg", "Cd44", "Cd7", "Klra1", "Klra6", "Klra7", "Fcer1g",
                 "Lef1", "Tcf7", "Klf2", "Satb1", "Cd8a", "Ccl5", "Cxcr6", "Gzma", "Gzmb", "Gzmk", "Ifng",
                 "Tnf", "Pdcd1", "Havcr2", "Tnfrsf9", "Ctla4", "Il2ra", "Foxp3", "Tcrg-C1", "Trdc",
                 "Trdv4", "Il17a", "Mki67", "Top2a", "Hist1h1b", "Hist1h2ae", "Hist1h4d")

seurat.sub.T@meta.data$granular_named_clusters <- factor(seurat.sub.T@meta.data$granular_named_clusters,
                                                                        levels = c("Cycling", "Gamma.delta", "Tregs",
                                                                                   "Late.stage.effector",
                                                                                   "CD8.effector", "CD8",
                                                                                   "CD4.effector", "CD4",
                                                                                   "Memory", "Naive"))
DotPlot(seurat.sub.T, features = T_dot_genes,
        group.by="granular_named_clusters", scale = TRUE, cols = c("gray80","black")) + RotatedAxis()


