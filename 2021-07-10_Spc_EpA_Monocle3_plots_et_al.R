load('/Users/ellorens/Downloads/Troy_UI_Inj_after_analysis_pipe_TOM_regTom_ENRIC_VER_2021_07_10.Robj')
load('/Users/ellorens/Downloads/Troy_UI_Inj_spc_all_cells_2105_ENRIC_VERSION_save_21-07-10 (1).Robj')
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(dplyr)

# Monocle3 on epc.combined 210710
cds_ep <- as.cell_data_set(epc.combined)
cds_ep <- cluster_cells(cds_ep, resolution = 0.0025)
cds_ep <- learn_graph(cds_ep, learn_graph_control=list(ncenter=60),  close_loop = FALSE, use_partition = TRUE)
plot_cells(cds_ep, color_cells_by="cluster", label_groups_by_cluster = FALSE, label_leaves = TRUE, label_branch_points = TRUE, cell_size = 1)
plot_cells(cds_ep)

cds_ep <- order_cells(cds_ep)

plot_cells(cds_ep,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size= 0, cell_size = 1.2, trajectory_graph_segment_size = 1)





# Remove likely remaining doublets based on mis-clustering. 19 cells in total.
troy <- spc.combined
plot <- DimPlot(troy, pt.size = 0.2)
select.cells <- CellSelector(plot = plot)
Idents(troy, cells = select.cells) <- "probdb"

troy_clean  <- subset(troy, idents = c("probdb"), invert = TRUE)
DimPlot(troy_clean, split.by = 'orig.ident')
FeaturePlot(troy_clean, features = c('Sox9','Foxj1','Rarres2','tdTomatoWRPE'), cols = c("lightgray", "red"))

cds <- as.cell_data_set(troy_clean)

plot_cells(cds, show_trajectory_graph = FALSE, cell_size = 1)
cds <- cluster_cells(cds, partition_qval = 0.5, resolution = 0.0025)
cds <- learn_graph(cds, learn_graph_control=list(ncenter=75),  close_loop = FALSE, use_partition = TRUE)
plot_cells(cds, color_cells_by="partition", label_groups_by_cluster = FALSE, graph_label_size= 0,label_leaves = FALSE, label_branch_points = FALSE, cell_size = 1, trajectory_graph_segment_size = 1)
plot_cells(cds)

cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size= 0, cell_size = 1.2, trajectory_graph_segment_size = 1)

# Subset EC and Fib/peri clusters
perivascular  <- subset(troy_clean, idents = c("EC", "Fib/peri"), invert = FALSE)
perivascular <- FindVariableFeatures(perivascular, selection.method = "vst", nfeatures = 1000)
perivascular <- ScaleData(perivascular, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(perivascular))
perivascular <- RunPCA(perivascular, features = VariableFeatures(perivascular), nfeatures.print = 10)
ElbowPlot(perivascular, ndims = 30)
DimHeatmap(perivascular, dims = 1:10, cells = 100, balanced = TRUE)
perivascular <- RunUMAP(perivascular, dims = 1:7)
perivascular <- FindNeighbors(perivascular, dims = 1:7)
perivascular <- FindClusters(perivascular, resolution = 0.2)
DimPlot(perivascular, group.by = 'orig.ident')
markers_pv <- FindAllMarkers(perivascular, only.pos = TRUE, logfc.threshold = 0.25)
FeaturePlot(perivascular, c('Cldn5','Tek','Egfl7','Ly6a'), cols = c("lightgray", "red"))
FeaturePlot(perivascular, c('Dcn','Col1a2','Pdpn','Lum'), cols = c("lightgray", "red"))
FeaturePlot(perivascular, c('Fam167b','Col4a1','Plvap','Mcam'), cols = c("lightgray", "red"))
FeaturePlot(perivascular, c('Mki67','Mcm2','Aurkb','Ccna2'), cols = c("lightgray", "red"))
# Monocle3 on perivascular cells
cds_pv <- as.cell_data_set(perivascular)
cds_pv <- cluster_cells(cds_pv, partition_qval = 0.5,resolution = 0.5)
cds_pv <- learn_graph(cds_pv, learn_graph_control=list(ncenter=50),  close_loop = FALSE, use_partition = TRUE)
plot_cells(cds_pv, color_cells_by="partition", label_groups_by_cluster = FALSE, label_leaves = TRUE, label_branch_points = TRUE, cell_size = 1)
plot_cells(cds_pv)

cds_pv <- order_cells(cds_pv)

plot_cells(cds_pv,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size= 0, cell_size = 1.2, trajectory_graph_segment_size = 1)


