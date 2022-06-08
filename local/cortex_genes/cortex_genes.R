#' ---
#' title: "R Notebook"
#' output: html_notebook
#' ---
#' 
## -----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggpubr)
library(org.Hs.eg.db)
library(clusterProfiler)

parkinsons_ctx <- readRDS("/Users/steinacker/Masters_thesis_work/data/parkinsons/cortex/parkinsons_ctx.rds")
DefaultAssay(parkinsons_ctx) <- "RNA"
outdir <- "/Users/steinacker/Masters_thesis_work/code_final/cortex_genes/output/"
if(!dir.exists(outdir)) dir.create(outdir)

#' 
#' ## Calculating mitochondrial RNA content
## -----------------------------------------------------------------------------
parkinsons_ctx <- PercentageFeatureSet(parkinsons_ctx, "^MT-", col.name = "percent_mito")

#' 
#' ## Violin plots for QC features
## -----------------------------------------------------------------------------
feats <- c("nFeature_RNA", "nCount_RNA")

VlnPlot(parkinsons_ctx, group.by = "orig.ident", features = "nFeature_RNA", ncol = 1, pt.size = 0) +
  ggtitle(label = "Number of detected genes per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_nfeat_cortex.png"), width = 5, height = 5, dpi = 600)

VlnPlot(parkinsons_ctx, group.by = "orig.ident", features = "nCount_RNA", ncol = 1, pt.size = 0) +
  ggtitle(label = "Number of RNA transcripts per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_ncount_cortex.png"), width = 5, height = 5, dpi = 600)

VlnPlot(parkinsons_ctx, group.by = "orig.ident", features = "percent_mito", ncol = 1, pt.size = 0) +
  ggtitle(label = "Percent mitochondrial reads per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_mito_sn.png"), width = 5, height = 5, dpi = 600)

#' 
#' # Writing function to add cell types to samples
## -----------------------------------------------------------------------------
AddCellTypes <- function(seurat_obj, cell_types_list) {
  cell_types_df <- seurat_obj %>% 
    FetchData("seurat_clusters") %>% 
    mutate(cell_type = as.factor(cell_types_list[seurat_clusters])) %>%
    dplyr::select(-seurat_clusters)
  seurat_obj <- AddMetaData(seurat_obj, cell_types_df, col.name = "cell_type")
  return(seurat_obj)
}

#' 
#' # Assigning cell types
## -----------------------------------------------------------------------------
neuronal <- c("MAP2", "DCX", "RBFOX3")
excitatory <- c("GRIN1", "HS3ST2", "CAMK2A")
inhibitory <- c("GAD1", "GAD2", "CALB2", "CNR1")
astrocyte <- c("GFAP", "AQP4", "GJA1", "SLC1A3")
oligodendrocyte <- c("PLP1", "MOG", "MBP")
opc <- c("COL9A1", "VCAN", "PDGFRA")
microglia <- c("P2RY12", "FYB1")
endothelial <- c("FLT1", "PECAM1")
all_markers <- c(
  microglia,
  excitatory,
  inhibitory,
  astrocyte,
  oligodendrocyte,
  opc,
  neuronal,
  endothelial
)
DotPlot(parkinsons_ctx, features = unique(all_markers)) + 
  theme(axis.text.x = element_text(angle = 90))

rep_markers <- c(
  "CAMK2A",
  "GAD1",
  "GFAP",
  "PLP1",
  "VCAN",
  "FYB1",
  "FLT1"
)

FeaturePlot(parkinsons_ctx, features = rep_markers, ncol = 4)
ggsave(filename = paste0(outdir, "biomarkers_pd_ctx.png"), width = 20, height = 10, dpi = 600)

cell_types_list <- c(
  "0" = "Oligodendrocytes",
  "1" = "Oligodendrocytes",
  "2" = "Oligodendrocytes",
  "3" = "Inhib_neurons",
  "4" = "Excit_neurons",
  "5" = "Excit_neurons",
  "6" = "Excit_neurons",
  "7" = "Excit_neurons",
  "8" = "Inhib_neurons",
  "9" = "Astrocytes",
  "10" = "Astrocytes",
  "11" = "Astrocytes",
  "12" = "Excit_neurons",
  "13" = "Excit_neurons",
  "14" = "Astrocytes",
  "15" = "OPCs",
  "16" = "Microglia",
  "17" = "OPCs",
  "18" = "Inhib_neurons",
  "19" = "Inhib_neurons",
  "20" = "Excit_neurons",
  "21" = "Excit_neurons",
  "22" = "Inhib_neurons",
  "23" = "Inhib_neurons",
  "24" = "Endothelial_cells"
)

parkinsons_ctx <- AddCellTypes(parkinsons_ctx, cell_types_list)
DimPlot(parkinsons_ctx, label = T) +
  ggtitle("Clusters") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(outdir, "clusters_pd_ctx_umap_all.png"), width = 4, height = 4, dpi = 500)

DimPlot(parkinsons_ctx, group.by = "cell_type", label = T, label.size = 3, label.box = T, repel = T) +
  ggtitle("Cell type") +
  theme(legend.position = "none")
ggsave(filename = paste0(outdir, "cell_types_pd_ctx.png"), width = 4, height = 4, dpi = 600)

#' 
#' # Adding diagnosis metadata
## -----------------------------------------------------------------------------
parkinsons_ctx <- AddMetaData(parkinsons_ctx, as.factor(str_split_fixed(parkinsons_ctx$orig.ident, "_", 3)[,2]), "diagnosis")

#' 
#' # Gene expression analysis
## -----------------------------------------------------------------------------
fetched_metadata <- c(
  "seurat_clusters",
  "diagnosis",
  "cell_type",
  "UMAP_1",
  "UMAP_2",
  "orig.ident"
)
fetched_data_gene_microglia <- c(
  "FTL",
  "SPP1",#
  "APOE",#
  "CD74",#
  "FCGR3A",
  "CST3",
  "CSF1R",
  "PTPRC"
)
fetched_data_gene_astrocytes <- c(
  "CHI3L1",
  "C3",
  "S100B",
  "CRYAB",
  "MAOB",
  "NFAT5",
  "HSPB1",
  "MT2A"
)

fetched_data_gene <- c(
  fetched_data_gene_microglia,
  fetched_data_gene_astrocytes
  ) %>% 
  unique()

cell_data <- parkinsons_ctx %>%
  FetchData(c(fetched_metadata, fetched_data_gene)) %>%
  mutate(diagnosis = ifelse(diagnosis == "Ctl", "NC", "PD")) %>%
  mutate(diagnosis = factor(diagnosis, levels = c("PD", "NC")))

#' 
#' # Activated microglia biomarkers
## -----------------------------------------------------------------------------
for(i in fetched_data_gene_microglia) {
  ggplot(data = filter(cell_data, cell_type == "Microglia"), mapping = aes_string(x = "diagnosis", y = i)) +
    geom_violin(fill = "gray") +
    stat_summary(geom = "point", fun = "mean", shape = 4, color = "red") +
    ggtitle(i) +
    stat_compare_means(label = "p.signif", vjust = 1, hjust = -2) +
    theme(axis.title = element_blank(),
          axis.text.x.bottom = element_text(angle = 0, size = 12),
          axis.text.y.left = element_text(size = 12))
  ggsave(filename = paste0(outdir, i, "_violin_microglia.png"), width = 2, height = 2)
}

#' 
#' # Reactive astrocyte biomarkers
## -----------------------------------------------------------------------------
for(i in fetched_data_gene_astrocytes) {
  ggplot(data = filter(cell_data, cell_type == "Astrocytes"), mapping = aes_string(x = "diagnosis", y = i)) +
    geom_violin(fill = "gray") +
    stat_summary(geom = "point", fun = "mean", shape = 4, color = "red") +
    ggtitle(i) +
    stat_compare_means(label = "p.signif", vjust = 1, hjust = -1) +
    theme(axis.title = element_blank(),
          axis.text.x.bottom = element_text(angle = 0, size = 12),
          axis.text.y.left = element_text(size = 12))
  ggsave(filename = paste0(outdir, i, "_violin_astrocytes.png"), width = 2, height = 2)
}

#' 
#' # Gene ontology analysis
## -----------------------------------------------------------------------------
addmargins(table(cell_data$diagnosis, cell_data$cell_type))

#' 
## -----------------------------------------------------------------------------
set.seed(2468)
parkinsons_ctx <- SetIdent(parkinsons_ctx, value="cell_type")
celltypes <- as.character(unique(Idents(parkinsons_ctx)))
deas_celltype <- sapply(as.character(celltypes), function(x) NULL)

for(celltype in as.character(celltypes)) deas_celltype[[celltype]] <- FindMarkers(parkinsons_ctx, group.by = "diagnosis", ident.1 = "PD", subset.ident = celltype, logfc.threshold = 0.05)

for(celltype in as.character(celltypes)) deas_celltype[[celltype]]$gene <- rownames(deas_celltype[[celltype]])

for(celltype in as.character(celltypes)) rownames(deas_celltype[[celltype]]) <- NULL

#for(celltype in as.character(celltypes)) deas_celltype[[celltype]] <- deas_celltype[[celltype]][which(deas_celltype[[celltype]]$p_val_adj < 0.01),]

gse_dotplot <- function(dea){
  genelist <- bitr(rownames(dea), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  genelist <- merge(genelist, dea[,"avg_log2FC", drop=F], by.x="SYMBOL", by.y="row.names")
  genelist <- genelist[order(genelist$avg_log2FC, decreasing = T),]

  genelist_FC <- genelist$avg_log2FC
  names(genelist_FC) <- genelist$ENTREZID
  gse <- gseGO(geneList=genelist_FC,
               ont ="ALL",
               keyType = "ENTREZID",
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               seed = T,
               nPermSimple = 100000,
               verbose = TRUE,
               OrgDb = org.Hs.eg.db,
               pAdjustMethod = "BH",
               eps = 0)
  return(gse)
}

#tmp <-bitr(c("3824","498","29990","23457","7415"), fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
	


gses_celltype <- sapply(as.character(celltypes), function(x) NULL)

for(celltype in as.character(celltypes)) rownames(deas_celltype[[celltype]]) <- deas_celltype[[celltype]]$gene

for(celltype in as.character(celltypes)) gses_celltype[[celltype]] <- gse_dotplot(deas_celltype[[celltype]])

gses_celltype_ctx <- gses_celltype

#' 
## -----------------------------------------------------------------------------
dotplot(gses_celltype[["Microglia"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 18)) +
  ggtitle("Microglia")
ggsave(filename = paste0(outdir, "microglia_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Astrocytes"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 18)) +
  ggtitle("Astrocytes")
ggsave(filename = paste0(outdir, "astrocytes_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Oligodendrocytes"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(#axis.text.y = element_text(size = 9),
        strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 18)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
  ggtitle("Oligodendrocytes")
ggsave(filename = paste0(outdir, "oligodendrocytes_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["OPCs"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 18)) +
  ggtitle("OPCs")
ggsave(filename = paste0(outdir, "opcs_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Endothelial_cells"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 18)) +
  ggtitle("Endothelial cells")
ggsave(filename = paste0(outdir, "endothelial_cells_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Excit_neurons"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 18)) +
  ggtitle("Excitatory neurons")
ggsave(filename = paste0(outdir, "excitatory_neurons_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Inhib_neurons"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12),
        plot.title = element_text(size = 18)) +
  ggtitle("Inhibitory neurons")
ggsave(filename = paste0(outdir, "inhibitory_neurons_go.png"), width = 7, height = 5, dpi = 600)

