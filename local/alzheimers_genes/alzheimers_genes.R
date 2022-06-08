#' ---
#' title: "R Notebook"
#' output:
#'   html_document:
#'     df_print: paged
#' ---
#' 
## ---- warning=FALSE-----------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

theme_set(
  theme(
    axis.text.x.bottom = element_text(angle = 90),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid = element_line(color = "lightgray")
  )
)

#' 
#' # Loading data
## -----------------------------------------------------------------------------
alzheimers <- readRDS("/Users/steinacker/Masters_thesis_work/data/alzheimers/alzheimers.rds")
outdir <- "/Users/steinacker/Masters_thesis_work/code_final/alzheimers_genes/output/"
DefaultAssay(alzheimers) <- "RNA"

#' 
#' # Quality control
#' 
#' ## Calculating mitochondrial RNA content
## -----------------------------------------------------------------------------
alzheimers <- PercentageFeatureSet(alzheimers, "^MT-", col.name = "percent_mito")

#' 
#' ## Violin plots for QC features
## -----------------------------------------------------------------------------
feats <- c("nFeature_RNA", "nCount_RNA")

VlnPlot(alzheimers, group.by = "orig.ident", features = "nFeature_RNA", ncol = 1, pt.size = 0) +
  ggtitle(label = "Number of detected genes per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_nfeat_alz.png"), width = 5, height = 5, dpi = 600)

VlnPlot(alzheimers, group.by = "orig.ident", features = "nCount_RNA", ncol = 1, pt.size = 0) +
  ggtitle(label = "Number of RNA transcripts per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_ncount_alz.png"), width = 5, height = 5, dpi = 600)

VlnPlot(alzheimers, group.by = "orig.ident", features = "percent_mito", ncol = 1, pt.size = 0) +
  ggtitle(label = "Percent mitochondrial reads per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_mito_alz.png"), width = 5, height = 5, dpi = 600)

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
neuronal <- c("MAP2", "DCX")
excitatory <- c("RBFOX3", "GRIN1", "HS3ST2", "CAMK2A")
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
DotPlot(alzheimers, features = unique(all_markers)) + 
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

FeaturePlot(alzheimers, features = rep_markers, ncol = 4)
ggsave(filename = paste0(outdir, "biomarkers_alz.png"), width = 20, height = 10, dpi = 600)

cell_types_list <- c(
  "0" = "Oligodendrocytes",
  "1" = "Excit_neurons",
  "2" = "Astrocytes",
  "3" = "OPCs",
  "4" = "Inhib_neurons",
  "5" = "Inhib_neurons",
  "6" = "Excit_neurons",
  "7" = "Excit_neurons",
  "8" = "Microglia",
  "9" = "Excit_neurons",
  "10" = "Excit_neurons",
  "11" = "Excit_neurons",
  "12" = "Endothelial_cells",
  "13" = "Excit_neurons"
)

alzheimers <- AddCellTypes(alzheimers, cell_types_list)

DimPlot(alzheimers, group.by = "seurat_clusters", label = T) +
  ggtitle("Cluster") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(outdir, "clusters_alz_umap_all.png"), width = 4, height = 4, dpi = 500)

DimPlot(alzheimers, group.by = "cell_type", label = T, label.size = 3, label.box = T, repel = T) +
  ggtitle("Cell type") +
  theme(legend.position = "none")
ggsave(filename = paste0(outdir, "cell_type_alz_umap_all.png"), width = 4, height = 4, dpi = 500)

#' 
#' # Adding diagnosis metadata
## -----------------------------------------------------------------------------
alzheimers <- AddMetaData(alzheimers, as.factor(substr(alzheimers$orig.ident, 1, 2)), "diagnosis")

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
  "FTL", # Kenkhius 2021
  "SPP1", # Ochocka 2021
  "APOE", # Ochocka 2021
  "CD74", # Ochocka 2021
  "FCGR3A", 
  "CST3", # Keren-Shaul 2017
  "CSF1R", # Keren-Shaul 2017
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

cell_data <- alzheimers %>%
  FetchData(c(fetched_metadata, fetched_data_gene))

#' 
#' # Activated microglia biomarkers
## -----------------------------------------------------------------------------
for(i in fetched_data_gene_microglia) {
  ggplot(data = filter(cell_data, cell_type == "Microglia"), mapping = aes_string(x = "diagnosis", y = i)) +
    geom_violin(fill = "gray") +
    stat_summary(geom = "point", fun = "mean", shape = 4, color = "red") +
    ggtitle(i) +
    stat_compare_means(label = "p.signif", vjust = 1, hjust = -1) +
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
alzheimers <- SetIdent(alzheimers, value="cell_type")
celltypes <- as.character(unique(Idents(alzheimers)))
deas_celltype <- sapply(as.character(celltypes), function(x) NULL)
for(celltype in as.character(celltypes)) deas_celltype[[celltype]] <- FindMarkers(alzheimers, group.by = "diagnosis", ident.1 = "AD", subset.ident = celltype, logfc.threshold = 0.05)
for(celltype in as.character(celltypes)) deas_celltype[[celltype]]$gene <- rownames(deas_celltype[[celltype]])
for(celltype in as.character(celltypes)) rownames(deas_celltype[[celltype]]) <- NULL
# for(celltype in as.character(celltypes)) deas_celltype[[celltype]] <- deas_celltype[[celltype]][which(deas_celltype[[celltype]]$p_val_adj < 0.01),]

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
               pAdjustMethod = "BH")
  return(gse)
}

gses_celltype <- sapply(as.character(celltypes), function(x) NULL)
for(celltype in as.character(celltypes)) rownames(deas_celltype[[celltype]]) <- deas_celltype[[celltype]]$gene
for(celltype in as.character(celltypes)) gses_celltype[[celltype]] <- gse_dotplot(deas_celltype[[celltype]])
gses_celltype_alz <- gses_celltype

#' 
## -----------------------------------------------------------------------------
dotplot(gses_celltype[["Microglia"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(#axis.text.y = element_text(size = 9),
        strip.text.x = element_text(size = 12)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
  ggtitle("Microglia")
ggsave(filename = paste0(outdir, "microglia_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Astrocytes"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12)) +
  ggtitle("Astrocytes")
ggsave(filename = paste0(outdir, "astrocytes_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Endothelial_cells"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(#axis.text.y = element_text(size = 9),
        strip.text.x = element_text(size = 12)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
#  theme(axis.text.y = element_text(size = 8)) +
  ggtitle("Endothelial cells")
ggsave(filename = paste0(outdir, "endothelial_GO.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["OPCs"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12)) +
  ggtitle("OPCs")
ggsave(filename = paste0(outdir, "opcs_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Oligodendrocytes"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        strip.text.x = element_text(size = 12)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  ggtitle("Oligodendrocytes")
ggsave(filename = paste0(outdir, "oligodendrocytes_go.png"), width = 7, height = 5, dpi = 600)
#Neither of the lower two cell types showed any upregulated or downregulated GO terms
# dotplot(gses_celltype[["Inhib_neurons"]], showCategory=5, split=".sign") +
#   facet_grid(.~.sign) +
#   ggtitle("Inhibitory neurons")
# 
# dotplot(gses_celltype[["Excit_neurons"]], showCategory=5, split=".sign") +
#   facet_grid(.~.sign) +
#   ggtitle("Excitatory neurons")

