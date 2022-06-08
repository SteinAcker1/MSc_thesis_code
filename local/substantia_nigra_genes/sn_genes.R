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

theme_set(
  theme(
    axis.text.x.bottom = element_text(angle = 90),
    panel.background = element_rect(fill = "white", colour = "grey50"),
    panel.grid = element_line(color = "lightgray")
  )
)

parkinsons_sn <- readRDS("/Users/steinacker/Masters_thesis_work/data/parkinsons/substantia_nigra/parkinsons.rds")
DefaultAssay(parkinsons_sn) <- "RNA"
outdir <- "~/Masters_thesis_work/code_final/substantia_nigra_genes/output/"
if(!dir.exists(outdir)) dir.create(outdir)
DimPlot(parkinsons_sn, group.by = "orig.ident")
DimPlot(parkinsons_sn, group.by = "seurat_clusters", label = T)

#' 
#' # Quality control
#' ## Calculating mitochondrial RNA content
## -----------------------------------------------------------------------------
parkinsons_sn <- PercentageFeatureSet(parkinsons_sn, "^MT-", col.name = "percent_mito")

#' 
#' ## Violin plots for QC features
## -----------------------------------------------------------------------------
feats <- c("nFeature_RNA", "nCount_RNA")

VlnPlot(parkinsons_sn, group.by = "orig.ident", features = "nFeature_RNA", ncol = 1, pt.size = 0) +
  ggtitle(label = "Number of detected genes per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_nfeat_sn.png"), width = 5, height = 5, dpi = 600)

VlnPlot(parkinsons_sn, group.by = "orig.ident", features = "nCount_RNA", ncol = 1, pt.size = 0) +
  ggtitle(label = "Number of RNA transcripts per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_ncount_sn.png"), width = 5, height = 5, dpi = 600)

VlnPlot(parkinsons_sn, group.by = "orig.ident", features = "percent_mito", ncol = 1, pt.size = 0) +
  ggtitle(label = "Percent mitochondrial reads per cell") +
  theme(axis.text.x.bottom = element_text(angle = 90),
        legend.position = "none")
ggsave(filename = paste0(outdir, "qc_vlnplot_mito_sn.png"), width = 5, height = 5, dpi = 600)

#' 
#' # Loading biomarkers
## -----------------------------------------------------------------------------
neuronal <- c("MAP2", "DCX")
excitatory <- c("RBFOX3", "GRIN1", "HS3ST2", "CAMK2A")
inhibitory <- c("GAD1", "GAD2", "CALB2", "CNR1")
astrocyte <- c("GFAP", "AQP4", "GJA1", "SLC1A3")
oligodendrocyte <- c("PLP1", "MOG", "MBP")
opc <- c("COL9A1", "VCAN", "PDGFRA")
microglia <- c("P2RY12", "FYB1", "CD74")
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

#' 
#' # Assigning cell types
## -----------------------------------------------------------------------------
DotPlot(parkinsons_sn, features = unique(all_markers)) + 
  theme(axis.text.x = element_text(angle = 90))

rep_markers <- c(
  "RBFOX3",
  "GFAP",
  "PLP1",
  "VCAN",
  "FYB1"
)

FeaturePlot(parkinsons_sn, features = rep_markers, ncol = 3)
ggsave(filename = paste0(outdir, "biomarkers_pd_sn.png"), width = 15, height = 10, dpi = 600)

cell_types_df <- FetchData(parkinsons_sn, "seurat_clusters")
cell_types_list <- c(
  "0" = "Oligodendrocytes",
  "1" = "Oligodendrocytes",
  "2" = "Oligodendrocytes",
  "3" = "Oligodendrocytes",
  "4" = "Oligodendrocytes",
  "5" = "Oligodendrocytes",
  "6" = "Astrocytes",
  "7" = "OPCs",
  "8" = "Oligodendrocytes",
  "9" = "Microglia",
  "10" = "Oligodendrocytes",
  "11" = "Astrocytes",
  "12" = "Oligodendrocytes",
  "13" = "Microglia",
  "14" = "OPCs",
  "15" = "Neurons",
  "16" = "Microglia",
  "17" = "Microglia"
)

cell_types_df <- cell_types_df %>%
  mutate(
    cell_type = cell_types_list[seurat_clusters] %>%
      as.factor()
  ) %>%
  dplyr::select(-seurat_clusters)

parkinsons_sn <- AddMetaData(parkinsons_sn, cell_types_df, col.name = "cell_type")

DimPlot(parkinsons_sn, group.by = "cell_type", label = T, label.size = 3, label.box = T, repel = T) +
  ggtitle("Cell type") +
  theme(legend.position = "none")
ggsave(filename = paste0(outdir, "cell_types_pd_sn.png"), width = 4, height = 4, dpi = 700)

DimPlot(parkinsons_sn, group.by = "seurat_clusters", label = T) +
  ggtitle("Clusters")+
  theme(legend.position = "none")
ggsave(filename = paste0(outdir, "clusters_pd_sn.png"), width = 4, height = 4, dpi = 700)

#' 
#' # Adding diagnosis metadata
## -----------------------------------------------------------------------------
parkinsons_sn <- AddMetaData(parkinsons_sn, as.factor(str_split_fixed(parkinsons_sn$orig.ident, "_", 3)[,2]), "diagnosis")

#' 
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
  "SPP1",
  "APOE",
  "CD74",
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
)

cell_data <- parkinsons_sn %>%
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
parkinsons_sn <- SetIdent(parkinsons_sn, value="cell_type")
celltypes <- as.character(unique(Idents(parkinsons_sn)))
deas_celltype <- sapply(as.character(celltypes), function(x) NULL)
for(celltype in as.character(celltypes)) deas_celltype[[celltype]] <- FindMarkers(parkinsons_sn, group.by = "diagnosis", ident.1 = "PD", subset.ident = celltype, logfc.threshold = 0.05)
for(celltype in as.character(celltypes)) deas_celltype[[celltype]]$gene <- rownames(deas_celltype[[celltype]])
for(celltype in as.character(celltypes)) rownames(deas_celltype[[celltype]]) <- NULL
# for(celltype in as.character(celltypes)) deas_celltype[[celltype]] <- deas_celltype[[celltype]][which(deas_celltype[[celltype]]$p_val_adj < 0.01),]

gse_dotplot <- function(dea){
  genelist <- bitr(rownames(dea), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  genelist <- merge(genelist, dea[,"avg_log2FC", drop=F], by.x="SYMBOL", by.y="row.names")
  genelist <- genelist[order(genelist$avg_log2FC, decreasing = T),]
  
  genelist_FC <- sort(genelist$avg_log2FC, decreasing = T)
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

saveRDS(gses_celltype, file = paste0(outdir, "pd_sn_GO.rds"))
# 
# deas_celltype$Neurons

#' 
## -----------------------------------------------------------------------------
dotplot(gses_celltype[["Microglia"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12)) +
  ggtitle("Microglia")
ggsave(filename = paste0(outdir, "microglia_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Astrocytes"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12)) +
  ggtitle("Astrocytes")
ggsave(filename = paste0(outdir, "astrocytes_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Oligodendrocytes"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12)) +
  ggtitle("Oligodendrocytes")
ggsave(filename = paste0(outdir, "oligodendrocytes_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["OPCs"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12)) +
  ggtitle("OPCs")
ggsave(filename = paste0(outdir, "opcs_go.png"), width = 7, height = 5, dpi = 600)

dotplot(gses_celltype[["Neurons"]], showCategory=5, split=".sign") +
  facet_grid(.~.sign) +
  theme(strip.text.x = element_text(size = 12)) +
  ggtitle("Neurons")
ggsave(filename = paste0(outdir, "neurons_go.png"), width = 7, height = 5, dpi = 600)

#' 
