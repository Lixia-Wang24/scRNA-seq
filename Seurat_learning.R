###Step 1 ----------------------- Update Seurat to the latest version-------------------------------------------------###
# Uninstall old Seurat
remove.packages(c("Seurat", "SeuratObject"))
# Install old Seurat
install.packages('Seurat')

setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))

# Install the remotes package to install packages that are not on CRAN
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac') # analysis of single-cell chromatin data
remotes::install_github("satijalab/seurat-data", quiet = TRUE) # automatically load datasets pre-packaged as Seurat objects
remotes::install_github("satijalab/azimuth", quiet = TRUE) # local annotation of scRNA-seq and scATAC-seq queries across multiple organs and tissues
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE) # enables use of additional integration and differential expression methods

###Step 2 ----------------------- Guided tutorial — 2,700 PBMCs -------------------------------------------------###
# load Seurat
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)

###----------------------- Load data and create a Seurat object -------------------------------------------------###

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "/Users/lixia/Data/data/scRNAseq/Seurat_PBMC3k/filtered_gene_bc_matrices/hg19/") #Read10X_h5() to read h5 file format (output of the cellranger)
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# Filter gene existing in fewer than 3 cells, and the cell has fewer than 200 genes.
# The number of unique genes and total molecules are automatically calculated.

# view & check
str(pbmc.data) # overview
dim(pbmc.data) # view dimension, pbmc.data is sparse matrix without 0 value
pbmc.data[1:3, 1:3] # view the top3 colums and rows
head(rownames(pbmc.data)) # view gene ensemble name
head(colnames(pbmc.data)) # view cell ID

pbmc # overview, pmbmc is object of Seurat
dim(pbmc) #view dimension
head(rownames(pbmc)) 
head(colnames(pbmc))
LayerData(pbmc, layer = "counts")[1:3, 1:3] # view count slot
head(pbmc@meta.data) # view meta.data
View(pbmc@meta.data) # Open the complete metadata table in a new window
table(pbmc$orig.ident) # If data has multiple samples, use this command to see how many cells are in each sample:

###----------------------- Quality Control -------------------------------------------------------------------###

# Calculate the mitochondrial contamination
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") # add "percent.mt" columns to object metadata,all genes starting with MT-
# An excessively high proportion of mitochondria (e.g., >20%) usually indicates a damaged or dead cell.
# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # ncol arranges the 3 images in 3 columns (i.e., 1 row and 3 columns).
# Visualize feature-feature relationships
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# fFilter cells that have unique feature counts over 2,500 or less than 200, >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# Note: This 5% is an empirical value. For metabolically active tissues such as the heart and liver, or for frozen samples, 
# this threshold may need to be relaxed to 10% or even 20%.
# Regarding the differences in UMI numbers, our strategy is to "correct it" rather than "delete it".

###----------------------- Normalization -------------------------------------------------------------------###
# By default, Seurat normalizes the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# Normalized values are stored in pbmc[["RNA"]]$data.
# Global-scaling relies on an assumption that each cell originally contains the same number of RNA molecules.
# SCTransform() doesn't have this assumption.

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000) # same to pbmc <- NormalizeData(pbmc), by default
LayerData(pbmc, layer = "data")[1:3, 1:3] # view data slot

###----------------------- Feature selection --------------------------------------------------------------###
# The subset of features that exhibit high cell-to-cell variation is calculated by modeling the mean-variance relationship inherent in single-cell data.
# By default, 2,000 features per dataset is returned.

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000) # 2000 is empirical value.
# vst = Variance Stabilizing Transformation
# vst 不仅仅是算方差，因为基因表达量和方差通常呈正相关（表达量越高的基因，方差天然就越大）。
# VST 方法会建立一个模型，根据基因的平均表达量来预测其预期的方差，然后找出那些“实际方差”远大于“预期方差”的基因。这些才是真正生物学上有意义的差异基因。

# The number of variable features
top_genes <- VariableFeatures(pbmc)
length(top_genes)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# View detailed statistics for the top 10 hypervariable genes.
head(HVFInfo(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)  # Functions included in the Seurat package
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) # repel prevents text overlap.
plot1 + plot2

###----------------------- Scaling ----------------------------------------------------------------------###
# Scaling is a linear transformation （Z-score Transformation）, which is a standard pre-processing step prior to dimensional reduction.
# By default, only variable features are scaled.
# Shifts the expression of each gene, so that the mean expression across cells is 0, the variance across cells is 1,
# This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate.
# The results of this are stored in pbmc[["RNA"]]$scale.data.

pbmc <- ScaleData(pbmc)
LayerData(pbmc, layer = "scale.data")[1:3, 1:3] # view scale.data slot

#If you want to create a heatmap to show the expression differences of non-hypervariable genes later, you can use all.gene for scaling.
# pbmc <- ScaleData(pbmc, features = all.genes)

# remove unwanted sources of variation，SCTransform() is recomended for this step.
# pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt") #regress out heterogeneity associated with mitochondrial contamination (cell cycle stage,nCount_RNA).

###----------------------- Linear dimensional reduction----------------------------------------------------------------------###
# By default, only the previously determined variable features are used as input.
# PCA outputs a list of genes with the most positive and negative loadings, representing modules of genes that exhibit either correlation (or anti-correlation) across single-cells

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc)) #写不写object = 没影响
# The `VariableFeatures(object = pbmc)` function itself does not return a value, but rather a "name".
# However, the `RunPCA` function uses these names to find the corresponding "scaled value" in `scale.data` for calculation.

# Overview of the PCA object.
pbmc[["pca"]]
print(pbmc[["pca"]])
print(pbmc[["pca"]], dims = 1:3, nfeatures = 5) # Top5 genes with the largest absolute weight.
# View the PCA coordinates (Cell Embeddings) of the cells.
pca_coords <- Embeddings(pbmc, reduction = "pca")
head(pca_coords[, 1:5])
# Extract the first 10 PC coordinates from the first 5 cells.
Embeddings(pbmc, reduction = "pca")[1:5, 1:10]

# View gene weights/contributions (Feature Loadings)
gene_loadings <- Loadings(pbmc, reduction = "pca")
head(gene_loadings[, 1:5])
# View the standard deviation for each PC.
head(Stdev(pbmc, reduction = "pca"), 10)


# Examine and visualize PCA results a few different ways
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") # by default,nfeatures = 30.

#View the distribution of cells in the PCA space (scatter plot),(by default, PC1 vs PC2)
DimPlot(pbmc, reduction = "pca") + NoLegend() # remove figure legend

# Verify signal strength (heatmap),(is it a real signal or noise?)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) # random 500 cells, TOP30 gene, 50% pos-gene and 50% neg-gene
# DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE, nfeatures = 50) 
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) 

# The standard deviation of each principal component's contribution helps us find the "inflection point".
# Elbow point determines how many components should we choose to include?
ElbowPlot(pbmc)  #the first 10 PCs

###----------------------- Clster the cells ----------------------------------------------------------------------###
# Find cluster by a graph-based clustering approach
pbmc <- FindNeighbors(pbmc, dims = 1:10) # 1.Embeddings of 10 PCs;2.K-Nearest Neighbor (KNN);3.Weighted SNN Graph;
pbmc <- FindClusters(pbmc, resolution = 0.5) # 4.Modularity Optimization 
#resolution = 0.4-1.2 for datasets of around 3K cells

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

###----------------------- Run non-linear dimensional reduction (UMAP/tSNE) --------------------------------------###
# visualization,cells with very similar gene expression profiles co-localize, 
# UMAP remains partial global structure; 
pbmc <- RunUMAP(pbmc, dims = 1:10)# 1.PCA matrix,2.KNN graph,3.Fuzzy Graph(P),4.initial 2D,5.optimization,6.UMAP coords

DimPlot(pbmc, reduction = "umap", label = TRUE)

saveRDS(pbmc, file = "/Users/lixia/Data/data/scRNAseq/Seurat_PBMC3k/output/pbmc_tutorial.rds")

# t-SNE focuses only on local structure.
# pbmc <- RunTSNE(pbmc, dims = 1:10) 
# DimPlot(pbmc, reduction = "tsne", label = TRUE)

pbmc <- RunTSNE(
  pbmc, 
  dims = 1:10, 
  perplexity = 30,  # 默认是 30。细胞少(比如<200)时调小，细胞多时调大。
  max_iter = 1000,  # 迭代次数，默认 1000。没跑开的话可以设大一点。
  seed.use = 42     # 锁定随机种子，保证下次跑出来的图一模一样
)

DimPlot(pbmc, reduction = "tsne", 
        label = TRUE,       # 在图上直接标出 Cluster ID 数字
        pt.size = 0.5,      # 调整点的大小
        label.size = 5      # 调整标签文字大小
) + NoLegend()              # 如果图例太占位置，可以把图例关掉


###----------------------- Finding differentially expressed features (cluster biomarkers) ------------------------###
# Find all maker of all cluster (cluster vs all other cluster)
# 1.pct; 2.avg_log2FC; 3.p-val; 4.p_val_adj
# all_markers <- FindAllMarkers(pbmc, 
#                              only.pos = TRUE,        # positive gene only
#                              min.pct = 0.25,         # >25% cells have this gene
#                              logfc.threshold = 0.25) # ident.1/ident.2

# ?FindMarkers # view the min.pct (0.01) and logfc.threshold (0.1) parameters

# The top5 positive marker for each cluster
library(dplyr)
# top5_markers <- all_markers %>%
#  group_by(cluster) %>%
#  top_n(n = 5, wt = avg_log2FC)

# print(top5_markers)

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Examine and visualization
VlnPlot(pbmc, features = "CCR7") 
FeaturePlot(pbmc, features = "CCR7") # UMAP

#  Find markers of particular clusters
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2) # find all markers of cluster 2, by default ident.2 = all others.
head(cluster2.markers, n = 5)
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3)) # find all markers distinguishing cluster 5 from clusters 0 and 3
head(cluster5.markers, n = 5)

# Instead of using statistical tests (such as Wilcoxon) to calculate p-values, use ROC curve analysis to assess the gene's taxonomic ability.
# For genes with extremely low expression levels but significant differences, the AUC may not be high.
# ROC:the "cleanest" Marker; p-val:the best balance between speed and accuracy
# MAST: rigorous;
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
head(cluster0.markers)

# visualization of marker expression
# RidgePlot(), CellScatter(), and DotPlot() 

VlnPlot(pbmc, features = c("MS4A1", "CD79A")) # shows expression probability distributions across clusters
# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A")) # visualizes feature expression on a tSNE or PCA plot

pbmc.markers %>%     #expression heatmap
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

###----------------------- Assigning cell type identity to clusters -----------------------------------------###
#clusters annotation by canonical markers 
# Replace the numerical designations (Cluster 0, 1, 2...) with biologically meaningful names (such as "T cell", "B cell"...).
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc) # give name a cluster ID
pbmc <- RenameIdents(pbmc, new.cluster.ids) # replace the name of idents from cluster IDs with names
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

# figure optimization
library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 12)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot
ggsave(filename = "/Users/lixia/Data/data/scRNAseq/Seurat_PBMC3k/output/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)

