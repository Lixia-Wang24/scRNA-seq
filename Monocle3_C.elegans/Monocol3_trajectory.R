R.version
BiocManager::version()

########################################## Install basic dependencies #########################################################
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# install a few Bioconductor dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

# install monocle3
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

remotes::install_version("grr", version = "0.9.5")

devtools::install_github('cole-trapnell-lab/monocle3') # 3：none

########################################## Load data #########################################################
# Set work place
setwd("/Users/lixia/Data/data/scRNAseq/Monocle3_C.elegans")
# Creat reult directory
dir.create("R.results", showWarnings = FALSE)

library(monocle3)

# Load the pre-processed scRNA-seq data for PBMCs
expression_matrix <- readRDS("packer_embryo_expression.rds")
cell_metadata <- readRDS("packer_embryo_colData.rds")
gene_annotation <- readRDS("packer_embryo_rowData.rds")
print(expression_matrix)


cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

########################################## preprocessing #########################################################
# Normalize, feature selection, scaling, PCA, Residual Regression, batch effect correction
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "batch", 
                 residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + 
                 bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")

############################### Reduce dimensionality and visualize the results ###################################

# Reduce dimensionality and visualize the results (UMAP dimension)， which is different from seurat.
cds <- reduce_dimension(cds)  # max_components = 2， default
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")
# plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "batch")

# visualize how individual genes vary along the trajectory
ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6", "ceh-36", "ham-1")
plot_cells(cds,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# Visualize the overall performance of this gene group.
plot_genes_by_group(cds, 
                    markers = ciliated_genes, 
                    group_cells_by = "cell.type")

########################################## Cluster ################################################################
# Clusters and partitions
# If there is almost no "path" between the two large groups of cells, 
# they will be cut off and divided into different partitions.
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

########################################## Learn the trajectory graph ##############################################
# Based on intercellular similarity, a tree is drawn on UMAP linking transitions in cell states.
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

########################################## Order the cells in pseudotime ###########################################
# Market experiment time.
plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)

#  Manually pick the starting point (root) of the trajectory on map 
cds <- order_cells(cds)

# Visualize pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

#----------------------------------------------------————————————————————————————————————————————————————------------
# Pick the root of the trajectory programmatically,
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="130-170"){
  cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

# Rank and calculate pesudotime
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

# Visualize pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

########################################## Subset cells by branch ##################################################

# Subset by picking start rot and end rot manually.
cds_sub <- choose_graph_segments(cds)

#------------------------------------------------------------————————————————————————————————————————————————-------
# # Subset by picking start rot and end rot programmatically.
# Define the start and end nodes on the trajectory graph
start_node <- "Y_1"
end_node <- "Y_20"

# Extract the principal graph object (Minimum Spanning Tree) for the UMAP projection
mst <- principal_graph(cds)[["UMAP"]]

# Calculate the shortest path between the two nodes using the igraph library
# $vpath[[1]] extracts the sequence of vertices (nodes) along that path
path_nodes <- igraph::shortest_paths(mst, from = start_node, to = end_node)$vpath[[1]]

# Get the character names of the nodes (e.g., "Y_1", "Y_5", "Y_20")
path_node_names <- names(path_nodes)

# Retrieve the mapping of each cell to its nearest principal graph node
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex

# Identify indices of cells that project onto the nodes found in the shortest path
cells_in_path <- which(closest_vertex %in% path_node_names)

# Subset the original CDS object to keep only the cells belonging to this specific path
cds_sub2 <- cds[, cells_in_path]

########################################## Working with 3D trajectories ##############################################
# 3D trajectories analysis based 2D processed data.
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

# Visualize 3D trajectories
cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
cds_3d_plot_obj

########################################## Finding genes that change as a function of pseudotime  ##############################################

# Find the genes that are differentially expressed on the different paths through the trajectory
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4) # all gene, all partition, 2D principle graph,globle weight
#Select genes with significant dynamic changes
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05)) #Sort by WBGene number in ascending order.

# Sort by morans_I in descending order
library(tidyverse)
pr_deg_ids <- ciliated_cds_pr_test_res %>% 
  filter(q_value < 0.05) %>% 
  arrange(desc(morans_I)) %>% 
  row.names()

head(pr_deg_ids) 

# Convert gene ID into gene symble
gene_info <- rowData(cds)
pr_deg_names <- gene_info[pr_deg_ids, "gene_short_name"]
head(pr_deg_names)

# Visualize highly significant genes
plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

plot_cells(cds, genes=c("gcy-23", "ugt-56", "D2096.9", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)


########################################## Collect the trajectory-variable genes into modules  ##############################################

# Scattered differentially expressed genes are automatically categorized into several "functional modules" with similar expression patterns.
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
head(gene_module_df)

#------------------------------------------------------------————————————————————————————————————————————————-------

# Visualize  the aggregate module scores within each group of cell types

# Create an R data frame containing cell_ID and cell types.
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell.type) 

# Combined gene modules x combined cell clusters
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)

# Change the row names of the data frame
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

# Visualze  gene-module expression by heatmap
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

# Visualze gene-module expression on UMAP
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(27, 10, 7, 30)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

########################################## gene expression along pseudotime ##############################################

# Transforming static cell expression data into a dynamic process that flows over time.
# The dynamics of a small set of genes，in particular cell type， as a function of pseudotime

AFD_genes <- c("gcy-8", "dac-1", "oig-8")

AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,
                       colData(cds)$cell.type %in% c("AFD")]

AFD_lineage_cds <- order_cells(AFD_lineage_cds)


plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)

########################################## gene expression in Subset cells by branch ##################################################

# select a subset
cds_subset2 <- choose_cells(cds)
# Find the genes from subset cells that are differentially expressed on the different paths through the trajectory
subset_pr_test_res2 <- graph_test(cds_subset2, neighbor_graph="principal_graph", cores=4)
# significant genes
pr_deg_ids2 <- row.names(subset(subset_pr_test_res2, q_value < 0.05))
# Find gene modules
gene_module_df2 <- find_gene_modules(cds_subset2[pr_deg_ids2,], resolution=0.001)
# Modules expression
agg_mat2 <- aggregate_gene_expression(cds_subset2, gene_module_df2)

# Gene modules are sorted by expression similarity.
module_dendro2 <- hclust(dist(agg_mat2))
gene_module_df2$module <- factor(gene_module_df2$module, 
                                levels = row.names(agg_mat2)[module_dendro2$order])

# Output the sorted module numbers
row.names(agg_mat2)[module_dendro2$order]

# Visualze gene-module expression on UMAP  
plot_cells(cds_subset2,
           genes=gene_module_df2,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

#------------- label ranked module names ---------------
# Rename the modules according to the order you found them.
new_module_names <- c("1" = "M1", "4" = "M2", "2" = "M3", "7" = "M4", 
                      "6" = "M5", "3" = "M6", "8" = "M7", "5" = "M8", 
                      "9" = "M9", "10" = "M10", "11" = "M11")

gene_module_df2$module_renamed <- new_module_names[as.character(gene_module_df2$module)]
gene_module_df2$module_renamed <- factor(gene_module_df2$module_renamed, 
                                         levels = unname(new_module_names))

# Redraw the graph using the new column module_renamed
plot_cells(cds_subset2, genes=gene_module_df2 %>% mutate(module = module_renamed),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

