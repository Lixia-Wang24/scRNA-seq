# Single-Cell RNA-Seq Analysis Workflow

## Overview
This repository provides a complete workflow for **single-cell RNA sequencing (scRNA-seq)** data analysis, covering three major analytical approaches to unlock cellular insights:

* **Seurat Analysis** – Standard clustering and cell type identification.
* **Monocle3 Trajectory Analysis** – Pseudotime ordering and lineage reconstruction.
* **scVelo Velocity Analysis** – RNA velocity for directional dynamics.

These tools work together to provide comprehensive insights into cellular heterogeneity, developmental trajectories, and dynamic processes.

---

## Features

### 🧬 Seurat Analysis (`Seurat_PBMC3k/Seurat_analysis.R`)
Used for the fundamental processing of scRNA-seq data.
* ✅ Quality control (QC) and filtering
* ✅ Normalization and feature selection
* ✅ Dimensional reduction (PCA, UMAP, t-SNE)
* ✅ Graph-based clustering
* ✅ Differential expression analysis (DEA)
* ✅ Cell type annotation & visualization



### 🌊 Trajectory Analysis (`Monocle3_C.elegans/Monocle3_trajectory.R`)
Focuses on the transitions between cellular states.
* ✅ Pseudotime trajectory inference
* ✅ Branch point detection and lineage reconstruction
* ✅ Gene module identification
* ✅ 3D trajectory visualization
* ✅ Cell fate analysis and dynamic gene expression patterns



### 🚀 Velocity Analysis (`RNA_velocity/scVelo_velocity.py`)
Predicts the future state of individual cells using the ratio of unspliced to spliced mRNA.
* ✅ RNA velocity estimation (Dynamical & Stochastic models)
* ✅ Velocity graph construction
* ✅ Latent time/Pseudotime inference
* ✅ Driver gene identification
* ✅ Directional trajectory visualization & confidence assessment


---

## Installation & Dependencies

### R Dependencies

```r
# Install Seurat ecosystem
install.packages('Seurat')
install.packages(c("dplyr", "patchwork", "ggplot2"))

# Install Monocle3 via Bioconductor and GitHub
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))

devtools::install_github('cole-trapnell-lab/monocle3')
```

Python Environment (scVelo)
```bash
# Create and activate environment
conda create -n scrna python=3.10
conda activate scrna

# Install scVelo and Scanpy
pip install scvelo scanpy numpy pandas matplotlib
```

## Analysis Workflow

The modules are designed to be used sequentially or as independent pipelines:

| Phase | Module | Primary Objective |
| :--- | :--- | :--- |
| **I** | **Seurat Preprocessing** | QC filtering, normalization, and identifying cluster-specific markers. |
| **II** | **Trajectory Inference** | Mapping developmental "pseudotime" and identifying lineage-driver genes. |
| **III** | **RNA Velocity** | Determining the vector and speed of cell state transitions. |

---

## Summary of Expected Metrics

| Metric | Target / Expectation |
| :--- | :--- |
| **Cluster Separation** | Distinct UMAP/t-SNE groupings matching known cell types. |
| **Trajectory Path** | A continuous path from progenitor cells to terminal states. |
| **Velocity Stream** | Vectors aligning with known biological differentiation axes. |
| **Correlation** | High consistency between pseudotime (Monocle3) and latent time (scVelo). |
