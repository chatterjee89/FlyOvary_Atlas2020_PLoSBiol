# 2020_scAtlas_Dmel_ovary

Single-cell RNA-seq analysis pipeline for a *Drosophila melanogaster* ovary cell atlas. The script implements cell clustering with **Seurat 2.3.4** followed by pseudotime trajectory inference with **Monocle v2**.

---

## Pipeline overview

### 1. Seurat 2.3.4 — clustering

| Step | Details |
|------|---------|
| Data loading | 10X Chromium output via `Read10X` |
| QC filtering | nGene 775–2200 · nUMI < 18,000 · mitochondrial fraction < 1% |
| Normalisation | LogNormalize (scale factor 10,000) |
| Cell-cycle regression | S and G2/M scores computed; `CC.Difference = S.Score − G2M.Score` regressed out together with nUMI and percent.mito |
| Variable genes | `FindVariableGenes` with default VMR dispersion |
| Dimensionality reduction | PCA (100 PCs); elbow plot used to select PCs 1–30 |
| Clustering | SLM algorithm (algorithm 3), resolution 4.5; `clustree` used to evaluate stability across resolutions |
| Embedding | UMAP (`n_neighbors = 20`, `min_dist = 0.35`) |

### 2. Monocle v2 — pseudotime trajectory

| Step | Details |
|------|---------|
| Input | Seurat cluster subset imported via `importCDS` |
| Preprocessing | Size factor and dispersion estimation; genes expressed in ≥ 20 cells retained |
| Initial embedding | tSNE for density-based clustering |
| Clustering | `clusterCells` with user-defined ρ and δ thresholds (visualised with `plot_rho_delta`) |
| Ordering genes | Differential gene test across clusters (`~Cluster`); top genes selected after manual inspection |
| Trajectory | DDRTree dimensionality reduction → `orderCells` for pseudotime assignment |
| Visualisation | Pseudotime, cluster, and State overlays on trajectory; pseudotime on tSNE |

---

## Dependencies

```r
library(Seurat)       # v2.3.4
library(monocle)      # v2
library(ggplot2)
library(ggpubr)
library(dplyr)
library(R.utils)
library(clustree)
library(cowplot)
library(RColorBrewer)
library(corrplot)
library(org.Dm.eg.db)
library(topGO)
```

---

## Input files

| File | Description |
|------|-------------|
| `~/cellranger_output/` | 10X Cell Ranger output directory (barcodes, genes, matrix) |
| `~/cell_cycle_genes.txt` | Cell cycle gene list (provided in paper supplementary); lines 1–68 = G2/M genes, lines 69–124 = S-phase genes |

---

## Key parameters to adjust

- **`###` placeholders** in the Monocle section mark values that require manual tuning per dataset:
  - Seurat cluster IDs to subset for trajectory analysis
  - Number of tSNE dimensions (`num_dim`)
  - ρ and δ thresholds for density-based clustering
  - Number of ordering genes

---

## Citation

If you use this pipeline, please cite the associated publication and the tools:

- **Seurat**: Butler et al. *Nature Biotechnology* (2018)
- **Monocle v2**: Qiu et al. *Nature Methods* (2017)
- **clustree**: Zappia & Oshlack *GigaScience* (2018)
