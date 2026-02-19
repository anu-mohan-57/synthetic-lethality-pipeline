# Synthetic Lethality Prediction Pipeline

A computational pipeline for predicting cell type-specific synthetic lethal (SL) gene pairs in human cancers using single-cell RNA sequencing data.

---

## Overview

Synthetic lethality occurs when the simultaneous disruption of two genes causes cell death, while disrupting either gene alone is non-lethal. This property makes SL gene pairs promising therapeutic targets - disabling the partner gene of an already-mutated cancer gene can selectively kill tumour cells while sparing normal tissue. The BRCA/PARP inhibitor relationship is the best-known clinical example of this principle in action.

Existing computational methods for SL prediction rely heavily on bulk RNA-seq data, which averages gene expression across heterogeneous cell populations and misses cell-type-specific interactions. This pipeline addresses that gap by using **single-cell RNA sequencing (scRNA-seq)** data to predict SL pairs at cell-type resolution, integrating druggability and gene network information to prioritise clinically actionable candidates.

---

## Method

The pipeline runs in five stages:

1. **Binary thresholding** — scRNA-seq expression values are converted to presence (1) or absence (0) per gene per cell
2. **Co-absence matrix construction** — gene pairs that are never simultaneously absent across a cell type are flagged as SL candidates; parallelised with Python's `multiprocessing` library for genome-scale data
3. **Druggability filtering** — candidates are cross-referenced against [DGIdb](https://www.dgidb.org/) to retain only therapeutically actionable pairs
4. **Network proximity filtering** — only pairs with a direct interaction (network distance = 1) in [BioGRID](https://thebiogrid.org/) are retained
5. **Novelty filtering** — known SL pairs from [SynLethDB](https://synlethdb.sist.shanghaitech.edu.cn/), [SLorth](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6602458/), and [SLKB](https://slkb.osubmi.org/) are removed, leaving only novel candidates

---

## Data Sources

| Database | Purpose |
|---|---|
| [Tabula Sapiens](https://tabula-sapiens-portal.ds.czbiohub.org/) | scRNA-seq atlas across 24 human tissues |
| [Human Protein Atlas](https://www.proteinatlas.org/) | Transcript expression across 81 cell types |
| [DGIdb](https://www.dgidb.org/) | Drug-gene interaction annotations |
| [BioGRID](https://thebiogrid.org/) | Gene interaction network |
| [SynLethDB / SLorth / SLKB](https://synlethdb.sist.shanghaitech.edu.cn/) | Known SL pairs for novelty filtering |

---



## Repository Contents

- `matrix_generation_SL.ipynb` — full pipeline: binary thresholding, co-absence matrix construction 
- `hepatocyte.csv` — predicted SL candidate pairs for hepatocyte cell type
- `TS_pipeline.ipynb` - Tabula Sapiens–based preprocessing and preparation of scRNA-seq expression matrices.

- `HPA_pipeline.ipynb` - Processing and integration of Human Protein Atlas transcript expression data across 81 cell types.

- `druggable_genes.ipynb` - Identification and curation of druggable genes from DGIdb for downstream filtering.

- `GRN.ipynb` - Gene regulatory / interaction network construction and processing used for proximity-based filtering of candidate SL pairs.

---

## Requirements

```
pandas
numpy
scipy
h5py
matplotlib
multiprocessing
```

---

## About

This pipeline was developed as part of a Master's thesis at the Indian Institute of Science Education and Research (IISER) Thiruvananthapuram, April 2025.

**Title:** Computational Prediction of Cell Type-Specific Synthetic Lethality in Human Cancers  
**Author:** Anugreha Mohan  
**Supervisors:** Dr. Sanu Shameer and Dr. Kamalakannan Vijayan

Before running the pipeline, please update the g variables in the script to match your local file structure.


