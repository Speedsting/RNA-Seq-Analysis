# RNA-Seq differential expression analysis

**Author:** Elijah Nichols
**Date:** 2025-08-12

This repository contains an R-based RNA-Seq analysis pipeline that compares disease vs healthy samples and reports differential expression results using DESeq2, edgeR, and limma-voom. The analysis generates volcano plots, PCA, heatmaps, and gene lists. The original R script used to produce the results is included in this repository. :contentReference[oaicite:1]{index=1}

### Quick links
- Notebook (reproducible report): `RNA_Seq_Analysis.ipynb`
- Full instructions & reproduction steps: `README.md`
- Raw counts (GCT): `samples_raw.gct`

### Highlights
- Differential expression via DESeq2, edgeR, limma-voom
- Volcano plot, PCA on rlog-transformed counts, top variance gene heatmap
- Output: CSVs of DE gene lists

If you want to reproduce this analysis on your machine, open the README for installation and run instructions.