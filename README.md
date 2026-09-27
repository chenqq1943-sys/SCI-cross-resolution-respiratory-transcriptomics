# Cross-resolution transcriptomic analysis reveals scale-dependent mitochondrial respiratory transcriptional remodeling after spinal cord injury

Repository release: **v1.1**

## Authors and affiliations

- Zhaohui Chen (first author), email: chenqq1943@126.com
- Yong Li (corresponding author), email: liyongpuwaike@163.com

1. Jinan University, Guangzhou, China
2. Hunan Children’s Hospital, Changsha, Hunan, China

Funding: Natural Science Foundation of Hunan Province, China (Grant No. 2026JJ81836).

## Data sources

This study uses public GEO datasets GSE5296, GSE234774, and GSE319931. Large GEO raw files are not redistributed here; obtain them directly from GEO. The repository contains final analysis code, compact metadata, predefined gene programs, manuscript-supporting derived data, final figures, and supplementary source tables.

## Analytical scope

The repository supports bulk transcriptomics, single-nucleus RNA-seq library-level pseudobulk, spatial transcriptomics, and chronic regional bulk RNA-seq analyses. It evaluates predefined Complex I assembly and aerobic-respiration transcriptional programs using biological replicate-aware inference and prespecified sensitivity analyses. Program scores are transcriptional summaries and do not directly measure oxygen consumption, ATP production, or mitochondrial function.

## Recommended script order

1. `scripts/01_GSE5296/01_rma_qc.R`
2. `scripts/01_GSE5296/02_probe_collapse_and_programs.R`
3. `scripts/02_GSE234774_snRNA/01_build_fulltranscriptome_pseudobulk.py`
4. `scripts/02_GSE234774_snRNA/02_normalize_score_compare.R`
5. `scripts/03_GSE234774_spatial/01_spatial_score_compare.py`
6. `scripts/04_GSE319931/01_score_programs.R`
7. `scripts/04_GSE319931/02_mixed_model.R`
8. `scripts/06_figures/01_plot_locked_figures.R`

The two release-v1 supplementary-table utilities are retained under `audit_archive/legacy_scripts/05_supplementary_tables/`. They are historical packaging/checking utilities, not the recommended way to rebuild the final v1.1 tables. The archived verifier is used only with the archived v1 table package for the 44-item mechanical audit.

Inputs are GEO-derived expression/count matrices and metadata prepared from the three public accessions, together with `gene_programs/`. Large raw inputs must be downloaded separately. Outputs include compact tables under `derived_data/`, final figure assets under `figures/`, and final table source files under `supplementary_source_tables/`.

The scripts are preserved byte-for-byte from Repository Release v1 because this release does not change algorithms, parameters, statistics, or outputs. They use relative project paths. Some downloads and platform-specific preprocessing require manual preparation and are not represented as a single one-click workflow. The historical figure script retains the internal filename `FigureS3` for the asset now submitted as Figure S2; the packaged file under `figures/supplementary/` has the final S2 name and unchanged bytes.

## Derived data and provenance

`derived_data/` contains the compact outputs needed to support the manuscript and final Supplement. `metadata/` records dataset/sample provenance, and `environment/` records available software/session information. Historical development and audit material is isolated under `audit_archive/` and is not part of the recommended workflow.

## Licensing

- **Original analysis code/scripts:** MIT License; see `LICENSE-CODE`.
- **Eligible author-generated non-code materials:** Creative Commons Attribution 4.0 International (CC BY 4.0); see `LICENSE-DATA-DOCS`.

Unless otherwise indicated, author-generated non-code materials in this repository are licensed under CC BY 4.0 to the extent that the authors own the relevant rights and may legally grant that license. GEO raw data, other third-party data, third-party software/code, figures, annotations, database content, and any material outside the authors' relicensing authority are excluded and remain subject to their original terms. Large GEO raw data are not redistributed here.
