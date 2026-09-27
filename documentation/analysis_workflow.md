# Analysis workflow

The recommended computational sequence is listed in the root README. GSE5296 uses array-level QC/RMA and predefined respiratory-program scoring; GSE234774 snRNA uses full-transcriptome library-level pseudobulk, filterByExpr, TMM, logCPM, and post-normalization program extraction; the spatial analysis summarizes spots/sections at biological-sample level for inference; GSE319931 uses repeated regional measures within animal. Large GEO inputs are not included.

The final plotting script reads released derived data and final supplementary source tables. Historical v1 table-building/checking utilities and their legacy table package are isolated under `audit_archive/`.
