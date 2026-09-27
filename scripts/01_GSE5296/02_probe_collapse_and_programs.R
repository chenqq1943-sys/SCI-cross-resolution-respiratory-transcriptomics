## GSE5296 RMA replication test  (purely computational; no biological interpretation)
# Use the caller's existing R library paths.
suppressMessages({
  library(GSVA); library(org.Mm.eg.db); library(GO.db); library(AnnotationDbi)
})
options(stringsAsFactors=FALSE)
base <- if (length(commandArgs(trailingOnly=TRUE))) commandArgs(trailingOnly=TRUE)[[1]] else getwd(); setwd(base)
logf <- file("_replication_log.txt","wt"); sink(logf, type="output")
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
cat("START", ts(), "\n")

## ---------- 1. data integrity ----------
mat <- read.csv(gzfile("GSE5296_RMA_probe_expression.csv.gz","rt"), header=TRUE,
                check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
mat <- as.matrix(mat)
colnames(mat) <- gsub('"','',colnames(mat)); rownames(mat) <- gsub('"','',rownames(mat))
nprobe <- nrow(mat); nsamp <- ncol(mat)
cat("RMA matrix dim:", nprobe, "x", nsamp, "\n")
cat("sample IDs (96):", paste(colnames(mat), collapse=" "), "\n")
ann <- read.csv("GSE5296_RMA_sample_annotation.csv", check.names=FALSE, stringsAsFactors=FALSE)
matched <- all(colnames(mat) == ann$GSM)
cat("annotation 100% order-matched:", matched, "\n")
cat("duplicate samples:", anyDuplicated(colnames(mat)), "\n")
cat("duplicate probe IDs:", anyDuplicated(rownames(mat)), "\n")
cat("NA present:", anyNA(mat), " | Inf present:", any(is.infinite(mat)), "\n")
if (nsamp != 96 || !matched) { cat("INTEGRITY FAILED; stopping\n"); quit(save="no") }

## ---------- 2. probe->gene via GPL1261 ----------
conG <- file("GPL1261.annot.gz","rt"); gpl <- readLines(conG); close(conG)
st <- grep("^!platform_table_begin", gpl); en <- grep("^!platform_table_end", gpl)
rows <- gpl[(st+1):(en-1)]
tab <- read.table(text=rows, sep="\t", header=TRUE, quote="", comment.char="",
                  fill=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
cn <- colnames(tab)
sym_col <- cn[grepl("Gene", cn, ignore.case=TRUE) & grepl("ymbol", cn, ignore.case=TRUE)][1]
id_col  <- cn[grepl("Gene", cn, ignore.case=TRUE) & grepl(" ID$|ID", cn) & !grepl("ymbol",cn)][1]
cat("GPL1261 symbol col:", sym_col, "| entrez col:", id_col, "\n")
pg <- data.frame(probe=as.character(tab$ID), symbol=as.character(tab[[sym_col]]),
                 entrez=as.character(tab[[id_col]]), stringsAsFactors=FALSE)
pg$symbol <- ifelse(is.na(pg$symbol) | pg$symbol=="" | pg$symbol=="---", NA, pg$symbol)
## restrict to probes actually in matrix
mat_probes <- rownames(mat)
pg <- pg[pg$probe %in% mat_probes, ]
total_probe <- nrow(mat)
mapped <- !is.na(pg$symbol)
cat("total probes in matrix:", total_probe, "\n")
cat("probes with annotation rows:", nrow(pg), "\n")
cat("mapped probes (symbol present):", sum(mapped), "\n")
cat("unmapped probes:", sum(!mapped), "\n")

## mean expression per probe across samples (for tie-break / selection)
probe_mean <- rowMeans(mat)

## select representative probe per symbol = highest mean expression
pgm <- pg[mapped, ]
pgm$mean_expr <- as.numeric(probe_mean[as.character(pgm$probe)])
ord <- order(pgm$symbol, -pgm$mean_expr)
pgm <- pgm[ord, ]
sel <- pgm[!duplicated(pgm$symbol), ]   # first per symbol after desc mean
colnames(sel)[colnames(sel)=="mean_expr"] <- "Mean_expression_selected_probe"
unique_symbols <- length(unique(pgm$symbol))
cat("unique gene symbols mapped:", unique_symbols, "\n")
## multi-probe distribution
multi <- table(table(pgm$symbol))
cat("multi-probe distribution (n_probes -> n_genes):\n")
print(multi)
cat("1-probe genes:", multi["1"], " | 2-probe genes:", ifelse(is.na(multi["2"]),0,multi["2"]), "\n")

## build gene-level matrix (selected probes)
geneMat <- mat[sel$probe, ]
rownames(geneMat) <- sel$symbol
geneMat <- geneMat[!duplicated(rownames(geneMat)), ]

## probe selection table output
cnt <- as.integer(table(pgm$symbol)); names(cnt) <- rownames(table(pgm$symbol))
v1 <- as.character(sel$symbol)
v2 <- as.character(sel$probe)
v3 <- as.integer(cnt[as.character(sel$symbol)])
v4 <- round(rowMeans(mat[as.character(sel$probe), , drop=FALSE]), 4)
v5 <- as.character(sel$entrez)
cat("sel rows:", nrow(sel), "| v1:", length(v1), "v2:", length(v2), "v3:", length(v3), "v4:", length(v4), "v5:", length(v5), "\n")
stopifnot(length(v1)==length(v2), length(v2)==length(v3), length(v3)==length(v4), length(v4)==length(v5))
sel_out <- data.frame(GeneSymbol=v1, Selected_probe=v2,
                      Number_of_candidate_probes=v3,
                      Mean_expression_selected_probe=v4,
                      Entrez=v5, stringsAsFactors=FALSE)
write.csv(sel_out, "GSE5296_RMA_probe_to_gene_selected.csv", row.names=FALSE)
cat("wrote GSE5296_RMA_probe_to_gene_selected.csv\n")

## ---------- 4. Sirt6 single-gene ----------
cat("Sirt6 selected probe:", ifelse("Sirt6" %in% rownames(geneMat),
      sel$probe[sel$symbol=="Sirt6"], "NOT FOUND"), "\n")
timepts <- c("0.5h","4h","24h","72h","7d","28d")
regions <- c("A","B","I")
sirt <- matrix(nrow=0,ncol=7)
colnames(sirt) <- c("region","time","SCI_n","sham_n","SCI_mean","sham_mean","delta_log2")
for (rg in regions) for (tp in timepts) {
  iSCI <- ann$region==rg & ann$treatment=="I" & ann$time==tp
  iSH  <- ann$region==rg & ann$treatment=="S" & ann$time==tp
  mSCI <- mean(geneMat["Sirt6", ann$GSM[iSCI]])
  mSH  <- mean(geneMat["Sirt6", ann$GSM[iSH]])
  sirt <- rbind(sirt, c(rg,tp,sum(iSCI),sum(iSH),mSCI,mSH,mSCI-mSH))
}
sirt <- as.data.frame(sirt, stringsAsFactors=FALSE)
for (cc in c("SCI_n","sham_n","SCI_mean","sham_mean","delta_log2")) sirt[[cc]] <- as.numeric(sirt[[cc]])
write.csv(sirt, "GSE5296_RMA_Sirt6_time_space.csv", row.names=FALSE)
cat("wrote GSE5296_RMA_Sirt6_time_space.csv\n")
print(sirt)

## ---------- 5. GO modules ----------
goMap <- list(
  OXPHOS = "GO:0006119",
  ComplexI_assembly = "GO:0032981",
  Aerobic_respiration = "GO:0009060",
  Mito_translation = "GO:0032543",
  Mito_morphogenesis = "NOT_FOUND_IN_GO"
)
## resolve GO->genes
moduleList <- list(); moduleStats <- list()
for (nm in names(goMap)) {
  gid <- goMap[[nm]]
  if (gid == "NOT_FOUND_IN_GO") { moduleStats[[nm]] <- c(go=gid, raw=NA, testable=NA, final=NA, note="no exact GO term (closest: GO:0007005 organization / GO:0000266 fission / GO:0008053 fusion); excluded per instruction"); next }
  map <- AnnotationDbi::select(org.Mm.eg.db, keys=gid, keytype="GO", columns="SYMBOL")
  syms <- unique(map$SYMBOL); syms <- syms[!is.na(syms) & syms!=""]
  testable <- syms[syms %in% rownames(geneMat)]
  moduleList[[nm]] <- testable
  moduleStats[[nm]] <- c(go=gid, raw=length(syms), testable=length(testable), final=length(testable), note="")
}
cat("=== module build stats ===\n")
for (nm in names(goMap)) { s <- moduleStats[[nm]]; cat(nm, ": ", paste(s,collapse=" "), "\n") }

## ---------- 6. scoring: ssGSEA + eigengene ----------
mods <- names(moduleList)[lengths(moduleList)>0]
## ssGSEA (GSVA 2.6.6 ssgseaParam)
gsets <- lapply(mods, function(m) moduleList[[m]]); names(gsets) <- mods
param <- GSVA::ssgseaParam(geneMat, gsets)
ss <- GSVA::gsva(param)
ss <- as.matrix(ss)  # modules x samples
cat("ssGSEA score dim:", nrow(ss), "x", ncol(ss), "\n")

## eigengene
eg <- matrix(NA, nrow=length(mods), ncol=ncol(geneMat), dimnames=list(mods, colnames(geneMat)))
for (m in mods) {
  X <- geneMat[moduleList[[m]], , drop=FALSE]
  Z <- t(scale(t(X)))            # standardize each gene across samples
  pca <- prcomp(t(Z))            # samples x genes
  v <- pca$x[,1]
  mu <- colMeans(X)
  if (cor(v, mu) < 0) v <- -v
  eg[m,] <- v
}
cat("eigengene done\n")

scoreList <- list(ssgsea=ss, eigengene=eg)

## ---------- 7. time x space comparison ----------
res <- NULL
for (mm in names(scoreList)) {
  sc <- scoreList[[mm]]
  for (m in mods) for (rg in regions) for (tp in timepts) {
    iSCI <- ann$region==rg & ann$treatment=="I" & ann$time==tp
    iSH  <- ann$region==rg & ann$treatment=="S" & ann$time==tp
    d <- mean(sc[m, ann$GSM[iSCI]]) - mean(sc[m, ann$GSM[iSH]])
    res <- rbind(res, data.frame(method=mm, module=m, region=rg, time=tp,
                 SCI_n=sum(iSCI), sham_n=sum(iSH),
                 SCI_mean_score=mean(sc[m, ann$GSM[iSCI]]),
                 sham_mean_score=mean(sc[m, ann$GSM[iSH]]),
                 delta_score=d, stringsAsFactors=FALSE))
  }
}
res <- as.data.frame(res, stringsAsFactors=FALSE)
write.csv(res, "GSE5296_RMA_mito_modules_time_space.csv", row.names=FALSE)
cat("wrote GSE5296_RMA_mito_modules_time_space.csv (", nrow(res), " rows)\n")

## ---------- 8. impact replication summary ----------
imp <- NULL
for (mm in names(scoreList)) {
  sc <- scoreList[[mm]]
  for (m in mods) {
    d <- sapply(timepts, function(tp){
      iSCI <- ann$region=="I" & ann$treatment=="I" & ann$time==tp
      iSH  <- ann$region=="I" & ann$treatment=="S" & ann$time==tp
      mean(sc[m, ann$GSM[iSCI]]) - mean(sc[m, ann$GSM[iSH]])
    })
    imp <- rbind(imp, data.frame(method=mm, module=m, t(d), number_negative_of_6=sum(d<0), stringsAsFactors=FALSE))
  }
}
colnames(imp)[3:8] <- timepts
write.csv(imp, "GSE5296_RMA_impact_replication_summary.csv", row.names=FALSE)
cat("wrote GSE5296_RMA_impact_replication_summary.csv\n")

## ---------- 9. GO Jaccard ----------
g2 <- mods[lengths(moduleList[mods])>0]
jac <- NULL
for (i in seq_along(g2)) for (j in seq_along(g2)) {
  A <- moduleList[[g2[i]]]; B <- moduleList[[g2[j]]]
  inter <- length(intersect(A,B)); uni <- length(union(A,B))
  jac <- rbind(jac, data.frame(moduleA=g2[i], moduleB=g2[j],
              shared_genes=inter, union_genes=uni,
              Jaccard=round(inter/uni,4), stringsAsFactors=FALSE))
}
write.csv(jac, "GSE5296_mito_GO_Jaccard.csv", row.names=FALSE)
cat("wrote GSE5296_mito_GO_Jaccard.csv\n")

## ---------- 10. QC sensitivity ----------
qc <- read.csv("GSE5296_RMA_QC_summary.csv", check.names=FALSE, stringsAsFactors=FALSE)
outl <- qc$GSM[qc$obvious_outlier==TRUE]
cat("pre-specified QC outliers (", length(outl), "):", paste(outl, collapse=" "), "\n")
if (length(outl)==0) {
  writeLines("No pre-specified outliers; sensitivity analysis not performed.",
             "GSE5296_RMA_QC_sensitivity.csv")
  cat("wrote GSE5296_RMA_QC_sensitivity.csv (no outliers)\n")
} else {
  sens <- NULL
  for (mm in names(scoreList)) {
    sc <- scoreList[[mm]]
    for (m in mods) for (tp in timepts) {
      iSCI <- ann$region=="I" & ann$treatment=="I" & ann$time==tp
      iSH  <- ann$region=="I" & ann$treatment=="S" & ann$time==tp
      full <- mean(sc[m, ann$GSM[iSCI]]) - mean(sc[m, ann$GSM[iSH]])
      # remove outliers present in the group
      iSCI2 <- iSCI & !(ann$GSM %in% outl)
      iSH2  <- iSH  & !(ann$GSM %in% outl)
      red <- mean(sc[m, ann$GSM[iSCI2]]) - mean(sc[m, ann$GSM[iSH2]])
      if (is.na(full) || is.na(red)) { dch <- NA }
      else if (sign(full)==0 || sign(red)==0) { dch <- (full*red) < 0 }
      else { dch <- sign(full) != sign(red) }
      sens <- rbind(sens, data.frame(excluded_GSM=paste(outl, collapse=";"), method=mm,
                   module=m, time=tp, full_delta=full, sensitivity_delta=red,
                   direction_changed=dch, stringsAsFactors=FALSE))
    }
  }
  write.csv(sens, "GSE5296_RMA_QC_sensitivity.csv", row.names=FALSE)
  cat("wrote GSE5296_RMA_QC_sensitivity.csv (", nrow(sens), " rows)\n")
}

cat("DONE", ts(), "\n")
sink(); close(logf)
cat("all complete\n")
