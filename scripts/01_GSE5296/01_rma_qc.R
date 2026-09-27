# GSE5296 预处理：CEL完整性 -> QC -> 标准RMA -> 导出。仅做这些, 不做任何生物学分析。
# Use the caller's existing R library paths.
options(repos = c(CRAN = "https://cloud.r-project.org"))
suppressPackageStartupMessages({library(affy); library(affyPLM)})

base <- if (length(commandArgs(trailingOnly=TRUE))) commandArgs(trailingOnly=TRUE)[[1]] else getwd()
setwd(base)
raw_dir <- file.path(base, "GSE5296_RAW")
ann_path <- file.path(base, "GSE5296_sample_annotation.csv")

out <- file("_pipeline_log.txt", open = "wt")
sink(out, type = "output", split = TRUE)
sink(out, type = "message", append = TRUE)
on.exit({sink(type="message"); sink(type="output"); close(out)})

cat("START", format(Sys.time()), "\n")
cat("R:", R.version.string, "\n")

## ---------- 三、CEL 完整性检查 ----------
cel_files <- sort(list.files(raw_dir, pattern = "CEL\\.gz$", ignore.case = TRUE))
cat("CEL files found:", length(cel_files), "\n")
stopifnot(length(cel_files) == 96)

gsm_from <- function(f) sub("^((GSM\\d+)).*", "\\2", f, perl=TRUE)
cel_gsm <- sapply(cel_files, gsm_from)
names(cel_gsm) <- NULL

ann <- read.csv(ann_path, stringsAsFactors = FALSE, check.names = FALSE)
ann_gsm <- ann$sample_id
if (is.null(ann_gsm)) ann_gsm <- ann$GSM
stopifnot(length(ann_gsm) == 96)

dup <- cel_gsm[duplicated(cel_gsm)]
missing <- setdiff(ann_gsm, cel_gsm)
extra <- setdiff(cel_gsm, ann_gsm)

ann_map <- setNames(as.list(ann$time_after_injury), ann_gsm)   # placeholder, replaced below
# build annotation by GSM
ann_by <- split(ann, ann_gsm)
region_of <- function(g) ann_by[[g]]$region[1]
treat_of <- function(g) ann_by[[g]]$treatment[1]
time_of  <- function(g) ann_by[[g]]$time_after_injury[1]

integrity <- data.frame(
  GSM = cel_gsm,
  CEL_filename = cel_files,
  region = sapply(cel_gsm, region_of),
  treatment = sapply(cel_gsm, treat_of),
  time = sapply(cel_gsm, time_of),
  match_status = "MATCH",
  stringsAsFactors = FALSE
)
integrity$match_status[!integrity$GSM %in% ann_gsm] <- "NOT_IN_ANNOTATION"
integrity$match_status[duplicated(integrity$GSM)] <- "DUPLICATE"

write.csv(integrity, "GSE5296_CEL_integrity_check.csv", row.names = FALSE)

cat("duplicates:", length(dup), if(length(dup)) paste(dup, collapse=",") else "", "\n")
cat("missing (annotation not in CEL):", length(missing), if(length(missing)) paste(missing,collapse=",") else "", "\n")
cat("extra (CEL not in annotation):", length(extra), if(length(extra)) paste(extra,collapse=",") else "", "\n")

if (length(dup) > 0 || length(missing) > 0 || length(extra) > 0 || length(cel_files) != 96) {
  cat("INTEGRITY FAIL - STOP, no RMA\n")
  q("no", status = 1)
}
cat("INTEGRITY OK: 96/96 matched\n")

## ---------- 四、读取 CEL ----------
cat("Reading CEL files ...", format(Sys.time()), "\n")
affybatch <- ReadAffy(filenames = cel_files, celfile.path = raw_dir)
cat("AffyBatch samples:", length(sampleNames(affybatch)), "platform:", annotation(affybatch), "\n")
stopifnot(length(sampleNames(affybatch)) == 96)

# 记录 sessionInfo
si_file <- file("GSE5296_RMA_sessionInfo.txt", open = "wt")
sink(si_file, type = "output"); sink(si_file, type = "message", append = TRUE)
cat("GSE5296 RMA session info\n")
cat("run date:", format(Sys.time()), "\n")
cat("CEL count:", length(sampleNames(affybatch)), "\n")
cat("platform:", annotation(affybatch), "\n")
cat("affy:", as.character(packageVersion("affy")), "\n")
cat("affyPLM:", as.character(packageVersion("affyPLM")), "\n")
cat("Bioconductor:", as.character(BiocManager::version()), "\n")
cat("== sessionInfo ==\n")
print(sessionInfo())
sink(type = "message"); sink(type = "output"); close(si_file)

## ---------- 五、CEL级QC ----------
cat("QC ...", format(Sys.time()), "\n")
pdf("GSE5296_RMA_QC.pdf", width = 14, height = 10)
par(mar = c(7,4,2,1))
# 1) raw log2 intensity boxplot
boxplot(affybatch, main = "Raw log2 intensity distribution", col = rainbow(96), las = 2,
        names = gsub("\\.CEL\\.gz$", "", sampleNames(affybatch)))
# 2) RLE / NUSE (affyPLM)
cat("fitPLM ...", format(Sys.time()), "\n")
Pset <- fitPLM(affybatch, normalize = TRUE, background = TRUE)
RLE <- affyPLM::RLE(Pset, type = "stats")
NUSE <- affyPLM::NUSE(Pset, type = "stats")
par(mfrow = c(2,1))
boxplot(data.frame(RLE), main = "RLE (relative log expression)", las = 2, col = rainbow(96),
        names = gsub("\\.CEL\\.gz$", "", sampleNames(affybatch)))
boxplot(data.frame(NUSE), main = "NUSE (normalized unscaled SE)", las = 2, col = rainbow(96),
        names = gsub("\\.CEL\\.gz$", "", sampleNames(affybatch)))
dev.off()  # close first page set; we'll reopen
pdf("GSE5296_RMA_QC2.pdf", width = 12, height = 10)  # correlation + PCA go into separate pdf, renamed after
par(mfrow = c(1,1))

## ---------- 六、RMA ----------
cat("RMA ...", format(Sys.time()), "\n")
eset <- rma(affybatch)   # 标准 RMA: background correction + quantile norm + median polish
cat("RMA dim:", nrow(exprs(eset)), "x", ncol(exprs(eset)), "\n")

## 7) array-array correlation (RMA log2)
cm <- cor(exprs(eset))
heatmap(cm, symm = TRUE, distfun = function(x) as.dist(1 - x), main = "Array-array correlation (RMA log2)")

## 8) PCA (技术离群筛查, 不做生物学解释)
pc <- prcomp(t(exprs(eset)), center = TRUE, scale. = TRUE)
pvar <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
snames <- gsub("\\.CEL\\.gz$", "", sampleNames(eset))
plot(pc$x[,1], pc$x[,2], xlab = paste0("PC1 (", pvar[1], "%)"), ylab = paste0("PC2 (", pvar[2], "%)"),
     main = "PCA of RMA log2 (technical outlier screen)", pch = 16, col = 1)
text(pc$x[,1], pc$x[,2], labels = snames, pos = 3, cex = 0.5)
dev.off()

## ---------- QC summary ----------
RLE_stats <- affyPLM::RLE(Pset, type = "stats")   # matrix rows=c(median,MAD), cols=arrays
NUSE_stats <- affyPLM::NUSE(Pset, type = "stats") # matrix rows=c(median,IQR), cols=arrays
cat("RLE stats rownames:", rownames(RLE_stats), "\n")
cat("NUSE stats rownames:", rownames(NUSE_stats), "\n")
rle_med <- as.numeric(RLE_stats["median", ])
rle_iqr <- as.numeric(RLE_stats["IQR", ])
nuse_med <- as.numeric(NUSE_stats["median", ])
nuse_iqr <- as.numeric(NUSE_stats["IQR", ])
# 阈值: NUSE 中位数 > 1.05 或 RLE 中位数绝对值极端 -> 标记
nuse_flag <- nuse_med > 1.05
rle_global <- median(rle_med)
rle_flag <- abs(rle_med) > max(0.1, 2 * mad(rle_med))
outlier <- nuse_flag | rle_flag
remark <- ifelse(outlier, "POSSIBLE_TECHNICAL_OUTLIER (NUSE/RLE)", "OK")

# PCA 离群(简单: 距中心大)
cent <- colMeans(pc$x[,1:2]); d <- sqrt(rowSums((pc$x[,1:2] - cent)^2)); pc_flag <- d > mean(d) + 2*sd(d)
remark[pc_flag & !outlier] <- paste(remark[pc_flag & !outlier], "; PCA_DIST_OUTLIER")
outlier_final <- outlier | pc_flag

qcd <- data.frame(
  GSM = cel_gsm,
  CEL_filename = cel_files,
  region = sapply(cel_gsm, region_of),
  treatment = sapply(cel_gsm, treat_of),
  time = sapply(cel_gsm, time_of),
  RLE_median = round(rle_med, 4),
  RLE_IQR = round(rle_iqr, 4),
  NUSE_median = round(nuse_med, 4),
  NUSE_IQR = round(nuse_iqr, 4),
  obvious_outlier = outlier_final,
  remark = remark,
  stringsAsFactors = FALSE
)
# 统一排序到矩阵列序(与 sampleNames 一致)
qcd <- qcd[match(sampleNames(eset), cel_files), , drop = FALSE]
# 修正 GSM 顺序: sampleNames 即 cel_files 顺序
qcd$GSM <- gsm_from(sampleNames(eset)); qcd$CEL_filename <- sampleNames(eset)
qcd$region <- sapply(qcd$GSM, region_of); qcd$treatment <- sapply(qcd$GSM, treat_of); qcd$time <- sapply(qcd$GSM, time_of)

write.csv(qcd, "GSE5296_QC_summary.csv", row.names = FALSE)
write.csv(qcd, "GSE5296_RMA_QC_summary.csv", row.names = FALSE)

## ---------- 七、导出核心结果 ----------
# 1) 表达矩阵 (矩阵列序 = sampleNames)
expr <- exprs(eset)
colnames(expr) <- gsm_from(sampleNames(eset))
gexpr <- cbind(Probe_ID = rownames(expr), data.frame(expr, check.names = FALSE))
gz <- gzfile("GSE5296_RMA_probe_expression.csv.gz", "wt")
write.csv(gexpr, gz, row.names = FALSE)
close(gz)

# 2) 样本注释 (严格按矩阵96列顺序)
sample_ann <- data.frame(
  GSM = colnames(expr),
  CEL_filename = sampleNames(eset),
  region = sapply(colnames(expr), region_of),
  treatment = sapply(colnames(expr), treat_of),
  time = sapply(colnames(expr), time_of),
  title = sapply(colnames(expr), function(g) ann_by[[g]]$title[1]),
  stringsAsFactors = FALSE
)
write.csv(sample_ann, "GSE5296_RMA_sample_annotation.csv", row.names = FALSE)

cat("RMA matrix dim:", nrow(expr), "x", ncol(expr), "\n")
cat("probesets:", nrow(expr), "samples:", ncol(expr), "\n")
cat("outliers flagged:", sum(qcd$obvious_outlier), "\n")
cat("DONE", format(Sys.time()), "\n")
