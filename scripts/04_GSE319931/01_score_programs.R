## GSE319931 chronic anatomical validation - Phase 3
## Predefined programs only: ComplexI_assembly, Aerobic_respiration
## Units: animal (paired, n=4 SCI; n=4 UI). No new pathways. No SIRT6.
suppressPackageStartupMessages({
  library(edgeR)
})
base <- if (length(commandArgs(trailingOnly=TRUE))) commandArgs(trailingOnly=TRUE)[[1]] else getwd()
repo_root <- Sys.getenv("SCI_REPOSITORY_ROOT", unset=getwd())
gse <- file.path(base,"GSE319931")
out <- gse
dir.create(out, showWarnings=FALSE)
LOG <- character()
note <- function(m){ cat(m,"\n"); LOG <<- c(LOG, m) }
note(paste("GSE319931 Phase3 started:", format(Sys.time(),"%Y-%m-%d %H:%M:%S")))

# ---- load counts & design ----
cnt <- read.delim(file.path(gse,"GSE319931_counts.tsv.gz"), header=TRUE, sep="\t", check.names=FALSE)
design <- read.csv(file.path(repo_root,"metadata","sample_design","GSE319931_sample_design.csv"), stringsAsFactors=FALSE)
samples <- design$sample
stopifnot(all(samples %in% colnames(cnt)))
counts <- as.matrix(cnt[, samples, drop=FALSE])
rownames(counts) <- cnt$GeneName
note(sprintf("counts matrix: %d genes x %d samples", nrow(counts), ncol(counts)))

# ---- target gene list (locked version) ----
tl <- read.csv(file.path(repo_root,"gene_programs","locked_respiratory_gene_programs.csv"), stringsAsFactors=FALSE)
ci_syms <- tl$GeneSymbol[ tl$ComplexI_assembly=="TRUE" ]
ae_syms <- tl$GeneSymbol[ tl$Aerobic_respiration=="TRUE" ]
note(sprintf("target: ComplexI=%d (orig), Aerobic=%d (orig)", length(ci_syms), length(ae_syms)))

map_prog <- function(syms, tag){
  present <- syms %in% rownames(counts)
  matched <- syms[present]; unmatched <- syms[!present]
  note(sprintf("%s: orig=%d matched_in_counts=%d unmatched=%d [%s]",
               tag, length(syms), sum(present), sum(!present), paste(unmatched, collapse=",")))
  list(idx=which(rownames(counts) %in% matched), matched=matched, unmatched=unmatched)
}
ci <- map_prog(ci_syms,"ComplexI_assembly")
ae <- map_prog(ae_syms,"Aerobic_respiration")

# ---- edgeR TMM -> logCPM ----
dge <- DGEList(counts)
keep <- filterByExpr(dge)
note(sprintf("filterByExpr keeps %d/%d genes", sum(keep), length(keep)))
dge_keep <- dge[keep,,keep.lib.sizes=FALSE]
dge_keep <- normLibSizes(dge_keep, method="TMM")
dge$samples$norm.factors <- dge_keep$samples$norm.factors
lcpm <- cpm(dge, log=TRUE, prior.count=2)   # logCPM full matrix
note(sprintf("logCPM built; lib.size norm.factors (TMM) applied; prior.count=2"))

# ---- program scoring ----
scoreA <- function(prog_idx){
  X <- lcpm[prog_idx,,drop=FALSE]
  z <- t(scale(t(X)))
  z[is.na(z)] <- 0
  colMeans(z)
}
scoreB <- function(prog_idx){
  X <- lcpm[prog_idx,,drop=FALSE]
  z <- t(scale(t(X)))
  v <- apply(z,1,sd); v[is.na(v)]<-0
  keep <- v>1e-9
  if(sum(keep)==0) return(rep(NA, ncol(X)))
  Z <- z[keep,,drop=FALSE]
  pc <- prcomp(t(Z), center=FALSE, scale.=FALSE)
  pc1 <- pc$x[,1]
  ref <- colMeans(Z)
  if(cor(pc1, ref) < 0) pc1 <- -pc1
  pc1
}
progs <- list(ComplexI_assembly=ci$idx, Aerobic_respiration=ae$idx)
scores <- list()
for(pn in names(progs)){
  scores[[paste(pn,"A")]] <- scoreA(progs[[pn]])
  scores[[paste(pn,"B")]] <- scoreB(progs[[pn]])
  note(sprintf("%s methodA range %.3f..%.3f ; methodB range %.3f..%.3f", pn,
        range(scores[[paste(pn,"A")]])[1], range(scores[[paste(pn,"A")]])[2],
        range(scores[[paste(pn,"B")]])[1], range(scores[[paste(pn,"B")]])[2]))
}

# long scores table
sc_long <- data.frame(sample=rep(samples, each=4),
                      program=rep(rep(names(progs),each=2), length(samples)),
                      method=rep(c("A","B"), length(samples)*length(progs)),
                      score=NA)
k<-1
for(s in samples){
  for(pn in names(progs)){
    for(m in c("A","B")){ sc_long$score[k] <- scores[[paste(pn,m)]][s]; k<-k+1 }
  }
}
sc_long$animal <- design$animal[match(sc_long$sample, design$sample)]
sc_long$condition <- design$condition[match(sc_long$sample, design$sample)]
sc_long$region <- design$region[match(sc_long$sample, design$sample)]
write.csv(sc_long, file.path(out,"GSE319931_program_scores.csv"), row.names=FALSE)
note("program_scores written")

# ---- primary paired regional (SCI within, n=4 animals) ----
reg_order <- c("Above","Epicenter","Below")
paired_res <- list()
for(pn in names(progs)){
  for(m in c("A","B")){
    key <- paste(pn,m)
    val <- scores[[key]]
    per_animal <- sapply(1:4, function(a){
      g <- design$sample[design$animal==a & design$condition=="SCI"]
      setNames(val[g], design$region[design$animal==a & design$condition=="SCI"])
    })  # 3 regions x 4 animals
    rownames(per_animal) <- reg_order
    # paired differences
    for(pair in c("E-A","E-B","A-B")){
      rr <- strsplit(pair,"-")[[1]]
      r1 <- ifelse(rr[1]=="E","Epicenter", ifelse(rr[1]=="A","Above","Below"))
      r2 <- ifelse(rr[2]=="E","Epicenter", ifelse(rr[2]=="A","Above","Below"))
      d <- per_animal[r1,] - per_animal[r2,]
      t <- tryCatch(t.test(d), error=function(e) NULL)
      w <- tryCatch(wilcox.test(d, exact=FALSE), error=function(e) NULL)
      paired_res[[length(paired_res)+1]] <- data.frame(
        program=pn, method=m, comparison=pair, region1=r1, region2=r2,
        mean_diff=mean(d), median_diff=median(d),
        n_consistent_neg=sum(d<0), n_total=4,
        cohen_dz=ifelse(sd(d)>0, mean(d)/sd(d), NA),
        paired_t_p=ifelse(is.null(t), NA, t$p.value),
        wilcox_p=ifelse(is.null(w), NA, w$p.value))
    }
  }
}
paired_df <- do.call(rbind, paired_res)
write.csv(paired_df, file.path(out,"GSE319931_paired_regional.csv"), row.names=FALSE)
note("paired regional written")

# ---- secondary: SCI region vs UI (descriptive/exploratory) ----
sec_res <- list()
ui_idx <- design$sample[design$condition=="Uninjured"]
for(pn in names(progs)){
  for(m in c("A","B")){
    key<-paste(pn,m); val<-scores[[key]]
    for(rg in c("Above","Epicenter","Below")){
      inj_idx <- design$sample[design$condition=="SCI" & design$region==rg]
      a<-val[ui_idx]; b<-val[inj_idx]
      t<-tryCatch(t.test(b,a), error=function(e) NULL)
      sec_res[[length(sec_res)+1]] <- data.frame(
        program=pn, method=m, region=rg,
        UI_mean=mean(a), SCI_mean=mean(b), delta=mean(b)-mean(a),
        direction=ifelse(mean(b)-mean(a)<0,"down","up"),
        n_inj_below_UImedian=sum(b < median(a)), n_total=4,
        welch_p=ifelse(is.null(t), NA, t$p.value))
    }
  }
}
sec_df <- do.call(rbind, sec_res)
write.csv(sec_df, file.path(out,"GSE319931_secondary_vs_UI.csv"), row.names=FALSE)
note("secondary vs UI written")

# ---- per-animal line figure (Above -> Epicenter -> Below) ----
library(grDevices)
pdf(file.path(out,"GSE319931_per_animal_regional.pdf"), width=8, height=8)
par(mfrow=c(2,2), mar=c(4,4,2,1))
for(pn in names(progs)){
  for(m in c("A","B")){
    key<-paste(pn,m); val<-scores[[key]]
    X <- sapply(1:4, function(a){
      g<-design$sample[design$animal==a & design$condition=="SCI"]
      val[g][order(match(design$region[design$animal==a & design$condition=="SCI"], reg_order))]
    })
    matplot(1:3, X, type="b", lty=1, pch=19, col=c("black","red","blue","green3"),
            xaxt="n", ylab="score", main=sprintf("%s / method %s", pn, m))
    axis(1, at=1:3, labels=reg_order)
    ui_mean <- mean(val[ui_idx])
    abline(h=ui_mean, lty=2, col="gray60")
  }
}
dev.off()
note("per-animal figure written")

# ---- session info / log ----
si <- capture.output(sessionInfo())
writeLines(c("GSE319931 Phase3 environment", si), file.path(out,"GSE319931_sessionInfo.txt"))
writeLines(c("GSE319931 Phase3 log", LOG, paste("finished:", format(Sys.time(),"%Y-%m-%d %H:%M:%S"))),
           file.path(out,"GSE319931_analysis_log.txt"))
note("DONE")
