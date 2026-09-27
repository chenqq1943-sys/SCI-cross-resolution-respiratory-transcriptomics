# =====================================================================
# Phase 2C CORRECTED: full-transcriptome pseudobulk + edgeR TMM + logCPM
# GSE234774 (Tabulae Paralytica) snRNA, library x celltype
# Full-transcriptome library sizes are used for TMM/logCPM; locked
# respiratory target genes are extracted only AFTER normalization.
# =====================================================================
suppressMessages({ library(edgeR) })
base <- if (length(commandArgs(trailingOnly=TRUE))) commandArgs(trailingOnly=TRUE)[[1]] else getwd()
repo_root <- Sys.getenv("SCI_REPOSITORY_ROOT", unset=getwd())
outdir   <- file.path(base,"Phase2C_corrected_fulltranscriptome_pseudobulk")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
t_start <- Sys.time()
LOG <- character()
note <- function(...) { m <- paste0(collapse=" ", list(...)); cat(m, "\n"); LOG <<- c(LOG, m) }
note("Phase 2C CORRECTED started:", format(t_start))

# ---------------- 0. Inputs ----------------
tl     <- read.csv(file.path(repo_root,"gene_programs","locked_respiratory_gene_programs.csv"), stringsAsFactors=FALSE)
counts <- read.csv(gzfile(file.path(outdir,"fulltranscriptome_pseudobulk_counts.csv.gz")),
                   check.names=FALSE, stringsAsFactors=FALSE)
cold   <- read.csv(file.path(outdir,"fulltranscriptome_pseudobulk_coldata.csv"), stringsAsFactors=FALSE)

genes <- counts[[1]]
M     <- as.matrix(counts[,-1]); rownames(M) <- genes
samp  <- colnames(M)
note("FULL pseudobulk counts dims: ", nrow(M), " genes x ", ncol(M), " samples")
stopifnot(identical(samp, cold$pseudobulk_sample))
note("coldata order 100% matches counts columns: TRUE")
stopifnot(all(cold$full_transcriptome_total_UMI == cold$meta_nCount_RNA_sum))
note("full_transcriptome_total_UMI == meta_nCount_RNA_sum (all samples): PASS")

# ---------------- 1. Library-size verification ----------------
cs <- colSums(M)
ok <- all(cs == cold$full_transcriptome_total_UMI)
note("colSums(full pseudobulk) == full_transcriptome_total_UMI == nCount_RNA sum: ", ok,
     " ; max|diff|=", max(abs(cs - cold$full_transcriptome_total_UMI)))
libsize_QC <- data.frame(
  check=c("colSums==cold$full_total_UMI","cold$full==meta_nCount","n_samples","n_genes",
          "77 samples documented (80-3 missing Microglia 1d_2/1d_8/4d_6)"),
  status=c(ifelse(ok,"PASS","FAIL"),
           ifelse(all(cold$full_transcriptome_total_UMI==cold$meta_nCount_RNA_sum),"PASS","FAIL"),
           as.character(ncol(M)), as.character(nrow(M)), "OK"),
  detail=c(sprintf("max|diff|=%d", max(abs(cs-cold$full_transcriptome_total_UMI))),"",
           sprintf("%d samples", ncol(M)), sprintf("%d genes", nrow(M)),
           "77 = 20 libs x 4 celltypes (80) - 3 absent Microglia combos"))
write.csv(libsize_QC, file.path(outdir,"normalization_QC.csv"), row.names=FALSE)
note("normalization_QC written")

# ---------------- 2. edgeR DGEList + filterByExpr on FULL transcriptome ----------------
dge <- DGEList(counts=M)
design_group <- factor(cold$celltype)
keep <- filterByExpr(dge, group=design_group)
note("filterByExpr: genes before=", nrow(dge), " after=", sum(keep))
dge <- dge[keep, , keep.lib.sizes=TRUE]
dge <- normLibSizes(dge, method="TMM")
lcpm <- cpm(dge, log=TRUE, prior.count=3)
colnames(lcpm) <- samp
note("TMM + logCPM (prior.count=3) on FULL transcriptome. edgeR=",
     as.character(packageVersion("edgeR")), " limma=", as.character(packageVersion("limma")))
write.csv(cbind(GeneSymbol=rownames(lcpm), as.data.frame(lcpm)),
          gzfile(file.path(outdir,"corrected_pseudobulk_logCPM.csv.gz")), row.names=FALSE)
note("corrected logCPM written: ", nrow(lcpm), " x ", ncol(lcpm))

# ---------------- 3. Extract locked target genes AFTER normalization ----------------
complexi_genes <- tl$GeneSymbol[tl$ComplexI_assembly==TRUE]
aerobic_genes  <- tl$GeneSymbol[tl$Aerobic_respiration==TRUE]
sirt6_gene     <- tl$GeneSymbol[tl$Sirt6==TRUE]
stopifnot(length(complexi_genes)==30, length(aerobic_genes)==53)
ci_av <- intersect(complexi_genes, rownames(lcpm)); ci_na <- setdiff(complexi_genes, rownames(lcpm))
ae_av <- intersect(aerobic_genes, rownames(lcpm)); ae_na <- setdiff(aerobic_genes, rownames(lcpm))
s6_av <- intersect(sirt6_gene, rownames(lcpm))
# mtDNA check
all_av <- rownames(lcpm)
mt_hits <- all_av[grepl("^(mt-|Mt-|MT-)", all_av)]
note("target availability: ComplexI=", length(ci_av), "/30 ; Aerobic=", length(ae_av), "/53 ; Sirt6 in lcpm=", s6_av %in% rownames(lcpm))
if (length(ci_na)>0) note("  ComplexI unavailable: ", paste(ci_na,collapse=","))
if (length(ae_na)>0) note("  Aerobic unavailable: ", paste(ae_na,collapse=","))
note("mtDNA-encoded genes present in lcpm: ", sum(grepl("^(mt-|Mt-|MT-)", ci_av)),
     " in ComplexI ; ", sum(grepl("^(mt-|Mt-|MT-)", ae_av)), " in Aerobic ; total mt-* in lcpm=", length(mt_hits))
target_map <- data.frame(
  program=c(rep("ComplexI_assembly",length(ci_av)), rep("Aerobic_respiration",length(ae_av))),
  GeneSymbol=c(ci_av, ae_av), stringsAsFactors=FALSE)
write.csv(target_map, file.path(outdir,"corrected_target_gene_availability.csv"), row.names=FALSE)

# ---------------- 4. Program scoring (EXACT original definitions) ----------------
ctypes <- c("Neurons","Astrocytes","Oligodendrocyte_combined","Microglia")
times_level <- c("uninjured","1d","4d","7d","14d","1m","2m")
cold$time <- factor(cold$time, levels=times_level)
cold$library <- as.character(cold$library)
progs <- list(ComplexI_assembly=ci_av, Aerobic_respiration=ae_av)

score_one <- function(ct, prog, keep_libs=NULL){
  libs <- samp[cold$celltype==ct]
  if (!is.null(keep_libs)) libs <- intersect(libs, keep_libs)
  P <- lcpm[prog, libs, drop=FALSE]
  gsd <- apply(P,1,sd)
  usable <- !is.na(gsd) & gsd>0
  P2 <- P[usable,,drop=FALSE]
  Z <- scale(t(P2))
  A <- rowMeans(Z, na.rm=TRUE); names(A) <- colnames(P2)
  pc <- prcomp(t(P2), center=TRUE, scale.=TRUE)
  pc1 <- pc$x[,1]; names(pc1) <- colnames(P2)
  if (cor(pc1, colMeans(P2), use="pairwise") < 0) pc1 <- -pc1
  list(A=A, B=pc1, n_gene=length(prog), n_usable=nrow(P2), libs=colnames(P2))
}
build_scores <- function(keep=NULL){
  rows <- list()
  for (ct in ctypes) for (pn in names(progs)){
    sc <- score_one(ct, progs[[pn]], keep)
    for (met in c("A","B")) for (lib in names(sc[[met]])){
      rows[[length(rows)+1]] <- data.frame(
        pseudobulk_sample=lib, celltype=ct,
        library=cold$library[match(lib,cold$pseudobulk_sample)],
        time=as.character(cold$time[match(lib,cold$pseudobulk_sample)]),
        program=pn, method=met, score=as.numeric(sc[[met]][[lib]]),
        n_program_genes=sc$n_gene, n_usable_genes=sc$n_usable, stringsAsFactors=FALSE)
    }
  }
  do.call(rbind, rows)
}
scores_main <- build_scores()
write.csv(scores_main, file.path(outdir,"corrected_program_scores.csv"), row.names=FALSE)
note("corrected program_scores (main) written: ", nrow(scores_main), " rows")

# ---------------- 5. Comparisons (EXACT original) ----------------
cohend <- function(x,y){ nx<-length(x); ny<-length(y)
  s2 <- ((nx-1)*var(x)+(ny-1)*var(y))/(nx+ny-2)
  if (is.na(s2)||s2<=0) return(NA_real_); (mean(y)-mean(x))/sqrt(s2) }
compare <- function(scores_df, ct, pn, met, timeT){
  ss <- scores_df[scores_df$celltype==ct & scores_df$program==pn & scores_df$method==met,]
  ctrl_libs <- samp[cold$celltype==ct & cold$time=="uninjured"]
  inj_libs  <- samp[cold$celltype==ct & cold$time==timeT]
  cs <- ss$score[match(intersect(ctrl_libs,ss$pseudobulk_sample), ss$pseudobulk_sample)]
  is_ <- ss$score[match(intersect(inj_libs,ss$pseudobulk_sample), ss$pseudobulk_sample)]
  cs <- cs[!is.na(cs)]; is_ <- is_[!is.na(is_)]
  nC<-length(cs); nI<-length(is_)
  if (nC<2 || nI<2) return(NULL)
  mC<-mean(cs); mI<-mean(is_); dlt<-mI-mC
  dir <- ifelse(dlt<0,"down","up")
  es <- cohend(cs,is_)
  pnom <- tryCatch(t.test(is_,cs)$p.value, error=function(e) NA_real_)
  md <- median(cs); nbelow<-sum(is_<md); nabove<-sum(is_>md)
  cons <- sprintf("%d/%d", ifelse(dir=="down",nbelow,nabove), nI)
  data.frame(celltype=ct, program=pn, method=met, time=timeT, n_control=nC, n_injured=nI,
             control_mean=mC, injured_mean=mI, delta=dlt, direction=dir,
             effect_size=es, nominal_p=pnom, replicate_consistency=cons, stringsAsFactors=FALSE)
}
run_prim <- function(scores_df){
  rows <- list()
  for (ct in ctypes) for (pn in names(progs)) for (met in c("A","B"))
    for (tt in c("7d","14d","1m","2m")){
      r <- compare(scores_df, ct, pn, met, tt); if(!is.null(r)) rows[[length(rows)+1]] <- r
    }
  tab <- do.call(rbind, rows); tab$adjusted_p <- p.adjust(tab$nominal_p, method="BH"); tab
}
prim_main <- run_prim(scores_main)
write.csv(prim_main, file.path(outdir,"corrected_primary_results.csv"), row.names=FALSE)
note("corrected primary comparisons written: ", nrow(prim_main), " rows")

# ---------------- 6. Exploratory 1d/4d ----------------
expl_full_cols <- c("celltype","program","method","time","n_control","n_injured","control_mean",
                    "injured_mean","delta","direction","effect_size","nominal_p",
                    "replicate_consistency","status","note")
expl <- list()
for (ct in ctypes) for (pn in names(progs)) for (met in c("A","B")) for (tt in c("1d","4d")){
  ss <- scores_main[scores_main$celltype==ct & scores_main$program==pn & scores_main$method==met,]
  ctrl_libs <- samp[cold$celltype==ct & cold$time=="uninjured"]
  inj_libs  <- samp[cold$celltype==ct & cold$time==tt]
  nI <- sum(inj_libs %in% ss$pseudobulk_sample)
  if (nI < 2){
    d <- data.frame(celltype=ct, program=pn, method=met, time=tt, n_control=length(ctrl_libs),
                    n_injured=nI, control_mean=NA, injured_mean=NA, delta=NA, direction=NA,
                    effect_size=NA, nominal_p=NA, replicate_consistency=NA,
                    status="Not evaluable", note="fewer than 2 injured libraries with score",
                    stringsAsFactors=FALSE)
  } else {
    r <- compare(scores_main, ct, pn, met, tt); r$status<-"Evaluable"; r$note<-""
    d <- r
  }
  expl[[length(expl)+1]] <- d[, expl_full_cols]
}
expl_tab <- do.call(rbind, expl); expl_tab$adjusted_p <- NA_real_
write.csv(expl_tab, file.path(outdir,"corrected_exploratory_1d_4d.csv"), row.names=FALSE)
note("corrected exploratory 1d/4d written")

# ---------------- 7. Sensitivity (SAME 7 prespecified combos) ----------------
low_samples <- c("timecourse_1d_3__Microglia","timecourse_4d_5__Microglia",
                 "uninjured_uninjured_10__Microglia","timecourse_1d_8__Oligodendrocyte_combined",
                 "timecourse_4d_6__Oligodendrocyte_combined","timecourse_1d_3__Neurons",
                 "timecourse_4d_6__Neurons")
note("sensitivity: excluding same 7 prespecified combos: ", paste(low_samples, collapse=", "))
prim_sens <- NULL
keep_all <- samp[!(samp %in% low_samples)]
scores_sens <- build_scores(keep_all)
prim_sens <- run_prim(scores_sens)
sens_rows <- list()
for (ct in ctypes) for (pn in names(progs)) for (met in c("A","B")) for (tt in c("7d","14d","1m","2m")){
  rm_ <- prim_main[prim_main$celltype==ct & prim_main$program==pn & prim_main$method==met & prim_main$time==tt,]
  rs <- prim_sens[prim_sens$celltype==ct & prim_sens$program==pn & prim_sens$method==met & prim_sens$time==tt,]
  dir_chg <- NA; mag_chg <- NA; cons_chg <- NA
  if (nrow(rm_)==1 && nrow(rs)==1){
    dir_chg <- rm_$direction != rs$direction
    mag_chg <- !is.na(rm_$delta) && !is.na(rs$delta) && abs(rs$delta - rm_$delta) > 0.15*max(abs(rm_$delta),0.05)
    cons_chg <- rm_$replicate_consistency != rs$replicate_consistency
  }
  sens_rows[[length(sens_rows)+1]] <- data.frame(
    celltype=ct, program=pn, method=met, time=tt,
    main_direction=if(nrow(rm_)==1) rm_$direction else NA,
    sens_direction=if(nrow(rs)==1) rs$direction else NA,
    main_delta=if(nrow(rm_)==1) rm_$delta else NA,
    sens_delta=if(nrow(rs)==1) rs$delta else NA,
    main_consistency=if(nrow(rm_)==1) rm_$replicate_consistency else NA,
    sens_consistency=if(nrow(rs)==1) rs$replicate_consistency else NA,
    direction_changed=dir_chg, effect_magnitude_changed=mag_chg,
    replicate_consistency_changed=cons_chg, stringsAsFactors=FALSE)
}
sens_tab <- do.call(rbind, sens_rows)
write.csv(sens_tab, file.path(outdir,"corrected_sensitivity.csv"), row.names=FALSE)
note("corrected sensitivity written")

# ---------------- 9. Astrocyte Aerobic gene-driver audit ----------------
if (length(ae_av)>=1){
  ast_libs <- samp[cold$celltype=="Astrocytes"]
  uni_libs <- ast_libs[cold$time[match(ast_libs,samp)]=="uninjured"]
  P <- lcpm[ae_av, ast_libs, drop=FALSE]
  Z <- scale(t(P))                     # libs x genes, gene-scaled across astro libs (Method A scale)
  timesU <- c("7d","14d","1m","2m")
  rows <- list(); per_time <- list()
  for (tt in timesU){
    inj_libs <- ast_libs[cold$time[match(ast_libs,samp)]==tt]
    contrib <- colMeans(Z[inj_libs,,drop=FALSE]) - colMeans(Z[uni_libs,,drop=FALSE])  # per gene z contribution
    for (g in ae_av){
      rows[[length(rows)+1]] <- data.frame(gene=g, time=tt,
        uninjured_mean_logCPM=mean(P[g,uni_libs]), injured_mean_logCPM=mean(P[g,inj_libs]),
        delta_logCPM=mean(P[g,inj_libs])-mean(P[g,uni_libs]),
        direction=ifelse(mean(P[g,inj_libs])>=mean(P[g,uni_libs]),"up","down"),
        contribution_z=as.numeric(contrib[g]), stringsAsFactors=FALSE)
    }
    per_time[[tt]] <- data.frame(gene=ae_av, contribution_z=contrib[ae_av], stringsAsFactors=FALSE)
  }
  audit <- do.call(rbind, rows)
  write.csv(audit, file.path(outdir,"astrocyte_aerobic_gene_driver_audit.csv"), row.names=FALSE)
  sums <- do.call(rbind, lapply(timesU, function(tt){
    s <- audit[audit$time==tt,]
    data.frame(time=tt, n_genes=nrow(s), n_up=sum(s$direction=="up"),
               n_down=sum(s$direction=="down"),
               frac_up=round(sum(s$direction=="up")/nrow(s),3), stringsAsFactors=FALSE)
  }))
  write.csv(sums, file.path(outdir,"astrocyte_aerobic_gene_driver_summary.csv"), row.names=FALSE)
  for (tt in timesU){
    g <- per_time[[tt]]
    top <- head(g[order(-g$contribution_z),],5); bot <- head(g[order(g$contribution_z),],5)
    note(sprintf("Astro/Aerobic %s: top+ drivers: %s ; top- drivers: %s",
                 tt,
                 paste(sprintf("%s(z=%.2f)",top$gene,top$contribution_z),collapse=","),
                 paste(sprintf("%s(z=%.2f)",bot$gene,bot$contribution_z),collapse=",")))
  }
  for (tt in timesU) note(sprintf("Astro/Aerobic %s: n_up=%d n_down=%d frac_up=%.3f",
                                  tt, sums$n_up[sums$time==tt], sums$n_down[sums$time==tt], sums$frac_up[sums$time==tt]))
  for (g in c("Mup4","Mup5")){
    if (g %in% rownames(lcpm)){
      note(paste0("Mup",g," present in lcpm; astro uninjured mean logCPM=",
                 round(mean(lcpm[g, uni_libs]),3)))
    } else note(paste0("Mup",g," NOT available in corrected logCPM (filtered out or absent)"))
  }
} else {
  note("No aerobic genes available in lcpm; astrocyte audit skipped.")
}

# ---------------- 10. sessionInfo & log ----------------
writeLines(capture.output(sessionInfo()), file.path(outdir,"sessionInfo.txt"))
writeLines(c(paste0("Analysis log — Phase 2C CORRECTED (full-transcriptome normalization)"),
             paste0("Started: ", format(t_start)),
             LOG,
             paste0("Finished: ", format(Sys.time()))),
           file.path(outdir,"analysis_log.txt"))
writeLines(c("Phase 2C corrected analysis script",
             "Rebuild: _rebuild_full_pseudobulk.py  ;  Analysis: _phase2c_corrected.R"),
           file.path(outdir,"_script_manifest.txt"))
note("DONE in ", round(difftime(Sys.time(), t_start, units="mins"),1), " min")
