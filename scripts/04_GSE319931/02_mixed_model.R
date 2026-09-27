# GSE319931: four prespecified regional models. No score or sample is changed.
# Usage: Rscript 02_mixed_model.R [data_dir] [output_dir] [repository_root]
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args) >= 1L) args[[1L]] else file.path(getwd(), "derived_data", "GSE319931")
output_dir <- if (length(args) >= 2L) args[[2L]] else file.path(getwd(), "reproduction_outputs", "GSE319931")
repo_root <- if (length(args) >= 3L) args[[3L]] else getwd()
if (!dir.exists(data_dir) || !dir.exists(output_dir)) stop("Data/output directory absent")

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

read_input <- function(name) read.csv(file.path(data_dir, name), stringsAsFactors = FALSE)
scores <- read_input("GSE319931_program_scores.csv")
design <- read.csv(file.path(repo_root, "metadata", "sample_design", "GSE319931_sample_design.csv"), stringsAsFactors = FALSE)
paired <- read_input("GSE319931_paired_regional.csv")
stopifnot(!anyDuplicated(design$sample), !anyDuplicated(scores[c("sample", "program", "method")]))
ix <- match(scores$sample, design$sample)
if (anyNA(ix)) stop("Program score contains a sample absent from design")
for (field in c("animal", "condition", "region")) {
  left <- as.character(scores[[field]])
  right <- as.character(design[[field]][ix])
  if (!identical(left, right)) stop(sprintf("Score/design disagreement: %s", field))
}

programs <- c("ComplexI_assembly", "Aerobic_respiration")
methods <- c("A", "B")
regions <- c("Above", "Epicenter", "Below")
contrasts <- list("E-A" = c(-1, 1, 0), "E-B" = c(0, 1, -1),
                  "A-B" = c(1, 0, -1))
key <- function(d) paste(d$program, d$method, sep = " / ")
dat <- scores[scores$condition == "SCI" & scores$program %in% programs &
                scores$method %in% methods, ]
if (nrow(dat) != 48L || anyNA(dat$animal) || anyNA(dat$score) ||
    any(!is.finite(dat$score)) || !identical(sort(unique(dat$animal)), 1:4) ||
    any(!dat$region %in% regions)) stop("SCI design is not four complete animals")
dat$animal <- factor(dat$animal, levels = 1:4)
dat$region <- factor(dat$region, levels = regions)
if (any(table(dat$program, dat$method, dat$animal, dat$region) != 1L))
  stop("Each model requires exactly one score per animal and region")

direction <- function(estimate, label) {
  pieces <- switch(label, "E-A" = c("Epicenter", "Above"),
                   "E-B" = c("Epicenter", "Below"),
                   "A-B" = c("Above", "Below"))
  if (is.null(pieces) || !is.finite(estimate) || estimate == 0) return("tie_or_unavailable")
  paste0(pieces[[1]], if (estimate < 0) "<" else ">", pieces[[2]])
}

rm_anova <- function(d) {
  aa <- aov(score ~ region + Error(animal/region), data = d)
  strata <- summary(aa)[["Error: animal:region"]]
  if (is.null(strata)) stop("Repeated-measures ANOVA stratum absent")
  tt <- strata[[1L]]
  rownames(tt) <- trimws(rownames(tt))
  if (!all(c("region", "Residuals") %in% rownames(tt)))
    stop("Repeated-measures ANOVA table incomplete")
  n <- nlevels(d$animal)
  mse <- tt["Residuals", "Mean Sq"]
  se <- sqrt(2 * mse / n)
  df <- tt["Residuals", "Df"]
  means <- tapply(d$score, d$region, mean)[regions]
  estimates <- vapply(contrasts, function(v) sum(v * means), numeric(1))
  raw_p <- 2 * pt(-abs(estimates / se), df = df)
  ci <- qt(0.975, df = df) * se
  list(F = unname(tt["region", "F value"]),
       num_df = unname(tt["region", "Df"]), den_df = unname(df),
       p = unname(tt["region", "Pr(>F)"]), estimate = estimates,
       SE = rep(se, 3L), lower = estimates - ci, upper = estimates + ci,
       raw_p = raw_p)
}

omnibus <- list(); pairwise <- list(); consistency <- list(); diagnostics <- character()
diagnostics <- c(diagnostics,
  "GSE319931 chronic SCI regional models; 4 biological animals, 3 repeated regions per animal.",
  "Only SCI included. Four prespecified program x method models; no animal or score excluded.",
  "Primary: score ~ region + (1 | animal); lmerTest Satterthwaite omnibus; emmeans contrasts.",
  "Three raw contrast P values receive Holm correction within each program x method.",
  "The 95% confidence limits are unadjusted; Holm correction applies to P values only.",
  "With n=4, residual normality cannot be established by Shapiro testing.",
  "Leave-one-animal-out fits below are influence diagnostics, not replacement analyses.")

for (program in programs) for (method in methods) {
  d <- droplevels(dat[dat$program == program & dat$method == method,
                        c("score", "region", "animal")])
  d$region <- factor(d$region, levels = regions)
  label <- paste(program, method, sep = " / ")
  warnings <- character()
  fit <- withCallingHandlers(
    lmerTest::lmer(score ~ region + (1 | animal), data = d, REML = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  conv <- fit@optinfo$conv$lme4$messages
  if (length(conv)) warnings <- unique(c(warnings, conv))
  singular <- lme4::isSingular(fit, tol = 1e-4)
  fallback <- singular || length(warnings) > 0L
  rm <- rm_anova(d) # independent check and prespecified fallback, never an added model
  if (fallback) {
    model_used <- "repeated_measures_ANOVA_fallback"
    Fv <- rm$F; num_df <- rm$num_df; den_df <- rm$den_df; omnibus_p <- rm$p
    estimates <- rm$estimate; ses <- rm$SE
    lowers <- rm$lower; uppers <- rm$upper; raw_ps <- rm$raw_p
  } else {
    model_used <- "linear_mixed_effects"
    tab <- anova(fit, type = 3, ddf = "Satterthwaite")
    if (!"region" %in% rownames(tab)) stop("Mixed model region test absent")
    Fv <- unname(tab["region", "F value"])
    num_df <- unname(tab["region", "NumDF"])
    den_df <- unname(tab["region", "DenDF"])
    omnibus_p <- unname(tab["region", "Pr(>F)"])
    em <- emmeans(fit, ~ region, lmer.df = "satterthwaite")
    cc <- contrast(em, method = contrasts, adjust = "none")
    ct <- as.data.frame(summary(cc, infer = c(TRUE, TRUE), adjust = "none"))
    if (!identical(as.character(ct$contrast), names(contrasts)))
      stop("Unexpected contrast order")
    estimates <- setNames(ct$estimate, ct$contrast)
    ses <- ct$SE; lowers <- ct$lower.CL; uppers <- ct$upper.CL
    raw_ps <- ct$p.value
    if (abs(Fv - rm$F) > 1e-6 || abs(omnibus_p - rm$p) > 1e-6)
      stop("Mixed model and balanced repeated-measures ANOVA disagree")
  }
  if (any(!is.finite(c(Fv, num_df, den_df, omnibus_p, estimates, ses,
                       lowers, uppers, raw_ps)))) stop("Nonfinite model result")
  holm <- p.adjust(raw_ps, method = "holm")
  omnibus[[label]] <- data.frame(program, method, model_used,
    formula = if (fallback) "score ~ region + Error(animal/region)" else
      "score ~ region + (1 | animal)",
    region_omnibus_F = Fv, df_num = num_df, df_den = den_df,
    region_omnibus_rawP = omnibus_p, singular_fit = singular,
    convergence_warning = length(warnings) > 0L,
    model_warning = paste(warnings, collapse = " | "),
    n_obs = nrow(d), n_animals = nlevels(d$animal))
  one <- data.frame(program, method, contrast = names(contrasts),
    estimate = unname(estimates), SE = unname(ses),
    CI_lower = unname(lowers), CI_upper = unname(uppers),
    raw_P = unname(raw_ps), Holm_P = unname(holm))
  pairwise[[label]] <- one

  psub <- paired[paired$program == program & paired$method == method, ]
  if (nrow(psub) != 3L || !setequal(psub$comparison, names(contrasts)))
    stop("Paired-result coverage mismatch")
  psub <- psub[match(names(contrasts), psub$comparison), ]
  if (any(abs(psub$mean_diff - one$estimate) > 1e-8))
    stop("Mixed and paired mean differences disagree")
  mixed_dir <- mapply(direction, one$estimate, one$contrast)
  paired_dir <- mapply(direction, psub$mean_diff, psub$comparison)
  consistency[[label]] <- data.frame(program, method,
    contrast = names(contrasts), mixed_model_direction = mixed_dir,
    paired_direction = paired_dir,
    direction_consistent = mixed_dir == paired_dir,
    mixed_Holm_P = one$Holm_P, paired_P = psub$paired_t_p)

  res <- residuals(fit)
  scaled_res <- res / sigma(fit)
  loo <- do.call(rbind, lapply(levels(d$animal), function(animal_id) {
    reduced <- droplevels(d[d$animal != animal_id, ])
    reduced$region <- factor(reduced$region, levels = regions)
    means <- tapply(reduced$score, reduced$region, mean)[regions]
    effects <- vapply(contrasts, function(v) sum(v * means), numeric(1))
    data.frame(animal_removed = animal_id, contrast = names(contrasts),
      estimate = unname(effects), direction = mapply(direction, effects,
                                                   names(contrasts)),
      direction_flips = mapply(direction, effects, names(contrasts)) != mixed_dir,
      omnibus_P = rm_anova(reduced)$p)
  }))
  diagnostics <- c(diagnostics,
    paste0("\n--- ", label, " ---"),
    paste("model_used:", model_used),
    paste("singular_fit:", singular),
    paste("warnings:", if (length(warnings)) paste(warnings, collapse = " | ") else "none"),
    sprintf("residuals: range %.6g to %.6g; SD %.6g; max |residual/sigma| %.3f",
            min(res), max(res), sd(res), max(abs(scaled_res))),
    paste("max_scaled_residual_animal:", as.character(d$animal[which.max(abs(scaled_res))])),
    sprintf("RM-ANOVA cross-check: F(%d,%d)=%.6g, P=%.8g",
            rm$num_df, rm$den_df, rm$F, rm$p),
    paste("leave-one-animal-out contrast directions:",
          paste(capture.output(print(loo, row.names = FALSE)), collapse = "\n")),
    paste("any direction flip:", any(loo$direction_flips)),
    paste("direction-flip animal/contrast:",
          if (any(loo$direction_flips))
            paste(paste0(loo$animal_removed[loo$direction_flips], "/",
                         loo$contrast[loo$direction_flips]), collapse = ", ")
          else "none"),
    "Omnibus P can vary materially after removing one of only four animals; no animal was removed from primary results.")
}

write.csv(do.call(rbind, omnibus),
          file.path(output_dir, "GSE319931_mixed_model_omnibus.csv"), row.names = FALSE)
write.csv(do.call(rbind, pairwise),
          file.path(output_dir, "GSE319931_mixed_model_pairwise.csv"), row.names = FALSE)
write.csv(do.call(rbind, consistency),
          file.path(output_dir, "GSE319931_mixed_vs_paired_consistency.csv"), row.names = FALSE)
writeLines(diagnostics, file.path(output_dir, "GSE319931_mixed_model_diagnostics.txt"))
writeLines(c("GSE319931 mixed-model session info", capture.output(sessionInfo())),
           file.path(output_dir, "GSE319931_mixed_model_sessionInfo.txt"))
