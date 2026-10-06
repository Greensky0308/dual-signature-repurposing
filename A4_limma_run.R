#!/usr/bin/env Rscript
# A4 step 2: limma differential expression in the independent validation cohorts.
#
# limma is used for every cohort (main analysis and validation) so that the
# pipeline has a single differential-expression method. A2 showed that limma and
# the original Welch's t-test give identical signatures on the discovery
# cohorts, so this switch does not change any reported result.
suppressPackageStartupMessages(library(limma))

OUT <- Sys.getenv("HCC_HF_OUT", "output")   # run from the pipeline root

run <- function(tag, disease_level, control_level, quantile_normalise = FALSE,
                min_mean = NULL) {
  expr <- read.csv(file.path(OUT, paste0("A4_", tag, "_expr.csv")),
                   row.names = 1, check.names = FALSE)
  grp <- read.csv(file.path(OUT, paste0("A4_", tag, "_groups.csv")))
  stopifnot(identical(colnames(expr), grp$sample))
  expr <- as.matrix(expr)
  n_na <- sum(!complete.cases(expr))
  if (n_na > 0) {
    cat(sprintf("%s: dropping %d genes with missing values\n", tag, n_na))
    expr <- expr[complete.cases(expr), , drop = FALSE]
  }
  if (!is.null(min_mean)) {
    # drop genes that are not expressed: standard low-expression filter for RNA-seq
    keep <- rowMeans(expr, na.rm = TRUE) >= min_mean
    cat(sprintf("%s: low-expression filter (>%.1f) removed %d of %d genes\n",
                tag, min_mean, sum(!keep), length(keep)))
    expr <- expr[keep, , drop = FALSE]
  }
  if (quantile_normalise) {
    expr <- normalizeBetweenArrays(expr, method = "quantile")
  }
  f <- factor(grp$group, levels = c(control_level, disease_level))
  design <- model.matrix(~ f)
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  out <- data.frame(probe = rownames(tt), log2FC = tt$logFC, p = tt$P.Value,
                    padj = p.adjust(tt$P.Value, method = "BH"),
                    t_limma = tt$t, stringsAsFactors = FALSE)
  write.csv(out, file.path(OUT, paste0("A4_limma_", tag, ".csv")), row.names = FALSE)
  cat(sprintf("%s: probes %d  |log2FC|>1 & padj<0.05: %d\n", tag, nrow(out),
              sum(abs(out$log2FC) > 1 & out$padj < 0.05, na.rm = TRUE)))
}

run("GSE76427", "tumour", "normal", quantile_normalise = TRUE)
run("GSE141910", "disease", "normal", min_mean = 1)
