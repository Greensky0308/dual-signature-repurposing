#!/usr/bin/env Rscript
# A2 step 2: limma differential expression on the same two cohorts, for direct
# comparison against the manuscript's Welch's t-test pipeline.
suppressPackageStartupMessages(library(limma))

OUT <- Sys.getenv("HCC_HF_OUT", "output")   # run from the pipeline root

for (tag in c("GSE14520", "GSE57345")) {
  expr <- read.csv(file.path(OUT, paste0("A2_", tag, "_expr.csv")),
                   row.names = 1, check.names = FALSE)
  grp <- read.csv(file.path(OUT, paste0("A2_", tag, "_groups.csv")))
  stopifnot(identical(colnames(expr), grp$sample))

  # disease coded as the second level so coef2 = disease - control
  f <- factor(grp$group, levels = c("control", "disease"))
  design <- model.matrix(~ f)
  fit <- lmFit(as.matrix(expr), design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  tt$probe <- rownames(tt)
  tt$padj <- p.adjust(tt$P.Value, method = "BH")
  out <- tt[, c("probe", "logFC", "P.Value", "padj", "AveExpr", "t")]
  names(out) <- c("probe", "log2FC", "p", "padj", "mean_avg", "t_limma")
  write.csv(out, file.path(OUT, paste0("A2_limma_", tag, ".csv")), row.names = FALSE)
  sig <- sum(abs(out$log2FC) > 1 & out$padj < 0.05, na.rm = TRUE)
  cat(sprintf("%s: probes %d  |log2FC|>1 & padj<0.05: %d\n", tag, nrow(out), sig))
}
