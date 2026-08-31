suppressPackageStartupMessages({
    library(data.table)
})

dt <- fread("/home/jiayiwang/VarCallbench/out/happy_summary_collector/9938922d/happy_summary_merged.csv")
dt <- dt[Filter == "PASS" & Type == "SNP"]

label_fn <- function(x) fcase(
    x == "clair3_rna",   "Clair3-RNA",
    x == "deep_variant", "DeepVariant",
    x == "longcallR",    "longcallR",
    x == "longcallR_nn", "longcallR-nn",
    x == "isolaser",     "isoLASER",
    default = x
)
dt[, tool := label_fn(tool)]
dt[, dataset := fcase(
    sample == "H526_dRNA_ONT", "H526-dRNA004",
    sample == "H526_bulk_PB",  "H526-MasSeq",
    sample == "H526_bulk_ONT", "H526-cDNAxR10",
    sample == "H211_dRNA_ONT", "H211-dRNA004",
    sample == "H211_bulk_PB",  "H211-MasSeq",
    sample == "H211_bulk_ONT", "H211-cDNAxR10",
    default = sample
)]

# Per-sample (not pooled) F1 credible interval -- same Jeffreys-Beta method as before.
dt[, N := TRUTH.TOTAL + QUERY.TOTAL]
dt[, F1_lower := 2 * qbeta(0.025, TRUTH.TP + 0.5, (N - TRUTH.TP) + 0.5)]
dt[, F1_upper := 2 * qbeta(0.975, TRUTH.TP + 0.5, (N - TRUTH.TP) + 0.5)]
dt[, F1_fmt := sprintf("%.3f [%.3f, %.3f]", METRIC.F1_Score, F1_lower, F1_upper)]

tool_levels <- c("Clair3-RNA", "DeepVariant", "isoLASER", "longcallR", "longcallR-nn")
dt[, tool := factor(tool, levels = tool_levels)]

dataset_order <- c("H211-MasSeq", "H211-cDNAxR10", "H211-dRNA004",
                    "H526-MasSeq", "H526-cDNAxR10", "H526-dRNA004")
dt[, dataset := factor(dataset, levels = dataset_order)]

wide <- dcast(dt, dataset + coverage ~ tool, value.var = "F1_fmt")
setorder(wide, dataset, coverage)

fwrite(wide, "/home/jiayiwang/VarCallbench/out/f1_ci_by_dataset.csv")

## ---- LaTeX: rows grouped by dataset (multirow), coverage as sub-rows, callers as columns ----
## Wrapped in landscape (requires \usepackage{pdflscape} in the preamble) since
## 5 caller columns of "F1 [lo, hi]" strings don't fit a portrait page width.
latex_lines <- c(
    "\\begin{landscape}",
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    paste0("\\begin{tabular}{l l ", paste(rep("l", length(tool_levels)), collapse = " "), "}"),
    "\\hline",
    paste("Dataset & Coverage &", paste(tool_levels, collapse = " & "), "\\\\"),
    "\\hline"
)

for (ds in dataset_order) {
    sub <- wide[dataset == ds]
    n_rows <- nrow(sub)
    if (n_rows == 0) next
    latex_lines <- c(latex_lines, sprintf("\\multirow{%d}{*}{%s}", n_rows, ds))
    for (i in seq_len(n_rows)) {
        vals <- sapply(tool_levels, function(tl) {
            v <- sub[[tl]][i]
            if (is.na(v)) "--" else v
        })
        latex_lines <- c(latex_lines, paste(" &", sub$coverage[i], "&", paste(vals, collapse = " & "), "\\\\"))
    }
    latex_lines <- c(latex_lines, "\\hline")
}
latex_lines <- c(latex_lines,
    "\\end{tabular}",
    "\\bcaption{SNV F1 score with 95\\% credible intervals, by sample and coverage cutoff.}{F1 = 2*TP / (TRUTH.TOTAL + QUERY.TOTAL); 95\\% credible intervals computed per sample and coverage cutoff via a Jeffreys-style Beta(TP + 0.5, N - TP + 0.5) posterior on the underlying proportion, doubled to the F1 scale. Unlike the pooled summary table, these intervals are per-sample and therefore wider.}",
    "\\label{tab:f1_ci_by_dataset}",
    "\\end{table}",
    "\\end{landscape}"
)
writeLines(latex_lines, "/home/jiayiwang/VarCallbench/out/f1_ci_by_dataset.tex")
cat("Done -- out/f1_ci_by_dataset.tex\n")
