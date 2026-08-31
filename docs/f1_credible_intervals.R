suppressPackageStartupMessages({
    library(data.table)
})

dt <- fread("/home/jiayiwang/VarCallbench/out/happy_summary_collector/9938922d/happy_summary_merged.csv")
dt <- dt[Filter == "PASS"]

label_fn <- function(x) fcase(
    x == "clair3_rna",   "Clair3-RNA",
    x == "deep_variant", "DeepVariant",
    x == "longcallR",    "longcallR",
    x == "longcallR_nn", "longcallR-nn",
    x == "isolaser",     "isoLASER",
    default = x
)
dt[, tool := label_fn(tool)]
dt[, cell_line := sub("^([A-Za-z0-9]+)_.*$", "\\1", sample)]
dt[, platform := fcase(
    grepl("bulk_PB$", sample), "MasSeq",
    grepl("bulk_ONT$", sample), "cDNAxR10",
    grepl("dRNA_ONT$", sample), "dRNA004"
)]
dt[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]

# F1 = 2*TP / (TRUTH.TOTAL + QUERY.TOTAL); treat as 2x a binomial proportion of
# "successes" (TP) out of N = TRUTH.TOTAL + QUERY.TOTAL "trials", and use a
# Jeffreys-style Beta(TP+0.5, N-TP+0.5) credible interval on that proportion,
# doubled to get the F1 credible interval. Same approach as the reference code.
dt[, N := TRUTH.TOTAL + QUERY.TOTAL]

# longcallR/longcallR-nn produce no usable INDEL calls -- drop that block entirely.
dt <- dt[!(tool %in% c("longcallR", "longcallR-nn") & Type == "INDEL")]
# isoLASER's INDEL calls are only usable on the PacBio platform (MasSeq) --
# same restriction used throughout this analysis.
dt <- dt[!(tool == "isoLASER" & Type == "INDEL" & platform != "MasSeq")]

# Pooled across all matching LongBench samples per (tool, Type, coverage), since
# per-sample N is often small -- pooling gives a more stable, publication-
# appropriate interval than averaging per-sample CIs. n=6 samples (2 cell lines
# x 3 platforms) for every row except isoLASER INDEL, which pools n=2 (2 cell
# lines x MasSeq only).
pooled <- dt[, .(TP = sum(TRUTH.TP), N = sum(N)), by = .(tool, Type, coverage)]
pooled[, F1 := 2 * TP / N]
pooled[, F1_lower := 2 * qbeta(0.025, TP + 0.5, (N - TP) + 0.5)]
pooled[, F1_upper := 2 * qbeta(0.975, TP + 0.5, (N - TP) + 0.5)]

tool_levels <- c("Clair3-RNA", "DeepVariant", "isoLASER", "longcallR", "longcallR-nn")
pooled[, tool := factor(tool, levels = tool_levels)]
setorder(pooled, Type, tool, coverage)

pooled[, F1_fmt := sprintf("%.3f [%.3f, %.3f]", F1, F1_lower, F1_upper)]

cat("=== Pooled F1 with 95% credible interval (Jeffreys/Beta), by Type x Caller x Coverage ===\n")
print(pooled[, .(Type, tool, coverage, F1_fmt)], nrows = Inf)

fwrite(pooled[, .(Type, tool, coverage, TP, N, F1, F1_lower, F1_upper)],
       "/home/jiayiwang/VarCallbench/out/f1_credible_intervals.csv")

## ---- LaTeX table: rows = Type/Caller, columns = coverage tiers ----
wide <- dcast(pooled, Type + tool ~ coverage, value.var = "F1_fmt")
covs <- c("1", "5", "10", "30", "50", "100")

latex_lines <- c(
    "\\begin{table}[H]",
    "\\centering",
    "\\small",
    paste0("\\begin{tabular}{l l ", paste(rep("r", length(covs)), collapse = " "), "}"),
    "\\hline",
    paste("Type & Variant Caller &", paste(covs, collapse = " & "), "\\\\"),
    "\\hline"
)

for (ty in levels(pooled$Type)) {
    sub <- wide[Type == ty]
    sub <- sub[order(tool)]
    n_rows <- nrow(sub)
    latex_lines <- c(latex_lines, sprintf("\\multirow{%d}{*}{%s}", n_rows, ty))
    for (i in seq_len(n_rows)) {
        vals <- sapply(covs, function(cv) {
            v <- sub[[cv]][i]
            if (is.na(v)) "--" else v
        })
        latex_lines <- c(latex_lines, paste(" &", sub$tool[i], "&", paste(vals, collapse = " & "), "\\\\"))
    }
    latex_lines <- c(latex_lines, "\\hline")
}
latex_lines <- c(latex_lines,
    "\\end{tabular}",
    "\\bcaption{F1 score with 95\\% credible intervals across LongBench samples.}{F1 = 2*TP / (TRUTH.TOTAL + QUERY.TOTAL), pooled across samples at each coverage cutoff (n=6: H211, H526 x MasSeq, cDNAxR10, dRNA004), except isoLASER INDEL, pooled across n=2 (H211, H526 x MasSeq only), since isoLASER's INDEL calls are not usable on ONT platforms. longcallR and longcallR-nn produce no usable INDEL calls and are omitted from that block. 95\\% credible intervals computed via a Jeffreys-style Beta(TP + 0.5, N - TP + 0.5) posterior on the underlying proportion, doubled to the F1 scale.}",
    "\\label{tab:f1_credible_intervals}",
    "\\end{table}"
)
writeLines(latex_lines, "/home/jiayiwang/VarCallbench/out/f1_credible_intervals.tex")
cat("\nLaTeX written to out/f1_credible_intervals.tex\n")
