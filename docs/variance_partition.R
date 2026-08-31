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
dt[, chemistry := fcase(
    grepl("bulk_PB$", sample), "MasSeq",
    grepl("bulk_ONT$", sample), "cDNAxR10",
    grepl("dRNA_ONT$", sample), "dRNA004"
)]
dt[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]
dt[, coverage := factor(coverage, levels = c(1, 5, 10, 30, 50, 100))]
dt[, tool := factor(tool, levels = c("DeepVariant", "Clair3-RNA", "isoLASER", "longcallR", "longcallR-nn"))]
dt[, chemistry := factor(chemistry, levels = c("cDNAxR10", "dRNA004", "MasSeq"))]
dt[, cell_line := factor(cell_line)]
dt[, logTP := log10(TRUTH.TP + 1)]

m <- dt[, .(METRIC.F1_Score, logTP, tool, chemistry, Type, coverage, cell_line)]
m <- na.omit(m)
cat("n =", nrow(m), "observations across", uniqueN(dt$sample), "samples\n\n")

fit_f1  <- lm(METRIC.F1_Score ~ tool + chemistry + Type + coverage + cell_line, data = m)
fit_tp  <- lm(logTP ~ tool + chemistry + Type + coverage + cell_line, data = m)

cat("Full model R^2 -- F1:", round(summary(fit_f1)$r.squared, 3),
    " logTP:", round(summary(fit_tp)$r.squared, 3), "\n\n")

terms <- c("tool", "chemistry", "Type", "coverage", "cell_line")
term_labels <- c(
    tool = "Variant caller", chemistry = "Chemistry", Type = "Variant type (SNV/INDEL)",
    coverage = "Coverage cutoff", cell_line = "Cell line"
)

partial_r2 <- function(fit, data, terms) {
    full_sse <- sum(residuals(fit)^2)
    sapply(terms, function(tm) {
        f_reduced <- update(fit, as.formula(paste(". ~ . -", tm)))
        reduced_sse <- sum(residuals(f_reduced)^2)
        (reduced_sse - full_sse) / reduced_sse
    })
}

pr2_f1 <- partial_r2(fit_f1, m, terms)
pr2_tp <- partial_r2(fit_tp, m, terms)

tab1 <- data.table(
    Term = term_labels[terms],
    F1_partial_R2 = round(pr2_f1[terms], 3),
    logTP_partial_R2 = round(pr2_tp[terms], 3)
)
tab1 <- tab1[order(-F1_partial_R2)]
cat("=== Partial R^2 (Table 1) ===\n")
print(tab1)
cat("\n")

coef_table <- function(fit) {
    s <- summary(fit)$coefficients
    dt_c <- data.table(term = rownames(s), estimate = s[, "Estimate"], p = s[, "Pr(>|t|)"])
    dt_c[term != "(Intercept)"]
}

cf1 <- coef_table(fit_f1)
setnames(cf1, c("estimate", "p"), c("F1_estimate", "F1_p"))
cftp <- coef_table(fit_tp)
setnames(cftp, c("estimate", "p"), c("logTP_estimate", "logTP_p"))
tab2 <- merge(cf1, cftp, by = "term", sort = FALSE)

cat("=== Coefficients (Table 2) ===\n")
cat("Reference levels: caller =", levels(m$tool)[1], ", chemistry =", levels(m$chemistry)[1],
    ", type = SNV, coverage =", levels(m$coverage)[1], ", cell line =", levels(m$cell_line)[1], "\n")
print(tab2, digits = 3)

fwrite(tab1, "/home/jiayiwang/VarCallbench/out/variance_partition_table1.csv")
fwrite(tab2, "/home/jiayiwang/VarCallbench/out/variance_partition_table2.csv")
