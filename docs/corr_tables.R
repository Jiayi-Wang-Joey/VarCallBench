suppressPackageStartupMessages({
    library(data.table)
})

perf <- fread("/home/jiayiwang/VarCallbench/out/happy_summary_collector/9938922d/happy_summary_merged.csv")
qual <- fread("/home/jiayiwang/VarCallbench/out/alignment_qc_collector/alignment_qc_merged.csv")
base <- fread("/home/jiayiwang/VarCallbench/out/callable_bases.csv")
perf <- perf[Filter == "PASS"]

label_fn <- function(x) fcase(
    x == "clair3_rna",   "Clair3-RNA",
    x == "deep_variant", "DeepVariant",
    x == "longcallR",    "longcallR",
    x == "longcallR_nn", "longcallR-nn",
    x == "isolaser",     "isoLASER",
    default = x
)
perf[, tool := label_fn(tool)]
perf[, platform := fcase(
    grepl("bulk_PB$", sample), "MasSeq",
    grepl("bulk_ONT$", sample), "cDNAxR10",
    grepl("dRNA_ONT$", sample), "dRNA004"
)]
perf[, Type := factor(Type, levels = c("SNP", "INDEL"), labels = c("SNV", "INDEL"))]

callers_snv   <- c("Clair3-RNA", "DeepVariant", "isoLASER", "longcallR", "longcallR-nn")
callers_indel <- c("Clair3-RNA", "DeepVariant", "isoLASER")   # longcallR/-nn have no usable INDEL calls, same exclusion as plot_LB-cor.R
covs <- c(1, 5, 10, 30, 50, 100)

spearman_or_na <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(c(rho = NA_real_, n = sum(ok)))
    ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
    c(rho = unname(ct$estimate), n = sum(ok))
}

## ---- Table 1: corr(1 - error_rate, F1) ----
dt1 <- perf[qual, on = .(sample = dataset_id)]
dt1[, one_minus_err := 1 - error_rate]

rows1 <- list()
for (ty in c("SNV", "INDEL")) {
    callers <- if (ty == "SNV") callers_snv else callers_indel
    for (cl in callers) {
        sub0 <- dt1[Type == ty & tool == cl]
        # isoLASER INDEL: restrict to PacBio (MasSeq), same restriction used elsewhere for this dataset
        if (ty == "INDEL" && cl == "isoLASER") sub0 <- sub0[platform == "MasSeq"]
        vals <- sapply(covs, function(cv) {
            sub <- sub0[coverage == cv]
            unname(spearman_or_na(sub$one_minus_err, sub$METRIC.F1_Score)["rho"])
        })
        row <- as.data.table(as.list(vals))
        setnames(row, as.character(covs))
        rows1[[paste(ty, cl)]] <- cbind(data.table(Type = ty, Caller = cl), row)
    }
}
tab1 <- rbindlist(rows1)

## ---- Table 2: corr(callable bases, TP) ----
dt2 <- perf[base, on = .(sample = dataset, coverage = coverage)]
dt2[, platform := fcase(
    grepl("bulk_PB$", sample), "MasSeq",
    grepl("bulk_ONT$", sample), "cDNAxR10",
    grepl("dRNA_ONT$", sample), "dRNA004"
)]

rows2 <- list()
for (ty in c("SNV", "INDEL")) {
    callers <- if (ty == "SNV") callers_snv else callers_indel
    for (cl in callers) {
        sub0 <- dt2[Type == ty & tool == cl]
        if (ty == "INDEL" && cl == "isoLASER") sub0 <- sub0[platform == "MasSeq"]
        vals <- sapply(covs, function(cv) {
            sub <- sub0[coverage == cv]
            unname(spearman_or_na(sub$bases, sub$TRUTH.TP)["rho"])
        })
        row <- as.data.table(as.list(vals))
        setnames(row, as.character(covs))
        rows2[[paste(ty, cl)]] <- cbind(data.table(Type = ty, Caller = cl), row)
    }
}
tab2 <- rbindlist(rows2)

cat("=== Table 1: Spearman corr(1 - error_rate, F1) ===\n")
cat("n per row = 6 (2 cell lines x 3 platforms), isoLASER INDEL n = 2 (MasSeq only, cell line only) -- unstable, small-n\n")
print(tab1, digits = 2)
cat("\n=== Table 2: Spearman corr(callable bases, TP) ===\n")
cat("n per row = 6 (2 cell lines x 3 platforms), isoLASER INDEL n = 2 (MasSeq only)\n")
print(tab2, digits = 2)

fwrite(tab1, "/home/jiayiwang/VarCallbench/out/corr_errorrate_F1.csv")
fwrite(tab2, "/home/jiayiwang/VarCallbench/out/corr_bases_TP.csv")
