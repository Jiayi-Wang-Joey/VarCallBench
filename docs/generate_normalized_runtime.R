suppressPackageStartupMessages({
    library(data.table)
})

## Only H211 is used for the normalized-runtime comparison.
perf <- fread("out/performances.tsv")
qc   <- fread("out/alignment_qc_collector/alignment_qc_merged.csv")

h211_datasets <- c("H211_bulk_ONT", "H211_bulk_PB", "H211_dRNA_ONT")

perf <- perf[dataset %in% h211_datasets]

## dataset -> cell_line / platform
perf[, platform := fcase(
    dataset == "H211_bulk_ONT", "cDNA",
    dataset == "H211_bulk_PB",  "MasSeq",
    dataset == "H211_dRNA_ONT", "dRNA"
)]
perf[, cell_line := "H211"]

## caller labels; isoLASER = isolaser_call ("isolaser run") only, not isolaser_annotate
perf[, caller := fcase(
    module == "clair3_rna",    "Clair3-RNA",
    module == "deep_variant",  "DeepVariant",
    module == "longcallR",     "longcallR",
    module == "longcallR_nn",  "longcallR-nn",
    module == "isolaser",      "isoLASER"
)]
perf <- perf[!is.na(caller)]

## max_pss (proportional set size), not max_rss: RSS double-counts shared
## memory across forked worker processes (e.g. isoLASER's workers sharing a
## loaded transcriptome/genome index), inflating multiprocess tools' memory
## relative to single-process tools like longcallR. plot_time.R already uses
## max_pss for its memory panel, so this keeps both consistent.
agg <- perf[, .(
    wall_time_s = sum(s),
    cpu_time_s  = sum(cpu_time),
    max_pss_MB  = max(max_pss)
), by = .(caller, cell_line, platform, dataset)]

agg <- merge(agg, qc[, .(dataset = dataset_id, mapped_reads, average_length)], by = "dataset")
agg[, dataset := NULL]
## mapped_reads/1e6 first to avoid integer overflow, and to keep this column
## human-sized (millions of bases) rather than a raw ~1e10 base count
agg[, mapped_Mbases := round((mapped_reads / 1e6) * average_length, 2)]

agg[, `:=`(
    runtime_per_M_reads_s  = round(wall_time_s / (mapped_reads / 1e6), 2),
    memory_per_M_reads_MB  = round(max_pss_MB  / (mapped_reads / 1e6), 2),
    cputime_per_M_reads_s  = round(cpu_time_s  / (mapped_reads / 1e6), 2),
    runtime_per_Gb_s  = round(wall_time_s / (mapped_Mbases / 1e3), 2),
    memory_per_Gb_MB  = round(max_pss_MB  / (mapped_Mbases / 1e3), 2),
    cputime_per_Gb_s  = round(cpu_time_s  / (mapped_Mbases / 1e3), 2)
)]
agg[, n_runs := NA_real_]

## per-caller average across the 3 H211 platforms
avg <- agg[, .(
    cell_line = "Average",
    platform  = "All platforms",
    mapped_reads = NA_integer_,
    mapped_Mbases = NA_real_,
    wall_time_s  = round(mean(wall_time_s), 2),
    max_pss_MB   = round(mean(max_pss_MB), 2),
    runtime_per_M_reads_s = round(mean(runtime_per_M_reads_s), 2),
    memory_per_M_reads_MB = round(mean(memory_per_M_reads_MB), 2),
    cputime_per_M_reads_s = round(mean(cputime_per_M_reads_s), 2),
    runtime_per_Gb_s = round(mean(runtime_per_Gb_s), 2),
    memory_per_Gb_MB = round(mean(memory_per_Gb_MB), 2),
    cputime_per_Gb_s = round(mean(cputime_per_Gb_s), 2),
    n_runs = .N
), by = caller]

agg[, `:=`(cpu_time_s = NULL, average_length = NULL)]
out <- rbind(avg, agg, fill = TRUE)
out <- out[, .(caller, cell_line, platform, mapped_reads, mapped_Mbases, wall_time_s, max_pss_MB,
               runtime_per_M_reads_s, memory_per_M_reads_MB, cputime_per_M_reads_s,
               runtime_per_Gb_s, memory_per_Gb_MB, cputime_per_Gb_s, n_runs)]

out[, is_avg := as.integer(cell_line == "Average")]
setorder(out, caller, -is_avg)
out[, is_avg := NULL]

fwrite(out, "out/normalized_runtime_v2.csv")
print(out)
