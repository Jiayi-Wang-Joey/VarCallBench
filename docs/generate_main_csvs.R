#!/usr/bin/env Rscript
# =============================================================================
# generate_main_csvs.R
# Generates corrected main-figure CSVs using DP>=5 filtered VCFs.
#
# Outputs (in OUTDIR):
#   cross_platform_detail_summary.csv  — bar plot (same format as existing)
#   cross_caller_incidence.csv         — upset incidence matrix
#
# Usage:
#   apptainer exec ~/VarCallbench/envs/r_plot.sif \
#       Rscript ~/VarCallbench/docs/generate_main_csvs.R
# =============================================================================
suppressPackageStartupMessages(library(data.table))

DP        <- 5
BENCH     <- path.expand("~/VarCallbench")
DP_DIR    <- file.path(BENCH, "out/dp_sweep/dp_filtered")
OUTDIR    <- file.path(BENCH, "out/main_csvs_dp5")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

CELL_LINES <- c("H211", "H526")
CALLERS    <- c("clair3_rna", "deep_variant", "longcallR", "longcallR_nn", "isolaser")
PLATFORMS  <- c("cDNA", "MasSeq", "dRNA")

# Clean caller name map
CALLER_LABELS <- c(
    clair3_rna   = "Clair3-RNA",
    deep_variant = "DeepVariant",
    longcallR    = "longcallR",
    longcallR_nn = "longcallR-nn",
    isolaser     = "isoLASER"
)

# Read variant IDs from a DP-filtered VCF
read_ids <- function(vcf_path) {
    if (!file.exists(vcf_path)) {
        warning("Missing: ", vcf_path)
        return(character(0))
    }
    cmd <- sprintf("bcftools query -f '%%CHROM\t%%POS\t%%REF\t%%ALT\n' %s", shQuote(vcf_path))
    dt  <- fread(cmd = cmd, sep = "\t", header = FALSE,
                 col.names = c("CHROM","POS","REF","ALT"), showProgress = FALSE)
    if (nrow(dt) == 0L) return(character(0))
    unique(paste(dt$CHROM, dt$POS, dt$REF, dt$ALT, sep = "_"))
}

vcf_path <- function(cl, caller, plat, dp = DP) {
    file.path(DP_DIR, sprintf("%s_%s_%s_dp%d.vcf.gz", cl, caller, plat, dp))
}

# =============================================================================
# 1. Cross-platform detail summary (bar plot)
# =============================================================================
message("=== Cross-platform detail summary ===")

xplat_rows <- list()

for (cl in CELL_LINES) {
    for (caller in CALLERS) {
        ids <- lapply(PLATFORMS, function(p) {
            # isolaser H526 dRNA excluded (tool failure)
            if (caller == "isolaser" && cl == "H526" && p == "dRNA") return(character(0))
            read_ids(vcf_path(cl, caller, p))
        })
        names(ids) <- PLATFORMS

        v_cdna   <- ids[["cDNA"]]
        v_masseq <- ids[["MasSeq"]]
        v_drna   <- ids[["dRNA"]]

        shared_all   <- Reduce(intersect, ids)
        shared_cm    <- setdiff(intersect(v_cdna,   v_masseq), v_drna)
        shared_cd    <- setdiff(intersect(v_cdna,   v_drna),   v_masseq)
        shared_md    <- setdiff(intersect(v_masseq, v_drna),   v_cdna)
        only_cdna    <- setdiff(v_cdna,   union(v_masseq, v_drna))
        only_masseq  <- setdiff(v_masseq, union(v_cdna,   v_drna))
        only_drna    <- setdiff(v_drna,   union(v_cdna,   v_masseq))

        counts <- c(
            cDNA_MasSeq_dRNA = length(shared_all),
            cDNA_MasSeq      = length(shared_cm),
            cDNA_dRNA        = length(shared_cd),
            MasSeq_dRNA      = length(shared_md),
            cDNA_only        = length(only_cdna),
            MasSeq_only      = length(only_masseq),
            dRNA_only        = length(only_drna)
        )
        total <- sum(counts)

        xplat_rows[[length(xplat_rows) + 1]] <- data.table(
            cell_line    = cl,
            caller       = CALLER_LABELS[caller],
            detail_class = names(counts),
            N            = as.integer(counts),
            proportion   = if (total > 0) counts / total else rep(NA_real_, length(counts))
        )
        message(sprintf("  %s | %s | shared_all=%d / total=%d (%.1f%%)",
                        cl, caller, counts["cDNA_MasSeq_dRNA"], total,
                        100 * counts["cDNA_MasSeq_dRNA"] / max(total, 1)))
    }
}

xplat_dt <- rbindlist(xplat_rows)
xplat_out <- file.path(OUTDIR, "cross_platform_detail_summary.csv")
fwrite(xplat_dt, xplat_out)
message("Wrote: ", xplat_out)

# =============================================================================
# 2. Cross-caller incidence matrix (upset plot)
# =============================================================================
message("\n=== Cross-caller incidence matrix ===")

incidence_rows <- list()

for (cl in CELL_LINES) {
    for (plat in PLATFORMS) {
        # Read all callers for this cell_line × platform
        id_list <- lapply(CALLERS, function(caller) {
            # isolaser H526 dRNA excluded (tool failure)
            if (caller == "isolaser" && cl == "H526" && plat == "dRNA") return(character(0))
            read_ids(vcf_path(cl, caller, plat))
        })
        names(id_list) <- CALLERS

        all_variants <- unique(unlist(id_list))
        if (length(all_variants) == 0) next

        dt <- data.table(
            cell_line  = cl,
            platform   = plat,
            variant_id = all_variants
        )
        for (caller in CALLERS) {
            dt[[CALLER_LABELS[caller]]] <- as.integer(all_variants %in% id_list[[caller]])
        }

        incidence_rows[[length(incidence_rows) + 1]] <- dt
        message(sprintf("  %s | %s | %d unique variants across %d callers",
                        cl, plat, length(all_variants), length(CALLERS)))
    }
}

incidence_dt  <- rbindlist(incidence_rows, fill = TRUE)
incidence_out <- file.path(OUTDIR, "cross_caller_incidence.csv")
fwrite(incidence_dt, incidence_out)
message("Wrote: ", incidence_out)

message("\nDone. All files in: ", OUTDIR)
