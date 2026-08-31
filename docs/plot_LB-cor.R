suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
})
setwd("~/Desktop/roland/VarCallbench/")
qual <- fread("out/alignment_qc_collector/alignment_qc_merged.csv")
perf <- fread("out/happy_summary_collector/9938922d/happy_summary_merged.csv")
base <- fread("out/callable_bases.csv")

dt1 <- perf[base, on = .(sample = dataset, coverage = coverage)]
pal <- c(
    "dRNA004" = "#fec44f",
    "cDNAxR10" = "#d95f0e",
    "MasSeq" = "#54278f"
)

dt1[, sample := fcase(
    sample == "H526_dRNA_ONT", "H526-dRNA004",
    sample == "H526_bulk_PB",  "H526-MasSeq",
    sample == "H526_bulk_ONT", "H526-cDNAxR10",
    sample == "H211_dRNA_ONT", "H211-dRNA004",
    sample == "H211_bulk_PB",  "H211-MasSeq",
    sample == "H211_bulk_ONT", "H211-cDNAxR10",
    default = sample
)]
dt1 <- dt1[Filter == "PASS"]
dt1[, tech := sub("^[^-]+-", "", sample)]
label_fn <- function(x) fcase(
    x == "clair3_rna",   "Clair3-RNA",
    x == "deep_variant", "DeepVariant",
    x == "longcallR",    "longcallR",     
    x == "longcallR_nn", "longcallR-nn",
    x == "isolaser",     "isoLASER",
    default = x
)

dt1[, tool := label_fn(tool)]
td <- dt1[!(grepl("longcallR", tool) & Type =="INDEL")]
td <- td[!(grepl("dRNA004|cDNA", tech) & 
               tool == "isoLASER" & Type == "INDEL")]
td[, Type := factor(Type, levels = c("SNP", "INDEL"), 
                    labels = c("SNV", "INDEL"))]
gg <- ggplot(td, aes(bases/1e6, TRUTH.TP, shape = tool, col = tech)) +
    geom_point(size = 2, alpha = 0.7) +
    facet_grid2(Type ~ coverage, scales = "free") +
    theme_classic() +
    scale_color_manual(values = pal) +
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.8
        ),
        strip.background = element_rect(
            fill = "white",
            color = "black",
            linewidth = 0.8),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        aspect.ratio = 1
    ) +
    labs(col = "Chemistry", x = "Callable Bases (Million)", y = "Number of TPs",
         shape = "Variant Caller")

dt <- dt1[Filter == "PASS" & !is.na(TRUTH.TP)]

# Overall
overall <- cor.test(dt$bases, dt$TRUTH.TP, method = "spearman")
cat("Overall: rho =", round(overall$estimate, 3), "p =", round(overall$p.value, 4), "\n")

# By variant type
for (vtype in c("SNP", "INDEL")) {
    sub <- dt[Type == vtype]
    ct <- cor.test(sub$bases, sub$TRUTH.TP, method = "spearman")
    cat(vtype, ": rho =", round(ct$estimate, 3), "p =", round(ct$p.value, 4), "\n")
}

ggsave("~/Desktop/roland/VarCallbench/plts/base_cor.pdf", gg, width = 6, height = 3.5)

dt2 <- perf[qual, on = .(sample = dataset_id)]
dt2[, tool := label_fn(tool)]
dt2[, sample := fcase(
    sample == "H526_dRNA_ONT", "H526-dRNA004",
    sample == "H526_bulk_PB",  "H526-MasSeq",
    sample == "H526_bulk_ONT", "H526-cDNAxR10",
    sample == "H211_dRNA_ONT", "H211-dRNA004",
    sample == "H211_bulk_PB",  "H211-MasSeq",
    sample == "H211_bulk_ONT", "H211-cDNAxR10",
    default = sample
)]
dt2[, tech := sub("^[^-]+-", "", sample)]
dt2 <- dt2[!(grepl("longcallR", tool) & Type =="INDEL")]
dt2 <- dt2[!(grepl("dRNA004|cDNA", tech) & 
               tool == "isoLASER" & Type == "INDEL")]
dt2[, Type := factor(Type, levels = c("SNP", "INDEL"), 
                    labels = c("SNV", "INDEL"))]


pp <- ggplot(dt2, aes((1-error_rate)*100, METRIC.F1_Score, shape = tool, 
                col = tech, group = tool)) +
    geom_point(size = 2, alpha = 0.7) +
    facet_grid2(Type ~ coverage, scales = "free") +
    theme_classic() +
    scale_color_manual(values = pal) +
    theme(
        panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.8
        ),
        strip.background = element_rect(
            fill = "white",
            color = "black",
            linewidth = 0.8),
        axis.line = element_line(color = "black", linewidth = 0.3),
        panel.spacing = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        legend.title = element_text(size = 11),
        aspect.ratio = 1
    ) +
    labs(col = "Chemistry", x = "1 - Error Rate (%)", y = "F1",
         shape = "Variant Caller")

p <- pp / gg + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold"))

dt <- dt2[Filter == "PASS" & coverage == 5]
ct_all <- cor.test(dt$mapped_reads, dt$TRUTH.TP, method = "spearman")


ggsave("~/Desktop/roland/VarCallbench/plts/base_qual_cor.pdf", p, width = 11, 
       height = 7)
