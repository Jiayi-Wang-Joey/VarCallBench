suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
    library(patchwork)
    library(scales)
})

dt <- fread("~/Desktop/roland/VarCallbench/out/happy_summary_collector/9938922d/happy_summary_merged.csv")
dt <- dt[Filter=="PASS"]
dt[, sample := fcase(
    sample == "H526_dRNA_ONT", "H526-dRNA004",
    sample == "H526_bulk_PB",  "H526-MasSeq",
    sample == "H526_bulk_ONT", "H526-cDNAxR10",
    sample == "H211_dRNA_ONT", "H211-dRNA004",
    sample == "H211_bulk_PB",  "H211-MasSeq",
    sample == "H211_bulk_ONT", "H211-cDNAxR10",
    default = sample
)]

label_fn <- function(x) fcase(
    x == "clair3_rna",   "Clair3-RNA",
    x == "deep_variant", "DeepVariant",
    x == "longcallR",    "longcallR",     
    x == "longcallR_nn", "longcallR-nn",
    x == "isolaser",     "isoLASER",
    default = x
)

dt[, tool := label_fn(tool)]
dt[, cell_line := sub("^([^-]+)-.*$", "\\1", sample)]
dt[, platform      := sub("^[^-]+-(.*)$", "\\1", sample)]
dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
cols <- c(
    "Clair3-RNA" = "#A6CEE3",
    "DeepVariant" = "#52AF43",
    "longcallR" = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER" = "#FDBF6F"
)
dt[, tool := factor(tool, levels = names(cols))]
td <- dt[!(grepl("longcallR", tool) & Type =="INDEL")]
td <- td[!(grepl("dRNA004|cDNA", platform) & tool == "isoLASER" & Type == "INDEL")]
td <- td[coverage == 5]

tt <- melt(
    td,
    id.vars = c("Type", "Filter", "aligner", "tool", "sample"),
    measure.vars = c("TRUTH.TP", "TRUTH.FN", "QUERY.FP"),
    variable.name = "Category",
    value.name = "Count"
)
tt[, Category := tstrsplit(Category, "\\.", keep = 2)]
tt[, Type := factor(Type, levels = c("SNP", "INDEL"))]

aes_number <- list(
    geom_col(stat = "identity"),
    geom_text(
        aes(label = Count),
        position = position_stack(vjust = 0.5),
        size = 1.5
    ),
    facet_grid2(Type ~ sample, scales="free_y"),
    theme_classic(),
    labs(y = "Total", x = "Variant caller"),
    scale_fill_brewer(palette = "Paired"),
    theme(
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor = element_blank(),
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
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
)

p1 <- ggplot(tt[grepl("H526", sample) ], aes(tool, Count, fill=Category)) +
    aes_number
ggsave("~/Desktop/roland/VarCallbench/plts/number_526.pdf", p1, width = 6, height = 4)

p2 <- ggplot(tt[grepl("H211", sample) ], aes(tool, Count, fill=Category)) +
    aes_number
ggsave("~/Desktop/roland/VarCallbench/plts/number_211.pdf", p2, width = 6, height = 4)
