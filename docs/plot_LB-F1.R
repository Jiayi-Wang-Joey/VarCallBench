suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
    library(patchwork)
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

truth_labels <- dt[Type=="SNP", .(TRUTH.TOTAL = TRUTH.TOTAL[1]),
                   by = .(coverage, cell_line, platform)]

aes <- list(
    theme_classic(),
    facet_grid2(cell_line ~ platform, scales = "free"),
    labs(
        x = "Coverage Cutoff (DP >= n)",
        y = "F1 Score",
        color = "Variant Caller",
    ),
    scale_color_manual(values = cols),
    scale_x_discrete(expand = expansion(add = c(0.8, 1))),
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
    )  
)

p1 <- ggplot(dt[Type=="SNP"], aes(coverage, METRIC.F1_Score,
                                  group = tool, col=tool)) +
    geom_point(alpha=0.8, size = 1.5) +
    geom_line(alpha=0.8, linewidth=0.8) +
    geom_text(data = truth_labels,
              aes(x = coverage, y = 1, label = paste0(TRUTH.TOTAL)),
              inherit.aes = FALSE,
              angle = 30, hjust = 0.5, vjust = 1.3,
              size = 2, color = "grey30") +
    aes 

td <- dt[Type=="INDEL"]
td <- td[!grepl("longcallR", tool)]
td <- td[!(grepl("dRNA004|cDNA", platform) & tool == "isoLASER")]
truth2 <- td[Type=="INDEL", .(TRUTH.TOTAL = TRUTH.TOTAL[1]),
                   by = .(coverage, cell_line, platform)]
p2 <- ggplot(td, aes(coverage, METRIC.F1_Score,
                                  group = tool, col=tool)) +
    geom_point(alpha=0.8, size = 1.5) +
    geom_line(alpha=0.8, linewidth=0.8) +
    geom_text(data = truth2,
              aes(x = coverage, y = 1, label = paste0(TRUTH.TOTAL)),
              inherit.aes = FALSE,
              angle = 30, hjust = 0.5, vjust = 1.3,
              size = 2, color = "grey30") +
    aes +
    theme(legend.position = "none")

gg <- p1 | p2 +
    plot_layout(guides = "collect")
ggsave("~/Desktop/roland/VarCallbench/plts/F1.pdf", gg, width = 12, height = 3.5)
