suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
    library(patchwork)
})

dt <- fread("~/Desktop/roland/VarCallbench/out/callable_bases.csv")
dt[, dataset := fcase(
    dataset == "H526_dRNA_ONT", "H526-dRNA004",
    dataset == "H526_bulk_PB",  "H526-MasSeq",
    dataset == "H526_bulk_ONT", "H526-cDNAxR10",
    dataset == "H211_dRNA_ONT", "H211-dRNA004",
    dataset == "H211_bulk_PB",  "H211-MasSeq",
    dataset == "H211_bulk_ONT", "H211-cDNAxR10",
    default = dataset
)]
dt[, cell_line := sub("^([^-]+)-.*$", "\\1", dataset)]
dt[, platform      := sub("^[^-]+-(.*)$", "\\1", dataset)]
dt[,coverage:=factor(coverage, levels = unique(sort(dt$coverage)))]
aes <- list(
    theme_classic(),
    facet_grid2(cell_line ~ platform, scales = "free"),
    labs(
        x = "Coverage Cutoff (DP >= n)",
        y = "Callable bases"
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


gg <- ggplot(dt, aes(coverage, bases)) +
    geom_bar(stat = "identity") + 
    aes

ggsave("~/Desktop/roland/VarCallbench/plts/callable_bases.pdf", gg, width = 6, height = 3.5)

