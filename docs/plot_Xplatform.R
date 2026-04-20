suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

dt <- fread("out/cross_platform_detail_summary.csv")
dt[, detail_class := factor(
    detail_class,
    levels = c(
        "cDNA_Kinnex_dRNA",
        "cDNA_Kinnex",
        "cDNA_dRNA",
        "Kinnex_dRNA",
        "cDNA_only",
        "Kinnex_only",
        "dRNA_only"
    )
)]

dt[, detail_label := fcase(
    detail_class == "cDNA_Kinnex_dRNA", "All 3",
    detail_class == "cDNA_Kinnex",      "cDNA + Kinnex",
    detail_class == "cDNA_dRNA",        "cDNA + dRNA",
    detail_class == "Kinnex_dRNA",      "Kinnex + dRNA",
    detail_class == "cDNA_only",        "cDNA only",
    detail_class == "Kinnex_only",      "Kinnex only",
    detail_class == "dRNA_only",        "dRNA only"
)]
dt[caller=="clair3_rna", caller:="Clair3-RNA"]
dt[caller=="deep_variant", caller:="DeepVariant"]
dt[, caller := factor(
    caller,
    levels = c("Clair3-RNA", "DeepVariant", "longcallR", "longcallR_nn")
)]
dt[,total:=sum(N), by = .(cell_line, caller)]
gg <- ggplot(dt, aes(x = caller, y = N, fill = detail_label)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(
        data = dt[detail_label == "All 3"],
        aes(
            y = 50000,
            label = scales::percent(proportion, accuracy = 0.1)
        ),
        size = 2.5
    ) +
    facet_grid(~cell_line) +
    labs(
        x = NULL,
        y = "Number of variants",
        fill = NULL, 
        label = NULL,
        size = NULL
    ) +
    theme_classic() +
    guides(size = "none") +
    scale_fill_brewer(palette = "Paired") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank()
    )

ggsave("plts/crossplatform.pdf", gg, width = 14, height = 8, units = "cm")
