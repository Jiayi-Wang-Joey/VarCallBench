suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
    library(ggh4x)
    library(tidytext)
})
setwd("/Volumes/jiayiwang/VarCallbench/")
dt <- fread("out/alignment_qc_collector/alignment_qc_merged.csv")
dt[, dataset_id := fcase(
    dataset_id == "H526_dRNA_ONT", "H526-dRNA004",
    dataset_id == "H526_bulk_PB",  "H526-Kinnex",
    dataset_id == "H526_bulk_ONT", "H526-cDNAxR10",
    dataset_id == "H211_dRNA_ONT", "H211-dRNA004",
    dataset_id == "H211_bulk_PB",  "H211-Kinnex",
    dataset_id == "H211_bulk_ONT", "H211-cDNAxR10",
    default = dataset_id
)]

dt[, cell_line := sub("-.*", "", dataset_id)]
dt[, tech := sub("^[^-]+-", "", dataset_id)]
metrics <- c(
    "total_sequences",
    "mapped_reads",
    "mapping_rate_pct",
    "average_length",
    "average_quality",
    "error_rate"
)

pal <- c(
    "dRNA004" = "#fec44f",
    "cDNAxR10" = "#d95f0e",
    "Kinnex" = "#54278f"
)

dt_long <- melt(
    dt,
    id.vars = c("dataset_id", "cell_line", "tech"),
    measure.vars = metrics,
    variable.name = "metric",
    value.name = "value"
)

dt_long[, value_scaled := fcase(
    metric %in% c("total_sequences", "mapped_reads"),
    value / 1e6,
    metric == "error_rate",
    value * 100,
    default = value
)]


dt_long[, metric := factor(metric,
                           levels = metrics,
                           labels = c("Total sequences", "Mapped reads", "Mapping rate", 
                                      "Read length", "Average quality", "Error rate")
)]


aes <- list(
    geom_bar(stat = "identity", position = position_dodge()),
    # scale_fill_identity(guide = "legend",
    #                     breaks = pal,
    #                     labels = names(pal)),
    scale_fill_manual(values=pal),
    facet_grid2(cell_line ~ metric, scales = "free", 
                axes="all", independent="all") ,
    theme_minimal(),
    theme(
        panel.border = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_blank()
    )
)

dt_long[, facet_group := interaction(metric, cell_line, drop = TRUE)]
gg <- ggplot(dt_long, aes(reorder_within(tech, value, facet_group), 
                          value_scaled, fill=tech)) + 
    aes +
    labs(x = NULL, y = "Value", fill = "Technology")

ggsave("plts/qual.pdf", gg, width=28, height=8, units="cm")
