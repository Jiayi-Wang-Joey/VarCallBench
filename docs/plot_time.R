suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
})

dt <- fread("/Volumes/jiayiwang/VarCallbench/out/performances.tsv")
dt <- dt[, .(
    s,
    max_rss,
    max_vms,
    max_uss,
    max_pss,
    io_in,
    io_out,
    mean_load,
    cpu_time,
    dataset,
    module
)]


vc <- c("clair3_rna", "deep_variant", "longcallR", "longcallR_nn")
dt[, dataset := fcase(
    dataset == "H526_dRNA_ONT", "H526-dRNA004",
    dataset == "H526_bulk_PB",  "H526-Kinnex",
    dataset == "H526_bulk_ONT", "H526-cDNAxR10",
    dataset == "H211_dRNA_ONT", "H211-dRNA004",
    dataset == "H211_bulk_PB",  "H211-Kinnex",
    dataset == "H211_bulk_ONT", "H211-cDNAxR10",
    default = dataset
)]

dt[, tech := sub("^[^-]+-", "", dataset)]
dt_sum <- dt[module %in% vc,
             .(
                 wall_time_min = mean(s) / 60,
                 cpu_time = mean(cpu_time) / 60,
                 memory_gb = mean(max_pss) / 1024
             ),
             by = .(tech, module)
]


dt_sum[, module := fifelse(module == "clair3_rna", "Clair3-RNA",
                           fifelse(module == "deep_variant", "DeepVariant",
                           fifelse(module == "longcallR_nn", "longcallR-nn",
                                           module)))]



cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928"
)

metrics <- c("wall_time_min", "cpu_time", "memory_gb")

td <- melt(
    dt_sum,
    id.vars = c("tech", "module"),
    measure.vars = metrics,
    variable.name = "metric",
    value.name = "value"
)
td[, metric := fcase(
    metric == "wall_time_min", "Elapsed time (min)",
    metric == "cpu_time", "CPU time (min)",
    metric == "memory_gb", "Memory (GB)"
)]

gg <- ggplot(td, aes(x = module, y = value, fill = module)) +
    geom_col() +
    facet_grid2(metric ~ tech, scales = "free") +
    theme_classic() +
    scale_fill_manual(values = cols) +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        legend.title = element_text(face = "bold")
    ) +
    labs(fill = "Variant caller", x = NULL, y = NULL)


ggsave("plts/time.pdf", gg, units = "cm", width = 15, height = 12)
