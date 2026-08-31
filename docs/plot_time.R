suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
})

dt <- fread("~/Desktop/roland/VarCallbench/out/performances.tsv")

dt <- dt[, .(s, max_rss, max_vms, max_uss, max_pss, io_in, io_out, mean_load, cpu_time, dataset, module)]

vc <- c("clair3_rna", "deep_variant", "longcallR", "longcallR_nn", "isolaser", "isolaser_annotate")

dt[, dataset := fcase(
    dataset == "H526_dRNA_ONT", "H526-dRNA004",
    dataset == "H526_bulk_PB",  "H526-MasSeq",
    dataset == "H526_bulk_ONT", "H526-cDNAxR10",
    dataset == "H211_dRNA_ONT", "H211-dRNA004",
    dataset == "H211_bulk_PB",  "H211-MasSeq",
    dataset == "H211_bulk_ONT", "H211-cDNAxR10",
    default = dataset
)]
dt <- dt[grepl("H211", dataset),]
dt[, tech := sub("^[^-]+-", "", dataset)]
dt <- dt[module %in% vc]

label_fn <- function(x) fcase(
    x == "clair3_rna",        "Clair3-RNA",
    x == "deep_variant",      "DeepVariant",
    x == "longcallR_nn",      "longcallR-nn",
    x == "isolaser",          "isoLASER",
    x == "isolaser_annotate", "Annotation",
    default = x
)

## Elapsed/CPU time are additive across sequential stages -> keep isolaser and
## isolaser_annotate as separate stacked segments
dt_time <- dt[, .(
    wall_time_min = mean(s) / 60,
    cpu_time      = mean(cpu_time) / 60
), by = .(tech, module)]
dt_time[, caller := fcase(module %in% c("isolaser", "isolaser_annotate"), "isolaser", default = module)]
dt_time[, step := module]

## Peak RSS/PSS is NOT additive across sequential processes -> collapse
## isolaser + isolaser_annotate into a single isoLASER bar via max(), no
## separate Annotation segment
dt_mem_stage <- dt[, .(mem = mean(max_pss) / 1024), by = .(tech, module)]
dt_mem_stage[, caller := fcase(module %in% c("isolaser", "isolaser_annotate"), "isolaser", default = module)]
dt_mem <- dt_mem_stage[, .(memory_gb = max(mem)), by = .(tech, caller)]
dt_mem[, step := caller]

dt_time[, `:=`(module = label_fn(module), caller = label_fn(caller), step = label_fn(step))]
dt_mem[,  `:=`(caller = label_fn(caller), step = label_fn(step))]

td_time <- melt(dt_time, id.vars = c("tech", "caller", "step"),
                measure.vars = c("wall_time_min", "cpu_time"),
                variable.name = "metric", value.name = "value")

td_mem <- melt(dt_mem, id.vars = c("tech", "caller", "step"),
               measure.vars = "memory_gb",
               variable.name = "metric", value.name = "value")

td <- rbind(td_time, td_mem)

td[, metric := fcase(
    metric == "wall_time_min", "Elapsed time (min)",
    metric == "cpu_time",      "CPU time (min)",
    metric == "memory_gb",     "Memory (GB)"
)]

td[, caller := factor(caller, levels = c("Clair3-RNA", "DeepVariant", "longcallR", "longcallR-nn", "isoLASER"))]
td[, step   := factor(step,   levels = c("Clair3-RNA", "DeepVariant", "longcallR", "longcallR-nn", "isoLASER", "Annotation"))]
td[, metric := factor(metric, levels = c("CPU time (min)", "Elapsed time (min)", "Memory (GB)"))]

cols <- c(
    "Clair3-RNA" = "#A6CEE3",
    "DeepVariant" = "#52AF43",
    "longcallR" = "#B294C7",
    "longcallR-nn" = "#B15928",
    "isoLASER" = "#FDBF6F",
    "Realign & Annotation" = "grey90"
)

gg <- ggplot(td, aes(x = caller, y = value, fill = step)) +
    geom_col(width = 0.8) +
    facet_grid2(metric ~ tech, scales = "free") +
    scale_fill_manual(values = cols, drop = FALSE) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_text(face = "bold")
    ) +
    labs(fill = "Variant caller", x = NULL, y = NULL)

gg


ggsave("plts/time_h211.pdf", gg, units = "cm", width = 15, height = 12)


td <- fread("~/Desktop/roland/VarCallbench/out/normalized_runtime_v2.csv")
td <- td[cell_line!="Average" & cell_line == "H211"]
td[, runtime_per_M_reads_min  := runtime_per_M_reads_s  / 60]
td[, cputime_per_M_reads_min  := cputime_per_M_reads_s  / 60]
td[, memory_per_M_reads_GB := memory_per_M_reads_MB / 1024]

metrics <- c("runtime_per_M_reads_min", "memory_per_M_reads_GB", "cputime_per_M_reads_min")

tdm <- melt(
    td,
    id.vars = c("caller", "platform"),
    measure.vars = metrics,
    variable.name = "metric",
    value.name = "value"
)

tdm[, metric := fcase(
    metric == "runtime_per_M_reads_min", "Elpased time (min)",
    metric == "memory_per_M_reads_GB",   "Memory (GB)",
    metric == "cputime_per_M_reads_min", "CPU time (min)"
)]
tdm[, platform := fcase(
    platform == "cDNA",  "cDNAxR10",
    platform == "dRNA",  "dRNA004",
    default = platform
)]
tdm[, caller := factor(
    caller,
    levels = c("Clair3-RNA", "DeepVariant", "longcallR", "longcallR-nn", "isoLASER")
)]
tdm[, metric := factor(
    metric,
    levels = c("CPU time (min)", "Elpased time (min)",  "Memory (GB)")
)]
gg <- ggplot(tdm, aes(x = caller, y = value, fill = caller)) +
    geom_col(width = 0.8) +
    facet_grid2(metric ~ platform, scales = "free") +
    scale_fill_manual(values = cols, drop = FALSE) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_text(face = "bold")
    ) +
    labs(fill = "Variant caller", x = NULL, y = NULL)

gg
ggsave("~/Desktop/roland/VarCallbench/plts/normalized_time_h211.pdf", gg, units = "cm", width = 15, height = 12)

td <- fread("~/Desktop/roland/VarCallbench/out/normalized_runtime_v2.csv")
td <- td[cell_line!="Average" & cell_line == "H211"]
td[, runtime_per_Gb_min  := runtime_per_Gb_s  / 60]
td[, cputime_per_Gb_min  := cputime_per_Gb_s  / 60]
td[, memory_per_Gb_GB := memory_per_Gb_MB / 1024]

metrics <- c("runtime_per_Gb_min", "memory_per_Gb_GB", "cputime_per_Gb_min")

tdm <- melt(
    td,
    id.vars = c("caller", "platform"),
    measure.vars = metrics,
    variable.name = "metric",
    value.name = "value"
)

tdm[, metric := fcase(
    metric == "runtime_per_Gb_min", "Elpased time (min/Gb)",
    metric == "memory_per_Gb_GB",   "Memory (GB/Gb)",
    metric == "cputime_per_Gb_min", "CPU time (min/Gb)"
)]
tdm[, platform := fcase(
    platform == "cDNA",  "cDNAxR10",
    platform == "dRNA",  "dRNA004",
    default = platform
)]
tdm[, caller := factor(
    caller,
    levels = c("Clair3-RNA", "DeepVariant", "longcallR", "longcallR-nn", "isoLASER")
)]
tdm[, metric := factor(
    metric,
    levels = c("CPU time (min/Gb)", "Elpased time (min/Gb)",  "Memory (GB/Gb)")
)]
gg <- ggplot(tdm, aes(x = caller, y = value, fill = caller)) +
    geom_col(width = 0.8) +
    facet_grid2(metric ~ platform, scales = "free") +
    scale_fill_manual(values = cols, drop = FALSE) +
    theme_classic() +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.title = element_text(face = "bold")
    ) +
    labs(fill = "Variant caller", x = NULL, y = NULL)

gg
ggsave("~/Desktop/roland/VarCallbench/plts/normalized_time_h211_per_base.pdf", gg, units = "cm", width = 15, height = 12)
