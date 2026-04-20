suppressPackageStartupMessages({
    library(ggplot2)
    library(ggh4x)
    library(data.table)
    library(dplyr)
})

setwd("/Volumes/jiayiwang/VarCallbench/")
dt <- fread("out/somatic_detection_collector/somatic_detection_merged.csv")
dt[, dataset := fcase(
    dataset == "H526_dRNA_ONT", "H526-dRNA004",
    dataset == "H526_bulk_PB",  "H526-Kinnex",
    dataset == "H526_bulk_ONT", "H526-cDNAxR10",
    dataset == "H211_dRNA_ONT", "H211-dRNA004",
    dataset == "H211_bulk_PB",  "H211-Kinnex",
    dataset == "H211_bulk_ONT", "H211-cDNAxR10",
    default = dataset
)]


cols <- c(
    "Clair3-RNA"   = "#A6CEE3",
    "DeepVariant"  = "#52AF43",
    "GATK"         = "#F06C45",
    "longcallR"    = "#B294C7",
    "longcallR-nn" = "#B15928"
)



gg <- ggplot(dt, aes(dataset, detection_rate, fill = caller)) +
    geom_bar(stat="identity", position = position_dodge(), alpha = 0.9) + 
    theme_classic() +
    scale_fill_manual(values = cols) +
    labs(x="Dataset", y="Somatic variant detection rate", fill="Variant caller")

ggsave("plts/detection_rate.pdf", gg, width = 22, height = 10, units = "cm")
