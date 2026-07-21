# ============================================================
# Fig. 3 — pathway enrichment bubble plot (re-rendered)
#
# Part of the code for the Data Descriptor
#   "Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq data from
#    heat-stressed lactating Holstein cows"
#
# Produces: Fig3_ab_complete.pdf (18x10 cm)
# Font rendering only: enrichment is NOT recomputed, so there is no version drift and no randomness.
#
# Inputs are expected under DATA_DIR (default "data/"), or set HS_DATA_DIR.
# Archived at figshare https://doi.org/10.6084/m9.figshare.32993768
#
# The original working comments are in Chinese and are left unchanged;
# an English line has been added below each of them.
# ============================================================

## Fig 3 重出 —— 只修字体渲染, 图的内容一个点不变
## 数据来自冻结的 CSV(<DATA_DIR>/*_final.csv, 2025-02 至今没动过),
## 不重跑富集分析 -> 无版本漂移、无随机性。
##
## 病根: 原 Fig3 那场 R 会话没开 showtext, cairo 自己排版把每个字形推进吸附到整数 pt
##       (实测 93% 的推进落在整数上), 字母间距忽松忽紧。Fig5/Fig6 那几场开了 showtext
##       所以没事 —— 主脚本由多次独立会话合并而成, 各场字体设置不一, 故只有部分图中招。
## 注意: 本图 base_family 是 "Helvetica"(不是 Arial), 在本机 fontconfig 解析到
##       Nimbus Sans(URW 的 Helvetica 克隆), 故 showtext 必须指同一个文件才能复刻原样。
## [EN] Re-render of Fig. 3 — font rendering only; not a single data point changes.
## Data come from the frozen CSVs (<DATA_DIR>/*_final.csv, untouched since 2025-02);
## the enrichment analysis is NOT re-run -> no package-version drift, no randomness.
## Cause: the original Fig. 3 session did not enable showtext, so cairo laid out the text itself and
## quantised every glyph advance to whole points (93% measured on integers), giving uneven letter spacing.
## The Fig. 5 / Fig. 6 sessions had showtext enabled and were unaffected — the master script was
## assembled from several independent sessions with differing font settings, so only some figures were hit.
## NOTE: base_family here is "Helvetica" (not Arial); on this machine fontconfig resolves it to
## Nimbus Sans (URW Helvetica clone), so showtext must point at that same file to reproduce the original.
# .libPaths(c("<site-specific R library path>", .libPaths()))
suppressMessages({library(ggplot2); library(dplyr); library(stringr); library(patchwork)
                  library(showtext); library(sysfonts)})
font_add("Helvetica", regular = "/usr/share/fonts/opentype/urw-base35/NimbusSans-Regular.otf")
showtext_auto()
showtext_opts(dpi = 600)
DATA_DIR <- Sys.getenv("HS_DATA_DIR", "data")
outdir   <- Sys.getenv("HS_OUT_DIR", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Fig 3 — 通路富集气泡图 (Nature红蓝配色, V7最终版)
# 输出: Fig3_ab_complete.pdf (18×10cm)
# [EN] Fig. 3 — pathway enrichment bubble plot (final v7)
# Output: Fig3_ab_complete.pdf (18x10 cm)
# ############################################################

library(stringr)
library(patchwork)

up_data <- read.csv(file.path(DATA_DIR, "up_data_final.csv"))
down_data <- read.csv(file.path(DATA_DIR, "down_data_final.csv"))

fix_organ <- function(x) {
  x <- gsub("mammary_gland", "Mammary gland", x)
  x <- gsub("^fat$", "Adipose", x)
  x <- gsub("^liver$", "Liver", x)
  x <- gsub("^rumen$", "Rumen", x)
  x <- gsub("^muscle$", "Muscle", x)
  return(x)
}
up_data$Organ <- fix_organ(up_data$Organ)
down_data$Organ <- fix_organ(down_data$Organ)
up_data$Organ <- factor(up_data$Organ,
  levels = c("Rumen", "Liver", "Mammary gland", "Adipose", "Muscle"))
down_data$Organ <- factor(down_data$Organ,
  levels = c("Rumen", "Liver", "Mammary gland", "Adipose", "Muscle"))

redundant_terms <- c(
  "mitotic sister chromatid segregation", "nuclear chromosome segregation",
  "sister chromatid segregation", "mitotic cell cycle checkpoint signaling",
  "regulation of cell cycle process", "regulation of cell cycle phase transition",
  "regulation of cell cycle", "cell cycle phase transition", "cell cycle checkpoint signaling",
  "mitotic nuclear division", "negative regulation of mitotic cell cycle phase transition",
  "negative regulation of mitotic cell cycle", "metaphase chromosome alignment",
  "regulation of mitotic cell cycle phase transition", "regulation of mitotic cell cycle",
  "negative regulation of cell cycle phase transition", "mitotic cell cycle phase transition",
  "negative regulation of cell cycle process", "negative regulation of cell cycle",
  "establishment of chromosome localization", "chromosome localization",
  "negative regulation of chromosome organization", "organelle fission",
  "regulation of sister chromatid segregation", "regulation of chromosome segregation",
  "regulation of G2/M transition of mitotic cell cycle", "regulation of cell cycle G2/M phase transition",
  "regulation of mitotic nuclear division", "G2/M transition of mitotic cell cycle",
  "cell cycle G2/M phase transition", "cytoskeleton-dependent cytokinesis",
  "mitotic cytokinesis", "cell morphogenesis involved in neuron differentiation",
  "neuron projection morphogenesis", "plasma membrane bounded cell projection morphogenesis",
  "cell projection morphogenesis", "neuron projection guidance",
  "generation of neurons", "neuron projection development",
  "negative regulation of cellular component organization",
  "negative regulation of organelle organization",
  "regulation of chromosome organization", "regulation of nuclear division",
  "monocarboxylic acid metabolic process", "actomyosin structure organization"
)

up_clean <- up_data[!up_data$Description %in% redundant_terms, ]
down_clean <- down_data[!down_data$Description %in% redundant_terms, ]

up_top <- up_clean %>% group_by(Description) %>%
  summarise(n_organs = n_distinct(Organ), min_p = min(pvalue), .groups = "drop") %>%
  arrange(desc(n_organs), min_p) %>% head(20)
down_top <- down_clean %>% group_by(Description) %>%
  summarise(n_organs = n_distinct(Organ), min_p = min(pvalue), .groups = "drop") %>%
  arrange(desc(n_organs), min_p) %>% head(20)

up_plot <- up_clean[up_clean$Description %in% up_top$Description, ]
down_plot <- down_clean[down_clean$Description %in% down_top$Description, ]
up_order <- up_top$Description[order(up_top$n_organs, -up_top$min_p)]
down_order <- down_top$Description[order(down_top$n_organs, -down_top$min_p)]
up_plot$Description <- factor(up_plot$Description, levels = up_order)
down_plot$Description <- factor(down_plot$Description, levels = down_order)

up_plot_dedup <- up_plot %>% group_by(Description, Organ) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup()
down_plot_dedup <- down_plot %>% group_by(Description, Organ) %>%
  slice_min(pvalue, n = 1, with_ties = FALSE) %>% ungroup()

beautify_term <- function(x) {
  x <- gsub(" metabolic process$", " metabolism", x)
  x <- gsub(" biosynthetic process$", " biosynthesis", x)
  x <- gsub(" catabolic process$", " catabolism", x)
  x <- gsub(" signaling pathway$", " signaling", x)
  x <- gsub("^regulation of ", "Reg. of ", x)
  x <- gsub("^negative regulation of ", "Neg. reg. of ", x)
  x <- gsub("^positive regulation of ", "Pos. reg. of ", x)
  x <- gsub(" process$", "", x)
  x <- gsub(" activity$", "", x)
  x <- gsub("anatomical structure formation involved in morphogenesis",
             "Anatomical structure morphogenesis", x)
  x <- gsub("purine-containing compound metabolic process",
             "Purine compound metabolism", x)
  x <- gsub("secondary alcohol metabolic process",
             "Secondary alcohol metabolism", x)
  x <- gsub("organic hydroxy compound metabolic process",
             "Organic hydroxy compound metabolism", x)
  x <- gsub("regulation of multicellular organismal development",
             "Reg. of multicellular development", x)
  x <- gsub("mitotic DNA integrity checkpoint signaling",
             "Mitotic DNA checkpoint signaling", x)
  x <- gsub("DNA integrity checkpoint signaling",
             "DNA checkpoint signaling", x)
  x <- gsub("Cytoskeleton in muscle cells", "Muscle cytoskeleton", x)
  x <- gsub("carboxylic acid metabolic process", "Carboxylic acid metabolism", x)
  x <- gsub("organic acid metabolic process", "Organic acid metabolism", x)
  x <- gsub("oxoacid metabolic process", "Oxoacid metabolism", x)
  x <- gsub("^([a-z])", "\\U\\1", x, perl = TRUE)
  return(x)
}

up_plot_dedup$Term_clean <- beautify_term(as.character(up_plot_dedup$Description))
down_plot_dedup$Term_clean <- beautify_term(as.character(down_plot_dedup$Description))
up_levels_clean <- unique(beautify_term(levels(up_plot_dedup$Description)))
down_levels_clean <- unique(beautify_term(levels(down_plot_dedup$Description)))
up_plot_dedup$Term_clean <- factor(up_plot_dedup$Term_clean, levels = up_levels_clean)
down_plot_dedup$Term_clean <- factor(down_plot_dedup$Term_clean, levels = down_levels_clean)

nature_compact_v2 <- theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5, color = "black"),
    axis.text.y = element_text(size = 5.5, color = "black"),
    axis.title = element_text(size = 7),
    axis.line = element_line(colour = "black", linewidth = 0.3),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

p_up <- ggplot(up_plot_dedup, aes(x = Organ, y = Term_clean, size = Count, color = pvalue)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#B2182B", high = "#FDDBC7", name = "P-value") +
  scale_size_continuous(range = c(1.5, 6), name = "Gene\nCount", breaks = c(5, 10, 20, 30)) +
  scale_x_discrete(expand = c(0.08, 0.08)) +
  scale_y_discrete(expand = c(0.03, 0.03)) +
  coord_cartesian(clip = "off") +
  ggtitle("Up-regulated") + labs(x = "", y = "") +
  nature_compact_v2

p_down <- ggplot(down_plot_dedup, aes(x = Organ, y = Term_clean, size = Count, color = pvalue)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#2166AC", high = "#D1E5F0", name = "P-value") +
  scale_size_continuous(range = c(1.5, 6), name = "Gene\nCount", breaks = c(5, 10, 20, 30)) +
  scale_x_discrete(expand = c(0.08, 0.08)) +
  scale_y_discrete(expand = c(0.03, 0.03)) +
  coord_cartesian(clip = "off") +
  ggtitle("Down-regulated") + labs(x = "", y = "") +
  nature_compact_v2

fig3 <- p_up | p_down
ggsave(file.path(outdir, "Fig3_ab_complete.pdf"), fig3,
       width = 18, height = 10, units = "cm", dpi = 600, device = cairo_pdf)

cat("FIG3 DONE\n")
