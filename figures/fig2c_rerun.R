# ============================================================
# Fig. 2 panel c — DEG jitter plot (re-rendered)
#
# Part of the code for the Data Descriptor
#   "Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq data from
#    heat-stressed lactating Holstein cows"
#
# Produces: Fig2_panel_c_DEG_jitter.pdf (8x8 cm)
# Re-rendered to fix glyph-advance quantisation; set.seed(42) added so the jitter is reproducible.
#
# Inputs are expected under DATA_DIR (default "data/"), or set HS_DATA_DIR.
# Archived at figshare https://doi.org/10.6084/m9.figshare.32993768
#
# The original working comments are in Chinese and are left unchanged;
# an English line has been added below each of them.
# ============================================================

## Fig 2 panel c (DEG jitter) 重出 —— 修字距
## 数据: <DATA_DIR>/DEG_results_all_tissues.rds (冻结)
##       DEG 口径 = 裸 P<0.05 & |log2FC|>=1, 与脚本原样一致 -> 计数不变
## ⚠ 原脚本【没有 set.seed】, geom_jitter 的横向偏移每跑一次都不同 =>
##   已发表那张的点云位置无法精确复现。此处补 set.seed(42): y值(log2FC)与计数
##   一个不变, 只有抖动噪声的随机横向位置与原图不同; 且从此可复现。
## 字体: 本图 base_family="Helvetica" -> 本机解析到 Nimbus Sans, showtext 指同一文件复刻原样。
## [EN] Re-render of Fig. 2 panel c (DEG jitter) — letter-spacing fix.
## Data: <DATA_DIR>/DEG_results_all_tissues.rds (frozen).
## DEG criteria = nominal P < 0.05 & |log2FC| >= 1, identical to the archived script -> counts unchanged.
## The original script had NO set.seed, so geom_jitter offsets differed on every run and the published
## point cloud cannot be reproduced exactly. set.seed(42) is added here: log2FC values and counts are
## unchanged; only the random horizontal jitter differs from the published panel, and it is now reproducible.
## Font: this panel uses base_family="Helvetica", which resolves to Nimbus Sans on this machine;
## showtext is pointed at the same file to reproduce the original exactly.
# .libPaths(c("<site-specific R library path>", .libPaths()))
suppressMessages({library(ggplot2); library(dplyr)
                  library(showtext); library(sysfonts)})
font_add("Helvetica", regular = "/usr/share/fonts/opentype/urw-base35/NimbusSans-Regular.otf")
showtext_auto()
showtext_opts(dpi = 600)
set.seed(42)
DATA_DIR <- Sys.getenv("HS_DATA_DIR", "data")
outdir   <- Sys.getenv("HS_OUT_DIR", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# Fig 2c — 纵向jitter散点图 (DEG)
# 输出: Fig2_panel_c_DEG_jitter.pdf (8×8cm)
# [EN] Fig. 2c — vertical jitter plot of DEG counts
# Output: Fig2_panel_c_DEG_jitter.pdf (8x8 cm)
# ############################################################

bulk_deg <- readRDS(file.path(DATA_DIR, "DEG_results_all_tissues.rds"))
results_combined_clean <- na.omit(bulk_deg)

tissue_labels <- c("fat_HS_vs_PF" = "Adipose", "liver_HS_vs_PF" = "Liver",
                   "mammary gland_HS_vs_PF" = "Mammary\ngland",
                   "muscle_HS_vs_PF" = "Muscle", "rumen_HS_vs_PF" = "Rumen")
results_combined_clean$Tissue <- tissue_labels[results_combined_clean$Comparison]
results_combined_clean$Tissue <- factor(results_combined_clean$Tissue,
  levels = c("Rumen", "Liver", "Mammary\ngland", "Adipose", "Muscle"))

results_combined_clean$Category <- "NS"
results_combined_clean$Category[results_combined_clean$pvalue < 0.05 &
  results_combined_clean$log2FoldChange >= 1] <- "Up"
results_combined_clean$Category[results_combined_clean$pvalue < 0.05 &
  results_combined_clean$log2FoldChange <= -1] <- "Down"

gray_data <- results_combined_clean[results_combined_clean$Category == "NS", ]
up_data_plot <- results_combined_clean[results_combined_clean$Category == "Up", ]
down_data_plot <- results_combined_clean[results_combined_clean$Category == "Down", ]

deg_n <- results_combined_clean %>%
  filter(Category != "NS") %>%
  group_by(Tissue, Category) %>%
  summarise(n = n(), .groups = "drop")

up_labels <- merge(deg_n[deg_n$Category == "Up", ],
  up_data_plot %>% group_by(Tissue) %>%
    summarise(ypos = max(log2FoldChange) + 1.5, .groups = "drop"), by = "Tissue")
down_labels <- merge(deg_n[deg_n$Category == "Down", ],
  down_data_plot %>% group_by(Tissue) %>%
    summarise(ypos = min(log2FoldChange) - 1.5, .groups = "drop"), by = "Tissue")

y_ceil <- max(up_data_plot$log2FoldChange, na.rm = TRUE) + 5
y_floor <- min(down_data_plot$log2FoldChange, na.rm = TRUE) - 5

bg_rects <- data.frame(
  xmin = seq(0.5, 4.5, by = 1), xmax = seq(1.5, 5.5, by = 1),
  ymin = y_floor, ymax = y_ceil,
  fill_col = c("grey90", "grey96", "grey90", "grey96", "grey90")
)

up_labels2 <- up_labels; up_labels2$ypos <- y_ceil - 1.5
down_labels2 <- down_labels; down_labels2$ypos <- y_floor + 1.5

p_d <- ggplot() +
  geom_rect(data = bg_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_rects$fill_col) +
  geom_jitter(data = gray_data, aes(x = Tissue, y = log2FoldChange),
              color = "#D0D0D0", size = 0.15, width = 0.35, alpha = 0.1) +
  geom_jitter(data = up_data_plot, aes(x = Tissue, y = log2FoldChange),
              color = "#D73027", size = 0.35, width = 0.35, alpha = 0.5) +
  geom_jitter(data = down_data_plot, aes(x = Tissue, y = log2FoldChange),
              color = "#4575B4", size = 0.35, width = 0.35, alpha = 0.5) +
  geom_text(data = up_labels2, aes(x = Tissue, y = ypos, label = n),
            color = "#D73027", size = 1.8, fontface = "bold", family = "Helvetica") +
  geom_text(data = down_labels2, aes(x = Tissue, y = ypos, label = n),
            color = "#4575B4", size = 1.8, fontface = "bold", family = "Helvetica") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "black") +
  labs(x = "", y = expression(Average~Log[2]~Fold~Change)) +
  scale_y_continuous(limits = c(y_floor, y_ceil), expand = c(0, 0)) +
  theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 7),   ## 归档脚本写的是 element_blank(), 但已发表面板上有纵轴标题 -> 以原版为准
    axis.text.y = element_text(size = 6, face = "bold", color = "black"),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.line.y = element_line(color = "black", linewidth = 0.3),
    axis.ticks.y = element_line(color = "black", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

ggsave(file.path(outdir, "Fig2_panel_c_DEG_jitter.pdf"), p_d,
       width = 8, height = 8, units = "cm", dpi = 600, device = cairo_pdf)

cat("FIG2C DONE\n")
