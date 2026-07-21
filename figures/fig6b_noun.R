# ============================================================
# Fig. 6 panel b — module x cell-type heatmap (re-rendered)
#
# Part of the code for the Data Descriptor
#   "Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq data from
#    heat-stressed lactating Holstein cows"
#
# Produces: Fig6_panel_b_module_heatmap.pdf (10x6 cm)
# showtext enabled; the unassigned cell type is filtered out on the hME table (not on the object).
#
# Inputs are expected under DATA_DIR (default "data/"), or set HS_DATA_DIR.
# Archived at figshare https://doi.org/10.6084/m9.figshare.32993768
#
# The original working comments are in Chinese and are left unchanged;
# an English line has been added below each of them.
# ============================================================

# .libPaths(c("<site-specific R library path>", .libPaths()))
suppressMessages({library(Seurat); library(hdWGCNA); library(ggplot2); library(reshape2)
                  library(showtext); library(sysfonts)})

## ⚠ 必须开 showtext: 原版面板(2026-04-10)就是这么出的 —— 文字全轮廓化、0 可抽字符。
## 不开时 cairo 自己排版, 会把每个字形推进吸附到整数 pt, 字母间距忽松忽紧。
## [EN] showtext MUST be enabled: the original panel (2026-04-10) was produced this way — text fully
## outlined, 0 extractable characters. Without it cairo quantises glyph advances and spacing goes uneven.
MSF <- "/usr/share/fonts/truetype/msttcorefonts/"
font_add("Arial", regular = paste0(MSF, "Arial.ttf"),
         bold = paste0(MSF, "Arial_Bold.ttf"),
         italic = paste0(MSF, "Arial_Italic.ttf"))
showtext_auto()
showtext_opts(dpi = 600)   # 必须与 ggsave(dpi=600) 一致, 否则字号会错
# [EN] dpi here MUST match ggsave(dpi=) below, otherwise font sizes come out wrong

DATA_DIR <- Sys.getenv("HS_DATA_DIR", "data")
outdir   <- Sys.getenv("HS_OUT_DIR", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
obj2 <- readRDS(file.path(DATA_DIR, "scRNA", "integrated_25cl_hdwgcna.rds"))
cat("loaded obj2:", ncol(obj2), "cells\n")

## ---- HSCs -> Unassigned ----
obj2$celltype <- as.character(obj2$celltype)
obj2$celltype[obj2$celltype == "HSCs"] <- "Unassigned"
## 定案(2026-07-17): 热图删掉 Unassigned 这一行 —— 只展示匹配到的细胞类型。
## ⚠ 不能在这里 subset 对象: GetMEs() 取的是 hdWGCNA 存档的模块特征值, 行数仍是 82013,
##   砍了对象再 hMEs$celltype <- meta 就会长度对不齐报错。要在 hMEs 表上滤(见下)。
## [EN] Decision (2026-07-17): the unassigned row is removed from the heatmap — matched cell types only.
## WARNING: do NOT subset the object here. GetMEs() returns the module eigengenes stored by hdWGCNA,
## which still have 82,013 rows; subsetting the object first makes the meta assignment length-mismatch.
## Filter on the hME table instead (see below).
cat("celltype AFTER relabel:\n"); print(table(obj2$celltype))

## ---- Fig 6b: module-celltype heatmap (原样, 只换标签) ----
## [EN] ---- Fig. 6b: module x cell-type heatmap (unchanged; labels only) ----
hMEs <- GetMEs(obj2, harmonized = TRUE)
meta_sc <- obj2@meta.data
hMEs$celltype <- meta_sc$celltype
hMEs <- hMEs[, !grepl("grey", colnames(hMEs)) | colnames(hMEs) == "celltype"]
cat("滤前:", nrow(hMEs), "行\n")
hMEs <- hMEs[hMEs$celltype != "Unassigned", ]
cat("滤掉 Unassigned 后:", nrow(hMEs), "行 (应为 81903)\n")

mean_hME <- aggregate(. ~ celltype, data = hMEs, FUN = mean)
rownames(mean_hME) <- mean_hME$celltype
mean_hME$celltype <- NULL
colnames(mean_hME) <- gsub("^ME", "", colnames(mean_hME))

mod_order_6b <- c("black","blue","brown","green","red","turquoise","yellow")
mod_order_6b <- mod_order_6b[mod_order_6b %in% colnames(mean_hME)]
mean_hME <- mean_hME[, mod_order_6b]

ct_order <- c("T cells","B cells","Monocytes","NK cells",
              "Plasma cells","Dendritic cells")
ct_order <- ct_order[ct_order %in% rownames(mean_hME)]
mean_hME <- mean_hME[ct_order, ]
cat("heatmap rows:", paste(rownames(mean_hME), collapse=" | "), "\n")
cat("heatmap cols:", paste(colnames(mean_hME), collapse=" | "), "\n")

cat("\n=== mean hME 矩阵 (行=celltype, 列=module) ===\n"); print(round(mean_hME,3))
cat("\n=== 每个 celltype 的最高模块 ===\n")
for (r in rownames(mean_hME)) { v <- unlist(mean_hME[r,]); cat(sprintf("%-16s top=%-10s (%.3f)\n", r, names(v)[which.max(v)], max(v))) }
cat("\n=== 每个 module 的最高 celltype ===\n")
for (cc in colnames(mean_hME)) { v <- mean_hME[[cc]]; names(v) <- rownames(mean_hME); cat(sprintf("%-10s top=%-16s (%.3f)\n", cc, names(v)[which.max(v)], max(v))) }

df_6b <- melt(as.matrix(mean_hME))
colnames(df_6b) <- c("CellType", "Module", "MeanhME")

p6b <- ggplot(df_6b, aes(x = Module, y = CellType, fill = MeanhME)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027",
                       midpoint = 0, name = "Mean hME") +
  labs(x = "", y = "") +
  theme_minimal(base_size = 7) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 6.5, color = "black"),
    axis.text.y  = element_text(size = 6.5, color = "black"),
    axis.line    = element_line(colour = "black", linewidth = 0.3),
    axis.ticks   = element_line(colour = "black", linewidth = 0.3),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 5.5),
    legend.text  = element_text(size = 5),
    legend.key.size = unit(0.2, "cm"),
    plot.margin  = margin(2, 2, 2, 2, "mm")
  )

## 宽度 10cm = 与原版 figures/Fig6b_module_celltype_heatmap_v2.pdf(283x170pt) 一致,
## 拼版时才不会比原面板窄。
## [EN] Width 10 cm matches the original panel (283x170 pt) so it drops straight into the
## assembled figure without being narrower than the panel it replaces.
ggsave(file.path(outdir, "Fig6_panel_b_module_heatmap.pdf"),
       p6b, width = 10, height = 6, units = "cm", dpi = 600, device = cairo_pdf)
cat("(showtext)\n")
cat("FIG6B DONE ->", outdir, "\n")
