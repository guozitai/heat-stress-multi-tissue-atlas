# ============================================================
# Fig. 5 panels a-d — scRNA-seq (re-rendered)
#
# Part of the code for the Data Descriptor
#   "Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq data from
#    heat-stressed lactating Holstein cows"
#
# Produces: Fig5_panel_a_tSNE / _b_cell_proportion / _c_DEG_counts / _d_marker_dotplot
# showtext enabled to match the original panels; the unassigned cluster is excluded (82,013 -> 81,903 cells).
#
# Inputs are expected under DATA_DIR (default "data/"), or set HS_DATA_DIR.
# Archived at figshare https://doi.org/10.6084/m9.figshare.32993768
#
# The original working comments are in Chinese and are left unchanged;
# an English line has been added below each of them.
# ============================================================

# .libPaths(c("<site-specific R library path>", .libPaths()))
suppressMessages({library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
                  library(showtext); library(sysfonts)})

## ⚠ 必须开 showtext: 原版面板(2026-04-10)就是这么出的 —— 文字全轮廓化、0 可抽字符。
## 不开 showtext 时 cairo 自己排版, 会把每个字形的推进吸附到整数 pt:
## 5pt 字号下 Arial 的 r(1.67)->2 撑宽20%, t(1.39)->1 压窄28%, i 和 t 宽度还并成一样,
## 字母间距忽松忽紧。showtext 用 FreeType 精确排版, 绕开量化。
## [EN] showtext MUST be enabled: the original panels (2026-04-10) were produced this way — all text
## outlined, 0 extractable characters. Without showtext cairo lays out the text itself and quantises
## each glyph advance to whole points: at 5pt, Arial r(1.67)->2 is 20% too wide, t(1.39)->1 is 28% too
## narrow, and i and t collapse to the same width, giving visibly uneven letter spacing.
## showtext uses FreeType for exact layout and avoids the quantisation.
MSF <- "/usr/share/fonts/truetype/msttcorefonts/"
font_add("Arial", regular = paste0(MSF, "Arial.ttf"),
         bold = paste0(MSF, "Arial_Bold.ttf"),
         italic = paste0(MSF, "Arial_Italic.ttf"))
showtext_auto()
showtext_opts(dpi = 600)   # 必须与下面 ggsave(dpi=600) 一致, 否则字号会错
# [EN] dpi here MUST match ggsave(dpi=) below, otherwise font sizes come out wrong

DATA_DIR <- Sys.getenv("HS_DATA_DIR", "data")
outdir   <- Sys.getenv("HS_OUT_DIR", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rdsdir <- file.path(DATA_DIR, "scRNA")

obj <- readRDS(file.path(rdsdir, "integrated_25cl_annotated.rds"))
if (!"treatment" %in% colnames(obj@meta.data)) obj$treatment <- ifelse(grepl("^HS", obj$sample), "HS", "PFTN")
cat("loaded:", ncol(obj), "cells\n")

## ---- HSCs -> Unassigned ----
obj$celltype <- as.character(obj$celltype)
obj$celltype[obj$celltype == "HSCs"] <- "Unassigned"

## 定案(2026-07-17): 图里【删掉】Unassigned 这一行/一列 —— unassigned 的定义就是
## "没匹配到已知细胞类型", 图只展示匹配到的, 正文照旧说明它的存在与去向。
## 实测: 分母 82013 -> 81903 对各比例的影响 <=0.07 个百分点(T 47.74 vs 47.79), 可忽略。
## Technical Validation 里那句事实句(CD34/PROM1/MPO 阳性率全 0)保留 = 不是藏, 是归位。
## [EN] Decision (2026-07-17): the unassigned row/column is REMOVED from the figures — "unassigned"
## means "matched no known cell type", so the figures show only matched types; the text still states
## its existence and fate. Measured: denominator 82,013 -> 81,903 shifts each proportion by <=0.07
## percentage points (T 47.74 vs 47.79) = negligible. The Technical Validation sentence
## (CD34/PROM1/MPO all 0%) is kept — this is relocation, not concealment.
n_before <- ncol(obj)
obj <- obj[, obj$celltype != "Unassigned"]
cat("删 Unassigned:", n_before, "->", ncol(obj), "细胞 (应为 82013 -> 81903)
")

celltype_colors <- c(
  "T cells" = "#E64B35", "B cells" = "#4DBBD5", "Monocytes" = "#00A087",
  "NK cells" = "#3C5488", "Plasma cells" = "#F39B7F",
  "Dendritic cells" = "#8491B4"
)
obj$celltype <- factor(obj$celltype,
  levels = c("T cells","B cells","Monocytes","NK cells",
             "Plasma cells","Dendritic cells"))
cat("celltype AFTER relabel:\n"); print(table(obj$celltype))

theme_nature <- theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.text = element_text(size = 6, color = "black"),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 5.5),
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.3, "cm")
  )

## ---- panel a: t-SNE ----
pA <- DimPlot(obj, reduction = "tsne", group.by = "celltype",
              cols = celltype_colors, pt.size = 0.05) +
  theme_nature + ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 2)))
cat("panel a built\n")

## ---- panel b: cell proportion stacked bar ----
prop_df <- obj@meta.data %>%
  group_by(treatment, celltype) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(treatment) %>% mutate(pct = n / sum(n) * 100)
pB_prop <- ggplot(prop_df, aes(x = treatment, y = pct, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = celltype_colors) +
  labs(x = "", y = "Proportion (%)") + theme_nature + theme(legend.position = "right")
cat("panel b built\n")

## ---- panel c: DEG counts (缓存, FindMarkers 慢) ----
## [EN] ---- panel c: DEG counts (cached; FindMarkers is slow) ----
degfile <- file.path(outdir, "deg_counts.rds")
if (file.exists(degfile)) {
  deg_counts <- readRDS(degfile); cat("panel c: 读缓存 deg_counts.rds\n")
  deg_counts <- deg_counts[deg_counts$celltype != "Unassigned", ]   # 缓存是含 Unassigned 时算的
} else {
  Idents(obj) <- "celltype"
  deg_counts <- data.frame(celltype = character(), direction = character(),
                           count = integer(), stringsAsFactors = FALSE)
  for (ct in levels(obj$celltype)) {
    sub <- subset(obj, celltype == ct)
    if (length(unique(sub$treatment)) < 2) next
    Idents(sub) <- "treatment"
    tryCatch({
      degs <- FindMarkers(sub, ident.1 = "HS", ident.2 = "PFTN",
                          logfc.threshold = 0.25, min.pct = 0.1)
      sig <- degs[degs$p_val_adj < 0.05, ]
      n_up <- sum(sig$avg_log2FC > 0); n_down <- sum(sig$avg_log2FC < 0)
      deg_counts <<- rbind(deg_counts,
        data.frame(celltype = ct, direction = "Up", count = n_up),
        data.frame(celltype = ct, direction = "Down", count = -n_down))
    }, error = function(e) cat("Skipped:", ct, "\n"))
  }
  saveRDS(deg_counts, degfile); cat("panel c: FindMarkers 算完并缓存\n")
}
print(deg_counts)
deg_counts <- deg_counts[deg_counts$celltype != "Unassigned", ]
deg_counts$celltype <- factor(deg_counts$celltype,
  levels = rev(c("T cells","B cells","Monocytes","NK cells",
                 "Plasma cells","Dendritic cells")))
pC_deg <- ggplot(deg_counts, aes(x = celltype, y = count, fill = direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  coord_flip() + labs(x = "", y = "Number of DEGs") +
  geom_hline(yintercept = 0, linewidth = 0.3) + theme_nature
cat("panel c built\n")

## ---- panel d: marker dot plot ----
markers_fix <- c("CD3D","CD3E","CD4","CD8A","CD8B","CD79A","MS4A1","CD19",
  "CD14","S100A8","S100A9","NKG7","GNLY","KLRD1","JCHAIN","MZB1","FCER1A")
markers_fix <- markers_fix[markers_fix %in% rownames(obj[["RNA"]])]
pD_dot <- DotPlot(obj, features = markers_fix, group.by = "celltype",
              assay = "RNA", cols = c("lightgrey","#B2182B")) +
  RotatedAxis() + theme_nature +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1)) + ggtitle("")
cat("panel d built\n")

## ---- 只出四张子图, 按【面板实际位置】命名 ----
## [EN] ---- Export the four panels, named by their ACTUAL position in the assembled figure ----
D <- 600
ggsave(file.path(outdir, "Fig5_panel_a_tSNE.pdf"),            pA,      width = 10, height = 8, units = "cm", dpi = D, device = cairo_pdf)
ggsave(file.path(outdir, "Fig5_panel_b_cell_proportion.pdf"), pB_prop, width = 8,  height = 8, units = "cm", dpi = D, device = cairo_pdf)
ggsave(file.path(outdir, "Fig5_panel_c_DEG_counts.pdf"),      pC_deg,  width = 10, height = 6, units = "cm", dpi = D, device = cairo_pdf)
ggsave(file.path(outdir, "Fig5_panel_d_marker_dotplot.pdf"),  pD_dot,  width = 16, height = 6, units = "cm", dpi = D, device = cairo_pdf)
cat("FIG5 PANELS DONE (showtext)
")
