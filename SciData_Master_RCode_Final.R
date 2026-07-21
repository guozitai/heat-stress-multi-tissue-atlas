# ============================================================
# Master analysis and figure script — Data Descriptor
#
#   "Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq data from
#    heat-stressed lactating Holstein cows"
#   Zitai Guo, Shengtao Gao, Lu Ma, Dengpan Bu
#
# WHAT THIS FILE COVERS
#   Count matrix -> DESeq2 differential expression            (Fig. 2)
#   GO / KEGG pathway enrichment                              (Fig. 3)
#   miRNA differential expression and target networks         (Fig. 4)
#   scRNA-seq integration, SingleR annotation, cell types     (Fig. 5)
#   hdWGCNA co-expression modules and module GO enrichment    (Fig. 6a-c)
#   CellChat cell-cell communication                          (Fig. 6d)
#
# WHAT THIS FILE DOES NOT COVER
#   Bulk fastq -> alignment -> raw count matrix. Read alignment and counting
#   were performed by the sequencing provider; the resulting count matrix is
#   the starting point here and is archived at figshare (DOI below).
#   scRNA-seq upstream processing is in run_dnbc4tools_all.sh.
#
# DATA
#   Sequence data      NCBI GEO GSE338587
#   Processed data     figshare https://doi.org/10.6084/m9.figshare.32993768
#
# ENVIRONMENT
#   R 4.4.2 on Ubuntu 22.04. Full package versions in sessionInfo.txt.
#
# NOTE ON LANGUAGE
#   The original working comments are in Chinese and have been left unchanged.
#   An English line has been added below each of them.
# ============================================================

# ============================================================
# Scientific Data 手稿 — 全部Figure R代码归档 (Master File)
# ============================================================
# 项目：Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq
#       data from heat-stressed lactating Holstein cows
# 作者：Zitai Guo
# 归档日期：2026-04-10
#
# 每个Figure只保留最终版。
# Fig 1 为Gemini生成的实验设计图（PNG/PDF），无R代码。
#
# ============================================================
# 运行环境
# Runtime environment
# ============================================================
#   OS: Linux (服务器)
#   R version: 4.4.2
#   额外library路径（每次新session必须运行）：
#     .libPaths(c(.libPaths(),
#                 "<additional R library path 1>",
#                 "<additional R library path 2>"))
#
#   关键R包：
#     Seurat 5.1.0, hdWGCNA 0.4.11, CellChat 2.2.0,
#     DESeq2, ComplexHeatmap, circlize, clusterProfiler, org.Bt.eg.db,
#     SingleR, celldex, ggplot2, patchwork, cowplot, ggraph, igraph,
#     readxl, dplyr, tidyr, stringr, reshape2
#
#   字体：Arial / Helvetica（已安装 ttf-mscorefonts-installer）
#
# ============================================================
# 数据路径
# Input data paths (see DATA_DIR below)
# ============================================================
#   Bulk count matrix:  <DATA_DIR>/20240806-guo-5tissues.xlsx (sheet 3)
#   Bulk DEG:           <DATA_DIR>/DEG_results_all_tissues.rds
#   通路富集(up):       <DATA_DIR>/up_data_final.csv
#   通路富集(down):     <DATA_DIR>/down_data_final.csv
#   miRNA热图数据:      <DATA_DIR>/heatmap0303.xlsx
#   miRNA靶基因:        <DATA_DIR>/targetgene_enrich.xlsx
#   scRNA对象(原始):    <DATA_DIR>/scRNA/integrated_data_with_singleR.rds
#
# ============================================================
# 输出文件清单
# List of output files
# ============================================================
#   Fig 2a   Fig2a_PCA_final.pdf               五组织PCA
#   Fig 2b   Fig2b_correlation_v7.png           Spearman相关性热图
#   Fig 2c   Fig2d_jitter_v2.pdf                DEG jitter散点图
#   Fig 2图例 Fig2_legend_combined.pdf           Tissue+Treatment+Spearman r
#   Fig 3    Fig3_final_v7.pdf                  GO/KEGG气泡图(up/down)
#   Fig 4a   Fig4a_miRNA_heatmap_h_v3.pdf       miRNA热图横版
#   Fig 4b   Fig4b_network_v10.pdf              miRNA-靶基因网络图(圆形布局)
#   Fig 4图例 Fig4_legend_combined_v2.pdf        Treatment+Tissue+Z-score+Regulation
#   Fig 5a   Fig5a_tsne_celltype.pdf            t-SNE by cell type
#   Fig 5b   Fig5c_cellprop_bar.pdf             Cell proportion bar
#   Fig 5c   Fig5d_deg_counts.pdf               DEG counts per cell type
#   Fig 5d   Fig5b_dotplot_markers_v2.pdf       Marker gene dot plot
#   Fig 6a   Fig6a_dendrogram.pdf               hdWGCNA聚类树
#   Fig 6b   Fig6b_module_celltype_heatmap_v2.pdf Module-celltype热图(无grey)
#   Fig 6c   Fig6c_module_GO_v3.pdf             Module GO气泡图
#   Fig 6d   Fig6e_signaling_comparison.pdf     信息流比较(HS红/PFTN蓝)
#
# 中间RDS文件：
# Intermediate RDS objects:
#   integrated_25cl_annotated.rds   25 cluster + celltype + t-SNE
#   integrated_25cl_hdwgcna.rds     上面 + hdWGCNA模块
#
# ============================================================
# 全局作图标准（Nature系）
# Global figure style (Nature-journal conventions)
# ============================================================
#   字体：Helvetica / Arial
#   轴线宽：0.3pt   轴文字：6pt   轴标题：7pt   标题：8pt bold
#   图例标题：5.5pt  图例文字：5pt
#   Panel标签：8pt bold 小写(a,b,c) — 在Illustrator中添加
#
# 配色方案：
# Colour scheme:
#   组织：Muscle=#66C2A5  Adipose=#FC8D62  Rumen=#8DA0CB
#         Liver=#E78AC3   Mammary=#A6D854
#   Treatment：HS=#D73027  PFTN=#4575B4
#   细胞类型：T cells=#E64B35  B cells=#4DBBD5  Monocytes=#00A087
#             NK cells=#3C5488  Plasma cells=#F39B7F
#             Dendritic cells=#8491B4  HSCs=#91D1C2
# ============================================================


# ============================================================
# 公共设置
# Common setup
# ============================================================

# ---- Input/output locations -------------------------------------------
# Place the input files listed in the header under DATA_DIR (default "data/"),
# or point HS_DATA_DIR at them:   export HS_DATA_DIR=/path/to/data
# The inputs are archived at figshare https://doi.org/10.6084/m9.figshare.32993768
DATA_DIR <- Sys.getenv("HS_DATA_DIR", "data")
# -----------------------------------------------------------------------

# Site-specific R library paths from the original analysis environment.
# Not required to reproduce the figures; edit or delete as appropriate.
# .libPaths(c(.libPaths(), "<additional R library path>"))

outdir <- Sys.getenv("HS_OUT_DIR", "figures")
rdsdir <- file.path(DATA_DIR, "scRNA")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ############################################################
# Fig 2a — 五组织PCA (纵向 6×12 inch)
# Fig 2a — PCA of five tissues (portrait, 6x12 inch)
# 输出: Fig2a_PCA_final.pdf
# Output: Fig2a_PCA_final.pdf
# ############################################################

library(ggplot2)
library(readxl)
library(DESeq2)
library(dplyr)

counts_raw <- read_excel(file.path(DATA_DIR, "20240806-guo-5tissues.xlsx"), sheet = 3)
gene_ids <- counts_raw[[1]]
keep <- !is.na(gene_ids) & !duplicated(gene_ids)
gene_ids <- gene_ids[keep]
counts_mat <- as.data.frame(lapply(counts_raw[keep, -1], function(x) as.numeric(as.character(x))))
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- colnames(counts_raw)[-1]
counts_mat[is.na(counts_mat)] <- 0
counts_mat <- round(as.matrix(counts_mat))

sample_names <- colnames(counts_mat)
meta_df <- data.frame(
  sample_id = sample_names,
  tissue = case_when(
    grepl("muscle|Muscle", sample_names) ~ "Muscle",
    grepl("fat", sample_names)           ~ "Adipose",
    grepl("rumen|Rumen", sample_names)   ~ "Rumen",
    grepl("LIVER", sample_names)         ~ "Liver",
    grepl("Mammary", sample_names)       ~ "Mammary"
  ),
  treatment = ifelse(grepl("^HS", sample_names), "HS", "PFTN"),
  row.names = sample_names
)

dds <- DESeqDataSetFromMatrix(counts_mat, colData = meta_df, design = ~ 1)
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("tissue", "treatment"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$Tissue <- factor(pca_data$tissue,
  levels = c("Muscle","Adipose","Rumen","Liver","Mammary"))
pca_data$Treatment <- factor(pca_data$treatment, levels = c("HS","PFTN"))

tissue_cols <- c(
  Muscle  = "#66C2A5", Adipose = "#FC8D62", Rumen = "#8DA0CB",
  Liver   = "#E78AC3", Mammary = "#A6D854"
)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Tissue, shape = Treatment)) +
  geom_point(size = 1.8) +
  scale_color_manual(values = tissue_cols) +
  scale_shape_manual(values = c(HS = 17, PFTN = 16)) +
  labs(x = paste0("PC1 (", percentVar[1], "%)"),
       y = paste0("PC2 (", percentVar[2], "%)")) +
  theme_classic() +
  theme(
    text            = element_text(family = "Arial"),
    axis.text       = element_text(size = 9, color = "black"),
    axis.title      = element_text(size = 11),
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(5, "mm"),
    axis.line       = element_line(linewidth = 0.5),
    axis.ticks      = element_line(linewidth = 0.5),
    legend.position = "right"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2, override.aes = list(
      size = 3, color = c("#D73027", "#4575B4")
    ))
  )

pdf(file.path(outdir, "Fig2a_PCA_final.pdf"), width = 6, height = 12)
print(p)
dev.off()
cat("Fig 2a done!\n")


# ############################################################
# Fig 2b — Spearman相关性热图 (棕黄色系)
# Fig 2b — Spearman correlation heatmap (brown-yellow palette)
# 输出: Fig2b_correlation_v7.png (6×6cm, 300dpi)
# Output: Fig2b_correlation_v7.png (6x6 cm, 300 dpi)
#        Fig2_panelB.pdf (8×8cm, 带legend版备用)
#         Fig2_panelB.pdf (8x8 cm, alternative version carrying the legend)
# ############################################################

library(ComplexHeatmap)
library(circlize)

# --- counts_mat 和 meta_df 已在 Fig 2a 中创建 ---
# --- counts_mat and meta_df were created in the Fig 2a section above ---

meta_df$tissue_display <- case_when(
  grepl("muscle|Muscle", meta_df$tissue, ignore.case = TRUE) ~ "Muscle",
  grepl("fat|Adipose", meta_df$tissue, ignore.case = TRUE)   ~ "Adipose",
  grepl("rumen|Rumen", meta_df$tissue, ignore.case = TRUE)   ~ "Rumen",
  grepl("liver|LIVER|Liver", meta_df$tissue, ignore.case = TRUE) ~ "Liver",
  grepl("mammary|Mammary", meta_df$tissue, ignore.case = TRUE)   ~ "Mammary gland",
  TRUE ~ meta_df$tissue
)

vsd_mat <- assay(vsd)
cor_mat <- cor(vsd_mat, method = "spearman")

tissue_order_new <- c("Muscle", "Adipose", "Rumen", "Liver", "Mammary gland")
sample_order_new <- rownames(meta_df)[order(
  match(meta_df$tissue_display, tissue_order_new), meta_df$treatment)]
cor_mat_new <- cor_mat[sample_order_new, sample_order_new]
meta_new <- meta_df[sample_order_new, ]

tissue_colors_hm <- c("Muscle" = "#66C2A5", "Adipose" = "#FC8D62",
                       "Rumen" = "#8DA0CB", "Liver" = "#E78AC3",
                       "Mammary gland" = "#A6D854")
treat_colors <- c("HS" = "#D73027", "PFTN" = "#4575B4")

col_fun_cor <- colorRamp2(c(0.6, 0.75, 0.9, 1), c("#FFF8E7", "#FEC44F", "#D95F0E", "#662506"))

ht_nolegend <- Heatmap(cor_mat_new,
  col = col_fun_cor, name = "r",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(
    Tissue = meta_new$tissue_display, Treatment = meta_new$treatment,
    col = list(Tissue = tissue_colors_hm, Treatment = treat_colors),
    show_annotation_name = FALSE, simple_anno_size = unit(2, "mm"),
    show_legend = FALSE
  ),
  left_annotation = rowAnnotation(
    Tissue = meta_new$tissue_display, Treatment = meta_new$treatment,
    col = list(Tissue = tissue_colors_hm, Treatment = treat_colors),
    show_annotation_name = FALSE, simple_anno_size = unit(2, "mm"),
    show_legend = FALSE
  ),
  show_heatmap_legend = FALSE,
  border = TRUE, use_raster = TRUE, raster_quality = 3
)

png(file.path(outdir, "Fig2b_correlation_v7.png"),
    width = 6, height = 6, units = "cm", res = 300)
draw(ht_nolegend, padding = unit(c(1, 1, 1, 1), "mm"))
dev.off()

pdf(file.path(outdir, "Fig2_panelB.pdf"), width = 8/2.54, height = 8/2.54)
draw(ht_nolegend, padding = unit(c(1, 1, 1, 1), "mm"))
dev.off()

cat("Fig 2b done!\n")


# ############################################################
# Fig 2c — 纵向jitter散点图 (DEG)
# Fig 2c — vertical jitter plot of DEG counts
# 输出: Fig2d_jitter_v2.pdf (8×8cm)
# Output: Fig2d_jitter_v2.pdf (8x8 cm)
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
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6, face = "bold", color = "black"),
    axis.title.x = element_text(size = 7),
    axis.text.x = element_text(size = 6, color = "black"),
    axis.line.y = element_line(color = "black", linewidth = 0.3),
    axis.ticks.y = element_line(color = "black", linewidth = 0.3),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

ggsave(file.path(outdir, "Fig2d_jitter_v2.pdf"), p_d,
       width = 8, height = 8, units = "cm", dpi = 600, device = cairo_pdf)

cat("Fig 2c done!\n")


# ############################################################
# Fig 2 组合图例 (Tissue + Treatment + Spearman r)
# Fig 2 combined legend (Tissue + Treatment + Spearman r)
# 输出: Fig2_legend_combined.pdf
# Output: Fig2_legend_combined.pdf
# ############################################################

library(cowplot)

p_pca_dummy <- ggplot(data.frame(
    x = 1:7, y = 1:7,
    Tissue = factor(c("Muscle","Adipose","Rumen","Liver","Mammary","Muscle","Adipose"),
                    levels = c("Muscle","Adipose","Rumen","Liver","Mammary")),
    Treatment = factor(c("HS","HS","HS","HS","HS","PFTN","PFTN"),
                       levels = c("HS","PFTN"))
  ), aes(x, y, color = Tissue, shape = Treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = tissue_cols) +
  scale_shape_manual(values = c(HS = 17, PFTN = 16)) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2, override.aes = list(
      size = 3, color = c("#D73027", "#4575B4")
    ))
  ) +
  theme_void() +
  theme(
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(5, "mm")
  )
leg_pca <- get_legend(p_pca_dummy)

p_sp_dummy <- ggplot(data.frame(x = 1, y = 1, r = 0.8), aes(x, y, fill = r)) +
  geom_tile() +
  scale_fill_gradient(low = "#FFF7EC", high = "#7F2704", name = "Spearman r",
                      limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  guides(fill = guide_colorbar(
    barwidth = unit(5, "mm"), barheight = unit(30, "mm"),
    title.position = "top", ticks = TRUE
  )) +
  theme_void() +
  theme(
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 9)
  )
leg_sp <- get_legend(p_sp_dummy)

combined_legend <- plot_grid(leg_pca, leg_sp, ncol = 1,
                             rel_heights = c(1, 0.6), align = "v")

pdf(file.path(outdir, "Fig2_legend_combined.pdf"), width = 2, height = 5)
print(combined_legend)
dev.off()
cat("Fig 2 legend done!\n")


# ############################################################
# Fig 3 — 通路富集气泡图 (Nature红蓝配色, V7最终版)
# Fig 3 — pathway enrichment bubble plot (red/blue palette, final v7)
# 输出: Fig3_final_v7.pdf (18×10cm)
# Output: Fig3_final_v7.pdf (18x10 cm)
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
ggsave(file.path(outdir, "Fig3_final_v7.pdf"), fig3,
       width = 18, height = 10, units = "cm", dpi = 600, device = cairo_pdf)

cat("Fig 3 done!\n")


# ############################################################
# Fig 4a — miRNA差异表达热图 (横版 v3)
# Fig 4a — differentially expressed miRNA heatmap (landscape, v3)
# 输出: Fig4a_miRNA_heatmap_h_v3.pdf (17×6cm)
# Output: Fig4a_miRNA_heatmap_h_v3.pdf (17x6 cm)
# ############################################################

hm_data <- read_excel(file.path(DATA_DIR, "heatmap0303.xlsx"), sheet = 1)
hm_data <- as.data.frame(hm_data)

count_cols <- grep("^HS|^PF", colnames(hm_data))
expr_mat <- as.matrix(hm_data[, count_cols])
rownames(expr_mat) <- hm_data$miRNA_name

log2_mat <- log2(expr_mat + 1)
scaled_mat <- t(scale(t(log2_mat)))
scaled_t <- t(scaled_mat)

tissue_vec <- ifelse(hm_data$tissue == "muscle", "Muscle", "Adipose")
names(tissue_vec) <- hm_data$miRNA_name
treat_vec <- rep(c("HS","PFTN"), each = 4)
names(treat_vec) <- rownames(scaled_t)

tissue_colors_mirna <- c("Muscle" = "#66C2A5", "Adipose" = "#FC8D62")
col_fun_hm <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

ha_top_mirna <- HeatmapAnnotation(
  Tissue = tissue_vec,
  col = list(Tissue = tissue_colors_mirna),
  show_annotation_name = FALSE,
  simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(
    Tissue = list(title_gp = gpar(fontsize = 5.5, fontfamily = "Helvetica"),
                  labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"))
  )
)

ha_left_mirna <- rowAnnotation(
  Treatment = treat_vec,
  col = list(Treatment = treat_colors),
  show_annotation_name = FALSE,
  simple_anno_size = unit(2.5, "mm"),
  annotation_legend_param = list(
    Treatment = list(title_gp = gpar(fontsize = 5.5, fontfamily = "Helvetica"),
                     labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"))
  )
)

do_heatmap_4a <- function() {
  ht <- Heatmap(scaled_t,
    col = col_fun_hm, name = "Z-score",
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 5.5, fontfamily = "Helvetica"),
    top_annotation = ha_top_mirna,
    left_annotation = ha_left_mirna,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 5.5, fontfamily = "Helvetica"),
      labels_gp = gpar(fontsize = 5, fontfamily = "Helvetica"),
      legend_height = unit(1.5, "cm")
    ),
    column_dend_height = unit(6, "mm"),
    row_split = treat_vec, row_gap = unit(1.5, "mm"),
    row_title_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
    row_title_rot = 0,
    border = TRUE, border_gp = gpar(col = "grey40", lwd = 0.5)
  )
  draw(ht, background = "white", padding = unit(c(3,3,3,3), "mm"),
       merge_legends = FALSE)
}

pdf(file.path(outdir, "Fig4a_miRNA_heatmap_h_v3.pdf"),
    width = 17/2.54, height = 6/2.54)
do_heatmap_4a()
dev.off()

cat("Fig 4a done!\n")


# ############################################################
# Fig 4b — miRNA-靶基因网络图 (圆形布局 v10, 统一配色)
# Fig 4b — miRNA-target gene network (circular layout v10, unified colours)
# 输出: Fig4b_network_v10.pdf (17×17cm正方形)
# Output: Fig4b_network_v10.pdf (17x17 cm, square)
# ############################################################

filter <- dplyr::filter

dat <- read_excel(file.path(DATA_DIR, "targetgene_enrich.xlsx")) %>%
  setNames(c("mirna","gene","regulation","tissue")) %>%
  mutate(
    tissue = case_when(
      str_detect(tissue, regex("fat|adipose", TRUE)) ~ "Adipose",
      str_detect(tissue, regex("muscle", TRUE))      ~ "Muscle"),
    regulation = case_when(
      str_detect(regulation, regex("up", TRUE))   ~ "Up",
      str_detect(regulation, regex("down", TRUE)) ~ "Down")
  ) %>% filter(!is.na(tissue))

TOP_N <- 4

make_polar <- function(data, cx, cy, r_mir = 1.6, r_gene = 4.2) {
  dat_t  <- data %>% group_by(mirna) %>% slice_head(n = TOP_N) %>% ungroup()
  mirnas <- unique(dat_t$mirna)
  n_mir  <- length(mirnas)
  mir_ang <- seq(pi/2, pi/2 + 2*pi, length.out = n_mir + 1)[1:n_mir]
  mir_nodes <- list(); gene_nodes <- list(); edges <- list()
  for (i in seq_along(mirnas)) {
    m <- mirnas[i]; ang <- mir_ang[i]
    mx <- cx + r_mir * cos(ang); my <- cy + r_mir * sin(ang)
    genes_i <- dat_t %>% filter(mirna == m)
    n_g <- nrow(genes_i)
    sw <- 2*pi / n_mir * 0.72
    g_angs <- if(n_g == 1) ang else seq(ang - sw/2, ang + sw/2, length.out = n_g)
    mir_nodes[[i]] <- data.frame(name = m, x = mx, y = my,
      size = 12 + n_g * 2, tissue = unique(data$tissue), stringsAsFactors = FALSE)
    for (j in seq_len(n_g)) {
      gx <- cx + r_gene * cos(g_angs[j]); gy <- cy + r_gene * sin(g_angs[j])
      gene_nodes[[length(gene_nodes) + 1]] <- data.frame(
        name = genes_i$gene[j], x = gx, y = gy,
        regulation = genes_i$regulation[j], stringsAsFactors = FALSE)
      edges[[length(edges) + 1]] <- data.frame(
        x0 = mx, y0 = my, x1 = gx, y1 = gy,
        regulation = genes_i$regulation[j], stringsAsFactors = FALSE)
    }
  }
  list(mir = bind_rows(mir_nodes), gene = bind_rows(gene_nodes), edge = bind_rows(edges))
}

adi <- make_polar(filter(dat, tissue == "Adipose"), cx = -4.5, cy = 0)
mus <- make_polar(filter(dat, tissue == "Muscle"),  cx =  4.5, cy = 0)

# 统一配色版（最终版）
col_adi_mir  <- "#FC8D62"
col_adi_gene <- "#FDE0D0"
col_adi_brd  <- "#E07040"
col_mus_mir  <- "#66C2A5"
col_mus_gene <- "#D4EDE4"
col_mus_brd  <- "#4AA88A"
col_up       <- "#C0392B"
col_dn       <- "#2471A3"

draw_ring <- function(cx, cy, r, col, lty = 3, lwd = 0.4) {
  theta <- seq(0, 2*pi, length.out = 300)
  lines(cx + r*cos(theta), cy + r*sin(theta), col = col, lty = lty, lwd = lwd)
}
draw_edges <- function(edges) {
  for (i in seq_len(nrow(edges))) {
    e <- edges[i, ]
    col <- if(e$regulation == "Up") col_up else col_dn
    arrows(e$x0, e$y0, e$x1, e$y1, length = 0.055, angle = 18, code = 2,
           col = adjustcolor(col, 0.7), lwd = 0.75)
  }
}
draw_genes <- function(genes, fill, border) {
  points(genes$x, genes$y, pch = 21, cex = 1.0, bg = fill, col = border, lwd = 0.7)
}
label_genes <- function(genes, cx, cy, cex = 0.36) {
  genes_u <- genes %>%
    mutate(dist = sqrt((x - cx)^2 + (y - cy)^2)) %>%
    group_by(name) %>% slice_max(dist, n = 1, with_ties = FALSE) %>% ungroup()
  for (i in seq_len(nrow(genes_u))) {
    ang <- atan2(genes_u$y[i] - cy, genes_u$x[i] - cx)
    ox <- cos(ang) * 0.42; oy <- sin(ang) * 0.42
    adj_x <- if(cos(ang) > 0.15) 0 else if(cos(ang) < -0.15) 1 else 0.5
    text(genes_u$x[i] + ox, genes_u$y[i] + oy, genes_u$name[i],
         cex = cex, col = "grey20", adj = c(adj_x, 0.5))
  }
}
draw_mirs <- function(mirs, fill, border) {
  points(mirs$x, mirs$y, pch = 21, cex = mirs$size/11,
         bg = fill, col = border, lwd = 1.6)
}
label_mirs <- function(mirs, cx, cy, cex = 0.46) {
  lbl <- sub("bta-", "", mirs$name)
  for (i in seq_len(nrow(mirs))) {
    ang <- atan2(mirs$y[i] - cy, mirs$x[i] - cx)
    tx <- mirs$x[i] - cos(ang) * 0.08; ty <- mirs$y[i] - sin(ang) * 0.08
    text(tx, ty, lbl[i], cex = cex, col = "black", font = 2, adj = c(0.5, 0.5))
  }
}

do_plot_4b <- function() {
  par(mar = c(1.5, 1, 2.5, 1), bg = "white", family = "Arial")
  plot(NULL, xlim = c(-9.5, 9.5), ylim = c(-5.8, 5.8),
       asp = 1, axes = FALSE, xlab = "", ylab = "")
  draw_ring(-4.5, 0, 1.6, adjustcolor(col_adi_mir, 0.25))
  draw_ring(-4.5, 0, 4.2, adjustcolor(col_adi_mir, 0.12))
  draw_ring( 4.5, 0, 1.6, adjustcolor(col_mus_mir, 0.25))
  draw_ring( 4.5, 0, 4.2, adjustcolor(col_mus_mir, 0.12))
  abline(v = 0, col = "grey82", lty = 2, lwd = 0.7)
  draw_edges(adi$edge); draw_edges(mus$edge)
  draw_genes(adi$gene, col_adi_gene, col_adi_brd)
  draw_genes(mus$gene, col_mus_gene, col_mus_brd)
  label_genes(adi$gene, -4.5, 0); label_genes(mus$gene, 4.5, 0)
  draw_mirs(adi$mir, col_adi_mir, "white"); draw_mirs(mus$mir, col_mus_mir, "white")
  label_mirs(adi$mir, -4.5, 0); label_mirs(mus$mir, 4.5, 0)
  text(-4.5, 5.5, "Adipose", cex = 1.0, font = 2, col = col_adi_brd, adj = c(0.5, 0.5))
  text( 4.5, 5.5, "Muscle",  cex = 1.0, font = 2, col = col_mus_mir, adj = c(0.5, 0.5))
  legend("bottomright", inset = c(0.005, 0.005),
         legend = c("Up-regulated","Down-regulated"),
         col = c(col_up, col_dn), lty = 1, lwd = 1.5, pch = NA,
         cex = 0.48, bty = "n", title = "Regulation",
         title.col = "grey20", title.cex = 0.52)
  legend("bottomleft", inset = c(0.005, 0.005),
         legend = c("Adipose miRNA","Muscle miRNA"),
         pt.bg = c(col_adi_mir, col_mus_mir),
         pch = 21, pt.cex = 1.1, col = "white", pt.lwd = 0.8,
         cex = 0.48, bty = "n", title = "Tissue",
         title.col = "grey20", title.cex = 0.52)
}

cairo_pdf(file.path(outdir, "Fig4b_network_v10.pdf"), width = 6.7, height = 6.7)
do_plot_4b()
dev.off()

cat("Fig 4b done!\n")


# ############################################################
# Fig 4 组合图例 (Treatment + Tissue + Z-score + Regulation)
# Fig 4 combined legend (Treatment + Tissue + Z-score + Regulation)
# 输出: Fig4_legend_combined_v2.pdf
# Output: Fig4_legend_combined_v2.pdf
# ############################################################

common_theme <- theme_void() +
  theme(
    legend.title      = element_text(size = 12, face = "bold"),
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(8, "mm"),
    legend.key.height = unit(6, "mm"),
    legend.spacing.y  = unit(1, "mm")
  )

leg_trt <- get_legend(
  ggplot(data.frame(x = 1:2, y = 1:2,
    g = factor(c("HS","PFTN"), levels = c("HS","PFTN"))),
    aes(x, y, fill = g)) +
    geom_tile() +
    scale_fill_manual(values = c(HS = "#D73027", PFTN = "#4575B4"), name = "Treatment") +
    common_theme
)

leg_tis <- get_legend(
  ggplot(data.frame(x = 1:2, y = 1:2,
    g = factor(c("Adipose","Muscle"), levels = c("Adipose","Muscle"))),
    aes(x, y, fill = g)) +
    geom_tile() +
    scale_fill_manual(values = c(Adipose = "#FC8D62", Muscle = "#66C2A5"), name = "Tissue") +
    common_theme
)

leg_zs <- get_legend(
  ggplot(data.frame(x = 1, y = 1, z = 0), aes(x, y, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, name = "Z-score", limits = c(-2, 2), breaks = c(-2,-1,0,1,2)) +
    guides(fill = guide_colorbar(barwidth = unit(8, "mm"), barheight = unit(30, "mm"),
      title.position = "top", ticks = TRUE)) +
    theme_void() +
    theme(legend.title = element_text(size = 12, face = "bold"),
          legend.text  = element_text(size = 11))
)

leg_reg <- get_legend(
  ggplot(data.frame(x = c(1,1), xend = c(2,2), y = c(2,1), yend = c(2,1),
    g = factor(c("Up-regulated","Down-regulated"),
               levels = c("Up-regulated","Down-regulated"))),
    aes(x = x, xend = xend, y = y, yend = yend, color = g)) +
    geom_segment(linewidth = 1.5) +
    scale_color_manual(values = c("Up-regulated" = "#C0392B", "Down-regulated" = "#2471A3"),
                       name = "Regulation") +
    common_theme + theme(legend.key.width = unit(8, "mm"))
)

combined <- plot_grid(leg_trt, leg_tis, leg_zs, leg_reg,
                      ncol = 1, rel_heights = c(1, 0.8, 1.4, 0.8), align = "v")

pdf(file.path(outdir, "Fig4_legend_combined_v2.pdf"), width = 2, height = 5.5)
print(combined)
dev.off()
cat("Fig 4 legend done!\n")


# ############################################################
# SingleR 细胞类型注释 + 25 cluster 重聚类 + 7种细胞类型映射
# SingleR cell-type annotation + re-clustering to 25 clusters + mapping to 7 cell types
# 输出: integrated_25cl_annotated.rds
# Output: integrated_25cl_annotated.rds
# ############################################################

library(Seurat)
library(SingleR)
library(celldex)

obj <- readRDS(file.path(rdsdir, "integrated_data_with_singleR.rds"))

# 重聚类到25 clusters
# Re-cluster into 25 clusters
obj <- FindClusters(obj, resolution = 0.85)

# SingleR (cluster-level, MonacoImmuneData)
ref <- celldex::MonacoImmuneData()
expr <- GetAssayData(obj, assay = "RNA", slot = "data")
clusters <- obj$seurat_clusters

pred_main <- SingleR(test = expr, ref = ref, labels = ref$label.main,
                     clusters = clusters, assay.type.test = "logcounts")
pred_fine <- SingleR(test = expr, ref = ref, labels = ref$label.fine,
                     clusters = clusters, assay.type.test = "logcounts")

# 手动映射到7种细胞类型（基于博士论文+marker验证）
# Manual mapping to 7 cell types (prior annotation + marker validation)
cluster_to_celltype <- c(
  "0"  = "T cells",      "1"  = "T cells",      "2"  = "B cells",
  "3"  = "B cells",      "4"  = "B cells",      "5"  = "Monocytes",
  "6"  = "T cells",      "7"  = "Monocytes",    "8"  = "T cells",
  "9"  = "Monocytes",    "10" = "Monocytes",    "11" = "T cells",
  "12" = "NK cells",     "13" = "Plasma cells",  "14" = "B cells",
  "15" = "Monocytes",    "16" = "B cells",      "17" = "T cells",
  "18" = "Dendritic cells", "19" = "T cells",   "20" = "B cells",
  "21" = "T cells",      "22" = "T cells",      "23" = "T cells",
  "24" = "HSCs"
)

celltype <- cluster_to_celltype[as.character(obj$seurat_clusters)]
names(celltype) <- Cells(obj)
obj <- AddMetaData(obj, metadata = celltype, col.name = "celltype")
obj$treatment <- ifelse(grepl("^HS", obj$sample), "HS", "PFTN")
obj <- RunTSNE(obj, dims = 1:20)

saveRDS(obj, file.path(rdsdir, "integrated_25cl_annotated.rds"))
cat("SingleR + 25cl annotation done!\n")


# ############################################################
# Fig 5a — t-SNE by cell type
# 输出: Fig5a_tsne_celltype.pdf (10×8cm)
# Output: Fig5a_tsne_celltype.pdf (10x8 cm)
# ############################################################

library(tidyr)

celltype_colors <- c(
  "T cells" = "#E64B35", "B cells" = "#4DBBD5", "Monocytes" = "#00A087",
  "NK cells" = "#3C5488", "Plasma cells" = "#F39B7F",
  "Dendritic cells" = "#8491B4", "HSCs" = "#91D1C2"
)

obj$celltype <- factor(obj$celltype,
  levels = c("T cells","B cells","Monocytes","NK cells",
             "Plasma cells","Dendritic cells","HSCs"))

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

pA <- DimPlot(obj, reduction = "tsne", group.by = "celltype",
              cols = celltype_colors, pt.size = 0.05) +
  theme_nature + ggtitle("") +
  guides(color = guide_legend(override.aes = list(size = 2)))

ggsave(file.path(outdir, "Fig5a_tsne_celltype.pdf"),
       pA, width = 10, height = 8, units = "cm")

cat("Fig 5a done!\n")


# ############################################################
# Fig 5b — Cell proportion stacked bar (HS vs PFTN)
# 输出: Fig5c_cellprop_bar.pdf (8×8cm)
# Output: Fig5c_cellprop_bar.pdf (8x8 cm)
# ############################################################

prop_df <- obj@meta.data %>%
  group_by(treatment, celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(treatment) %>%
  mutate(pct = n / sum(n) * 100)

pC <- ggplot(prop_df, aes(x = treatment, y = pct, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = celltype_colors) +
  labs(x = "", y = "Proportion (%)") +
  theme_nature + theme(legend.position = "right")

ggsave(file.path(outdir, "Fig5c_cellprop_bar.pdf"),
       pC, width = 8, height = 8, units = "cm")

cat("Fig 5b done!\n")


# ############################################################
# Fig 5c — DEG counts per cell type
# 输出: Fig5d_deg_counts.pdf (10×6cm)
# Output: Fig5d_deg_counts.pdf (10x6 cm)
# ############################################################

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
    n_up <- sum(sig$avg_log2FC > 0)
    n_down <- sum(sig$avg_log2FC < 0)
    deg_counts <- rbind(deg_counts,
      data.frame(celltype = ct, direction = "Up", count = n_up),
      data.frame(celltype = ct, direction = "Down", count = -n_down))
  }, error = function(e) cat("Skipped:", ct, "\n"))
}

deg_counts$celltype <- factor(deg_counts$celltype, levels = rev(levels(obj$celltype)))

pD <- ggplot(deg_counts, aes(x = celltype, y = count, fill = direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Up" = "#D73027", "Down" = "#4575B4")) +
  coord_flip() + labs(x = "", y = "Number of DEGs") +
  geom_hline(yintercept = 0, linewidth = 0.3) + theme_nature

ggsave(file.path(outdir, "Fig5d_deg_counts.pdf"),
       pD, width = 10, height = 6, units = "cm")

cat("Fig 5c done!\n")


# ############################################################
# Fig 5d — Marker gene dot plot
# 输出: Fig5b_dotplot_markers_v2.pdf (16×6cm)
# Output: Fig5b_dotplot_markers_v2.pdf (16x6 cm)
# ############################################################

markers_fix <- c(
  "CD3D","CD3E","CD4","CD8A","CD8B",
  "CD79A","MS4A1","CD19",
  "CD14","S100A8","S100A9",
  "NKG7","GNLY","KLRD1",
  "JCHAIN","MZB1",
  "FCER1A"
)
markers_fix <- markers_fix[markers_fix %in% rownames(obj[["RNA"]])]

pB <- DotPlot(obj, features = markers_fix, group.by = "celltype",
              assay = "RNA", cols = c("lightgrey","#B2182B")) +
  RotatedAxis() + theme_nature +
  theme(axis.text.x = element_text(size = 5, angle = 45, hjust = 1)) +
  ggtitle("")

ggsave(file.path(outdir, "Fig5b_dotplot_markers_v2.pdf"),
       pB, width = 16, height = 6, units = "cm")

cat("Fig 5d done!\n")


# ############################################################
# Fig 6a — hdWGCNA dendrogram
# Fig 6b — Module-celltype heatmap (去grey, v2)
# Fig 6b — module-celltype heatmap (grey module removed, v2)
# Fig 6c — Module GO enrichment (Fig3风格, v3)
# Fig 6c — module GO enrichment (Fig 3 style, v3)
# ############################################################

library(hdWGCNA)
library(reshape2)
library(clusterProfiler)
library(org.Bt.eg.db)

# --- 重建v5格式Assay (hdWGCNA要求) ---
# --- Rebuild v5-format Assay (required by hdWGCNA) ---
if (inherits(obj[["RNA"]], "Assay")) {
  counts_mat_sc <- GetAssayData(obj, assay = "RNA", slot = "counts")
  data_mat_sc   <- GetAssayData(obj, assay = "RNA", slot = "data")
  obj2 <- CreateSeuratObject(counts = counts_mat_sc, meta.data = obj@meta.data)
  obj2[["RNA"]]$data <- data_mat_sc
  obj2@reductions <- obj@reductions
  obj2@graphs     <- obj@graphs
  obj2$treatment <- ifelse(grepl("^HS", obj2$sample), "HS", "PFTN")
} else {
  obj2 <- obj
}

# hdWGCNA pipeline
obj2 <- SetupForWGCNA(obj2, gene_select = "fraction", fraction = 0.05, wgcna_name = "PBMC")
obj2 <- MetacellsByGroups(obj2, group.by = c("celltype","treatment"),
  k = 25, ident.group = "celltype", max_shared = 10)
obj2 <- NormalizeMetacells(obj2)
obj2 <- SetDatExpr(obj2,
  group_name = c("T cells","B cells","Monocytes","NK cells","Plasma cells","Dendritic cells"),
  group.by = "celltype", assay = "RNA")
obj2 <- TestSoftPowers(obj2, networkType = "signed")
obj2 <- ConstructNetwork(obj2, soft_power = 7, tom_name = "PBMC_TOM",
  networkType = "signed", minModuleSize = 30, mergeCutHeight = 0.25)

# Fig 6a: dendrogram
pdf(file.path(outdir, "Fig6a_dendrogram.pdf"), width = 6, height = 3)
PlotDendrogram(obj2, main = "")
dev.off()
cat("Fig 6a done!\n")

# ModuleEigengenes
obj2 <- ScaleData(obj2, features = rownames(obj2))
obj2 <- ModuleEigengenes(obj2)
obj2 <- ModuleConnectivity(obj2)

# Fig 6b: module-celltype heatmap (去grey)
# Fig 6b: module-celltype heatmap (grey module removed)
hMEs <- GetMEs(obj2, harmonized = TRUE)
meta_sc <- obj2@meta.data
hMEs$celltype <- meta_sc$celltype
hMEs <- hMEs[, !grepl("grey", colnames(hMEs)) | colnames(hMEs) == "celltype"]

mean_hME <- aggregate(. ~ celltype, data = hMEs, FUN = mean)
rownames(mean_hME) <- mean_hME$celltype
mean_hME$celltype <- NULL
colnames(mean_hME) <- gsub("^ME", "", colnames(mean_hME))

mod_order_6b <- c("black","blue","brown","green","red","turquoise","yellow")
mod_order_6b <- mod_order_6b[mod_order_6b %in% colnames(mean_hME)]
mean_hME <- mean_hME[, mod_order_6b]

ct_order <- c("T cells","B cells","Monocytes","NK cells",
              "Plasma cells","Dendritic cells","HSCs")
ct_order <- ct_order[ct_order %in% rownames(mean_hME)]
mean_hME <- mean_hME[ct_order, ]

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

ggsave(file.path(outdir, "Fig6b_module_celltype_heatmap_v2.pdf"),
       p6b, width = 8, height = 6, units = "cm", dpi = 600, device = cairo_pdf)
cat("Fig 6b done!\n")

# Fig 6c: Module GO enrichment (Fig3风格 v3)
# Fig 6c: module GO enrichment (Fig 3 style, v3)
modules_6c <- GetModules(obj2)
modules_6c <- modules_6c[modules_6c$module != "grey", ]

all_genes_6c <- unique(modules_6c$gene_name)
gene_map_6c <- bitr(all_genes_6c, fromType = "SYMBOL", toType = "ENTREZID",
                    OrgDb = org.Bt.eg.db)

go_list_6c <- list()
for (mod in unique(modules_6c$module)) {
  mod_genes <- modules_6c$gene_name[modules_6c$module == mod]
  mod_entrez <- gene_map_6c$ENTREZID[gene_map_6c$SYMBOL %in% mod_genes]
  if (length(mod_entrez) < 5) next
  ego <- enrichGO(gene = mod_entrez, OrgDb = org.Bt.eg.db, ont = "BP",
    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
  if (!is.null(ego) && nrow(ego@result[ego@result$p.adjust < 0.05, ]) > 0) {
    df_go <- ego@result[ego@result$p.adjust < 0.05, ]
    df_go$module <- mod
    go_list_6c[[mod]] <- head(df_go, 5)
  }
}

go_df_6c <- bind_rows(go_list_6c)
go_df_6c <- go_df_6c %>%
  group_by(Description, module) %>%
  slice_min(p.adjust, n = 1, with_ties = FALSE) %>% ungroup()

go_df_6c$Term_clean <- beautify_term(go_df_6c$Description)

long_fix <- c(
  "Immunoglobulin production involved in immunoglobulin-mediated immune response" =
    "Ig production in Ig-mediated immunity",
  "Establishment of protein localization to endoplasmic reticulum" =
    "Protein localization to ER (estab.)",
  "Protein localization to endoplasmic reticulum" = "Protein localization to ER",
  "Purine ribonucleoside triphosphate biosynthesis" = "Purine rNTP biosynthesis",
  "Purine nucleoside triphosphate biosynthesis" = "Purine NTP biosynthesis",
  "Proton motive force-driven ATP synthesis" = "PMF-driven ATP synthesis",
  "Production of molecular mediator of immune response" = "Immune mediator production",
  "Antigen receptor-mediated signaling" = "Antigen receptor signaling",
  "Ribonucleoprotein complex biogenesis" = "RNP complex biogenesis",
  "Ribosomal large subunit biogenesis" = "Ribosomal large subunit biogen.",
  "Actin cytoskeleton organization" = "Actin cytoskeleton org.",
  "Intracellular protein transport" = "Intracellular protein transport"
)
for (old in names(long_fix)) {
  go_df_6c$Term_clean <- gsub(beautify_term(old), long_fix[old],
                               go_df_6c$Term_clean, fixed = TRUE)
}
go_df_6c$Term_clean <- ifelse(nchar(go_df_6c$Term_clean) > 40,
  paste0(substr(go_df_6c$Term_clean, 1, 37), "..."), go_df_6c$Term_clean)

mod_order_6c <- c("blue","turquoise","brown","green","yellow","red","black")
mod_order_6c <- mod_order_6c[mod_order_6c %in% unique(go_df_6c$module)]
mod_labels_6c <- tools::toTitleCase(mod_order_6c)
go_df_6c$Module <- factor(go_df_6c$module, levels = mod_order_6c, labels = mod_labels_6c)
go_df_6c <- go_df_6c[!is.na(go_df_6c$Module), ]

term_order_6c <- go_df_6c %>%
  group_by(Term_clean) %>%
  summarise(mean_p = mean(-log10(p.adjust)), .groups = "drop") %>%
  arrange(mean_p) %>% pull(Term_clean)
go_df_6c$Term_clean <- factor(go_df_6c$Term_clean, levels = term_order_6c)

fig6c <- ggplot(go_df_6c, aes(x = Module, y = Term_clean,
                               size = Count, color = p.adjust)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#B2182B", high = "#FDDBC7", name = "P.adjust") +
  scale_size_continuous(range = c(1.5, 6), name = "Gene\nCount") +
  scale_x_discrete(expand = c(0.08, 0.08)) +
  scale_y_discrete(expand = c(0.03, 0.03)) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "") +
  nature_compact_v2

ggsave(file.path(outdir, "Fig6c_module_GO_v3.pdf"),
       fig6c, width = 12, height = 10, units = "cm", dpi = 600, device = cairo_pdf)
cat("Fig 6c done!\n")

# 保存hdWGCNA对象
# Save the hdWGCNA object
saveRDS(obj2, file.path(rdsdir, "integrated_25cl_hdwgcna.rds"))


# ############################################################
# Fig 6d — Signaling flow comparison (rankNet, HS红/PFTN蓝)
# Fig 6d — signalling flow comparison (rankNet; HS red / PFTN blue)
# 输出: Fig6e_signaling_comparison.pdf (16×12cm)
# Output: Fig6e_signaling_comparison.pdf (16x12 cm)
# ############################################################

library(CellChat)

# 分成HS和PFTN
# Split into HS and PFTN
obj_hs   <- subset(obj2, treatment == "HS")
obj_pftn <- subset(obj2, treatment == "PFTN")

cc_hs   <- createCellChat(object = obj_hs,   group.by = "celltype", assay = "RNA")
cc_pftn <- createCellChat(object = obj_pftn, group.by = "celltype", assay = "RNA")

CellChatDB <- CellChatDB.human
cc_hs@DB   <- CellChatDB
cc_pftn@DB <- CellChatDB

process_cellchat <- function(cc) {
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  return(cc)
}

cc_hs   <- process_cellchat(cc_hs)
cc_pftn <- process_cellchat(cc_pftn)

object.list <- list(HS = cc_hs, PFTN = cc_pftn)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

pE <- rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
pE <- pE + scale_fill_manual(values = c("HS" = "#D73027", "PFTN" = "#4575B4"))

ggsave(file.path(outdir, "Fig6e_signaling_comparison.pdf"),
       pE, width = 16, height = 12, units = "cm")

cat("Fig 6d done!\n")


# ============================================================
# Console summary printed after all figure code has run.
cat("\n=== 全部Figure代码运行完毕 ===\n")
cat("输出目录:", outdir, "\n")
cat("Fig 1:  Gemini生成 (无R代码)\n")
cat("Fig 2a: Fig2a_PCA_final.pdf\n")
cat("Fig 2b: Fig2b_correlation_v7.png\n")
cat("Fig 2c: Fig2d_jitter_v2.pdf\n")
cat("Fig 2图例: Fig2_legend_combined.pdf\n")
cat("Fig 3:  Fig3_final_v7.pdf\n")
cat("Fig 4a: Fig4a_miRNA_heatmap_h_v3.pdf\n")
cat("Fig 4b: Fig4b_network_v10.pdf\n")
cat("Fig 4图例: Fig4_legend_combined_v2.pdf\n")
cat("Fig 5a: Fig5a_tsne_celltype.pdf\n")
cat("Fig 5b: Fig5c_cellprop_bar.pdf\n")
cat("Fig 5c: Fig5d_deg_counts.pdf\n")
cat("Fig 5d: Fig5b_dotplot_markers_v2.pdf\n")
cat("Fig 6a: Fig6a_dendrogram.pdf\n")
cat("Fig 6b: Fig6b_module_celltype_heatmap_v2.pdf\n")
cat("Fig 6c: Fig6c_module_GO_v3.pdf\n")
cat("Fig 6d: Fig6e_signaling_comparison.pdf\n")
cat("============================================================\n")
