# ============================================================
# Scientific Data 手稿 Figure R Code 汇总
# 项目：Multi-tissue mRNA-seq, miRNA-seq and single-cell RNA-seq
#       data from heat-stressed lactating Holstein cows
# 生成日期：2026-04-10
# ============================================================
#
# 运行环境：
#   OS: Linux (服务器)
#   R version: 4.x
#   额外library路径（每个脚本开头必须加）：
#     .libPaths(c(.libPaths(),
#                 "/home/jiangjun/R/x86_64-pc-linux-gnu-library/4.4",
#                 "/home/jiangjun/R/x86_64-pc-linux-gnu-library/4.1"))
#
# 已安装关键R包及版本：
#   Seurat 5.1.0, hdWGCNA 0.4.11, CellChat 2.2.0,
#   DESeq2, clusterProfiler, org.Bt.eg.db,
#   ComplexHeatmap, ggplot2, patchwork, cowplot,
#   readxl, dplyr, stringr, reshape2
#
# 数据路径：
#   Bulk count matrix: /home/guozitai/bulkRNA/20240806-guo-5tissues.xlsx (sheet 3)
#   Bulk DEG: /home/jiangjun/1231/DEG_results_all_tissues.rds
#   通路富集(up): /home/guozitai/20250220/up_data_final.csv
#   通路富集(down): /home/guozitai/20250220/down_data_final.csv
#   miRNA热图: /home/guozitai/20250220/heatmap0303.xlsx
#   miRNA靶基因: /home/guozitai/20250220/targetgene_enrich.xlsx
#   scRNA对象: /home/guozitai/scRNA/h5ad/new/integrated_25cl_hdwgcna.rds
#   CellChat对象: /home/guozitai/scRNA/h5ad/new/cellchat_merged.rds (不再使用)
#
# 输出目录：
#   outdir <- "/home/guozitai/20260409_scientific_data/figures"
#
# 全局作图标准（Nature系）：
#   字体：Helvetica或Arial
#   轴线宽：0.3pt → ggplot linewidth = 0.3/2.835
#   轴文字：6pt | 轴标题：7pt | 标题：8pt bold
#   图例标题：5.5pt | 图例文字：5pt
#   Panel标签：8pt bold 小写(a,b,c) — 在Illustrator中添加
#
# 配色方案：
#   组织：Muscle=#66C2A5, Adipose=#FC8D62, Rumen=#8DA0CB,
#         Liver=#E78AC3, Mammary=#A6D854
#   Treatment：HS=#D73027, PFTN=#4575B4
#   细胞类型：T cells=#E64B35, B cells=#4DBBD5, Monocytes=#00A087,
#             NK cells=#3C5488, Plasma cells=#F39B7F, DCs=#8491B4,
#             HSCs=#91D1C2
# ============================================================

.libPaths(c(.libPaths(),
            "/home/jiangjun/R/x86_64-pc-linux-gnu-library/4.4",
            "/home/jiangjun/R/x86_64-pc-linux-gnu-library/4.1"))

outdir <- "/home/guozitai/20260409_scientific_data/figures"


# ============================================================
# 1. Fig 2a: PCA — 五组织全局PCA (纵向 6x12 inch)
# ============================================================
# 输入：bulk count matrix (xlsx sheet 3)
# 输出：Fig2a_PCA_final.pdf
# ============================================================

library(ggplot2)
library(readxl)
library(DESeq2)
library(dplyr)

counts_raw <- read_excel("/home/guozitai/bulkRNA/20240806-guo-5tissues.xlsx", sheet = 3)
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
  Muscle  = "#66C2A5",
  Adipose = "#FC8D62",
  Rumen   = "#8DA0CB",
  Liver   = "#E78AC3",
  Mammary = "#A6D854"
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
      size = 3,
      color = c("#D73027", "#4575B4")
    ))
  )

pdf(file.path(outdir, "Fig2a_PCA_final.pdf"), width = 6, height = 12)
print(p)
dev.off()
cat("Fig2a done!\n")


# ============================================================
# 2. Fig 2 组合图例 (Tissue + Treatment + Spearman r)
# ============================================================
# 输入：无（纯图例）
# 输出：Fig2_legend_combined.pdf
# ============================================================

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
cat("Fig2 legend done!\n")


# ============================================================
# 3. Fig 4 组合图例 (Treatment + Tissue + Z-score + Regulation)
# ============================================================
# 输入：无（纯图例）
# 输出：Fig4_legend_combined_v2.pdf
# ============================================================

common_theme <- theme_void() +
  theme(
    legend.title      = element_text(size = 12, face = "bold"),
    legend.text       = element_text(size = 11),
    legend.key.width  = unit(8, "mm"),
    legend.key.height = unit(6, "mm"),
    legend.spacing.y  = unit(1, "mm")
  )

leg_trt <- get_legend(
  ggplot(data.frame(x=1:2, y=1:2, g=factor(c("HS","PFTN"), levels=c("HS","PFTN"))),
         aes(x, y, fill=g)) +
    geom_tile() +
    scale_fill_manual(values=c(HS="#D73027", PFTN="#4575B4"), name="Treatment") +
    common_theme
)

leg_tis <- get_legend(
  ggplot(data.frame(x=1:2, y=1:2, g=factor(c("Adipose","Muscle"), levels=c("Adipose","Muscle"))),
         aes(x, y, fill=g)) +
    geom_tile() +
    scale_fill_manual(values=c(Adipose="#FC8D62", Muscle="#66C2A5"), name="Tissue") +
    common_theme
)

leg_zs <- get_legend(
  ggplot(data.frame(x=1, y=1, z=0), aes(x, y, fill=z)) +
    geom_tile() +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                         name="Z-score", limits=c(-2,2), breaks=c(-2,-1,0,1,2)) +
    guides(fill = guide_colorbar(
      barwidth=unit(8,"mm"), barheight=unit(30,"mm"),
      title.position="top", ticks=TRUE
    )) +
    theme_void() +
    theme(legend.title=element_text(size=12, face="bold"),
          legend.text=element_text(size=11))
)

leg_reg <- get_legend(
  ggplot(data.frame(x=c(1,1), xend=c(2,2), y=c(2,1), yend=c(2,1),
                    g=factor(c("Up-regulated","Down-regulated"),
                             levels=c("Up-regulated","Down-regulated"))),
         aes(x=x, xend=xend, y=y, yend=yend, color=g)) +
    geom_segment(linewidth=1.5) +
    scale_color_manual(values=c("Up-regulated"="#C0392B","Down-regulated"="#2471A3"),
                       name="Regulation") +
    common_theme +
    theme(legend.key.width=unit(8,"mm"))
)

combined <- plot_grid(leg_trt, leg_tis, leg_zs, leg_reg,
                      ncol=1, rel_heights=c(1, 0.8, 1.4, 0.8), align="v")

pdf(file.path(outdir, "Fig4_legend_combined_v2.pdf"), width = 2, height = 5.5)
print(combined)
dev.off()
cat("Fig4 legend done!\n")


# ============================================================
# 4. Fig 4b: miRNA-target 网络图（配色更新版）
# ============================================================
# 输入：/home/guozitai/20250220/targetgene_enrich.xlsx
# 输出：Fig4b_network_v10.pdf
# 注意：仅替换颜色变量，其余代码不变
# ============================================================
# 在原始Fig4b代码中，替换以下6行颜色定义：

# col_adi_mir  <- "#FC8D62"    # 原: "#E8713A"
# col_adi_gene <- "#FDE0D0"    # 原: "#FDDCCA"
# col_adi_brd  <- "#E07040"    # 原: "#C85820"
# col_mus_mir  <- "#66C2A5"    # 原: "#4A7BB5"
# col_mus_gene <- "#D4EDE4"    # 原: "#D0E0F0"
# col_mus_brd  <- "#4AA88A"    # 原: "#2C5EA0"

# 其余代码（make_polar, draw_ring, draw_edges等函数）均不变，
# 完整代码见原始 Fig4b_network_v10 脚本。


# ============================================================
# 5. Fig 6b: Module-celltype 热图（去grey）
# ============================================================
# 输入：integrated_25cl_hdwgcna.rds (obj已在内存)
# 输出：Fig6b_module_celltype_heatmap_v2.pdf
# ============================================================

library(Seurat)
library(hdWGCNA)
library(reshape2)

# obj <- readRDS("/home/guozitai/scRNA/h5ad/new/integrated_25cl_hdwgcna.rds")
# ↑ 如果obj不在内存才需要运行

hMEs <- GetMEs(obj, harmonized = TRUE)
meta <- obj@meta.data
hMEs$celltype <- meta$celltype

# 去掉grey
hMEs <- hMEs[, !grepl("grey", colnames(hMEs)) | colnames(hMEs) == "celltype"]

mean_hME <- aggregate(. ~ celltype, data = hMEs, FUN = mean)
rownames(mean_hME) <- mean_hME$celltype
mean_hME$celltype <- NULL
colnames(mean_hME) <- gsub("^ME", "", colnames(mean_hME))

mod_order <- c("black","blue","brown","green","red","turquoise","yellow")
mod_order <- mod_order[mod_order %in% colnames(mean_hME)]
mean_hME <- mean_hME[, mod_order]

ct_order <- c("T cells","B cells","Monocytes","NK cells","Plasma cells","Dendritic cells","HSCs")
ct_order <- ct_order[ct_order %in% rownames(mean_hME)]
mean_hME <- mean_hME[ct_order, ]

df <- melt(as.matrix(mean_hME))
colnames(df) <- c("CellType", "Module", "MeanhME")

p6b <- ggplot(df, aes(x = Module, y = CellType, fill = MeanhME)) +
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
cat("Fig6b done!\n")


# ============================================================
# 6. Fig 6c: Module GO 气泡图（Fig3 V7风格，term缩短）
# ============================================================
# 输入：integrated_25cl_hdwgcna.rds (obj已在内存)
# 输出：Fig6c_module_GO_v3.pdf
# ============================================================

library(clusterProfiler)
library(org.Bt.eg.db)

modules <- GetModules(obj)
modules <- modules[modules$module != "grey", ]

all_genes <- unique(modules$gene_name)
gene_map <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Bt.eg.db)

go_list <- list()
for (mod in unique(modules$module)) {
  mod_genes <- modules$gene_name[modules$module == mod]
  mod_entrez <- gene_map$ENTREZID[gene_map$SYMBOL %in% mod_genes]
  if (length(mod_entrez) < 5) next
  
  ego <- enrichGO(
    gene          = mod_entrez,
    OrgDb         = org.Bt.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  if (!is.null(ego) && nrow(ego@result[ego@result$p.adjust < 0.05, ]) > 0) {
    df <- ego@result[ego@result$p.adjust < 0.05, ]
    df$module <- mod
    go_list[[mod]] <- head(df, 5)
  }
}

go_df <- bind_rows(go_list)
go_df <- go_df %>%
  group_by(Description, module) %>%
  slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
  ungroup()

# 美化term
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
  x <- gsub("^([a-z])", "\\U\\1", x, perl = TRUE)
  return(x)
}
go_df$Term_clean <- beautify_term(go_df$Description)

# 手动缩短长term
long_fix <- c(
  "Immunoglobulin production involved in immunoglobulin-mediated immune response" = "Ig production in Ig-mediated immunity",
  "Establishment of protein localization to endoplasmic reticulum" = "Protein localization to ER (estab.)",
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
  go_df$Term_clean <- gsub(beautify_term(old), long_fix[old], go_df$Term_clean, fixed = TRUE)
}
go_df$Term_clean <- ifelse(nchar(go_df$Term_clean) > 40,
                           paste0(substr(go_df$Term_clean, 1, 37), "..."),
                           go_df$Term_clean)

mod_order <- c("blue", "turquoise", "brown", "green", "yellow", "red", "black")
mod_order <- mod_order[mod_order %in% unique(go_df$module)]
mod_labels <- tools::toTitleCase(mod_order)
go_df$Module <- factor(go_df$module, levels = mod_order, labels = mod_labels)
go_df <- go_df[!is.na(go_df$Module), ]

term_order <- go_df %>%
  group_by(Term_clean) %>%
  summarise(mean_p = mean(-log10(p.adjust)), .groups = "drop") %>%
  arrange(mean_p) %>%
  pull(Term_clean)
go_df$Term_clean <- factor(go_df$Term_clean, levels = term_order)

# Fig3 V7 同款主题
nature_compact <- theme_minimal(base_family = "Helvetica", base_size = 7) +
  theme(
    axis.text.x       = element_text(angle = 45, hjust = 1, size = 6.5, color = "black"),
    axis.text.y       = element_text(size = 5.5, color = "black"),
    axis.title        = element_text(size = 7),
    axis.line         = element_line(colour = "black", linewidth = 0.3),
    axis.ticks        = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.title        = element_text(size = 8, face = "bold", hjust = 0.5),
    legend.title      = element_text(size = 5.5),
    legend.text       = element_text(size = 5),
    legend.key.size   = unit(0.2, "cm"),
    legend.margin     = margin(0, 0, 0, 0),
    plot.margin       = margin(2, 2, 2, 2, "mm")
  )

fig6c <- ggplot(go_df, aes(x = Module, y = Term_clean,
                            size = Count, color = p.adjust)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#B2182B", high = "#FDDBC7", name = "P.adjust") +
  scale_size_continuous(range = c(1.5, 6), name = "Gene\nCount") +
  scale_x_discrete(expand = c(0.08, 0.08)) +
  scale_y_discrete(expand = c(0.03, 0.03)) +
  coord_cartesian(clip = "off") +
  labs(x = "", y = "") +
  nature_compact

ggsave(file.path(outdir, "Fig6c_module_GO_v3.pdf"),
       fig6c, width = 12, height = 10, units = "cm", dpi = 600, device = cairo_pdf)
cat("Fig6c done!\n")


# ============================================================
# Nature通用主题函数（供其他脚本调用）
# ============================================================

theme_nature <- theme_classic() +
  theme(
    text             = element_text(family = "Arial"),
    axis.text        = element_text(size = 6, color = "black"),
    axis.title       = element_text(size = 7),
    plot.title       = element_text(size = 8, face = "bold"),
    legend.title     = element_text(size = 5.5),
    legend.text      = element_text(size = 5),
    legend.key.size  = unit(2.5, "mm"),
    axis.line        = element_line(linewidth = 0.3 / 2.835),
    axis.ticks       = element_line(linewidth = 0.3 / 2.835)
  )


cat("\n=== 全部代码运行完毕 ===\n")
cat("输出文件：\n")
cat("  Fig2a_PCA_final.pdf\n")
cat("  Fig2_legend_combined.pdf\n")
cat("  Fig4_legend_combined_v2.pdf\n")
cat("  Fig6b_module_celltype_heatmap_v2.pdf\n")
cat("  Fig6c_module_GO_v3.pdf\n")
cat("\n注意：Fig4b网络图仅需替换6个颜色变量，完整代码见原始脚本。\n")
cat("注意：Fig 1(BioRender/Gemini), Fig 2b(Spearman热图), Fig 2c(jitter),\n")
cat("      Fig 3(通路气泡图), Fig 4a(miRNA热图), Fig 5(scRNA全部)\n")
cat("      均由之前的独立脚本生成，不在本汇总中。\n")
