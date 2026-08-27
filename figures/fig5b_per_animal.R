## ---------------------------------------------------------------------------
## Fig. 5b - cell-type proportions per animal
##
## Reviewer 2 asked that cell-type proportions be shown at the individual-animal level
## rather than for cells pooled within each treatment group. This script produces the
## submitted panel: one stacked bar per animal, grouped by treatment. It supersedes the
## Fig. 5b section of fig5_panels_noun.R, which drew the pooled two-bar version.
##
## Fig5b_source_data.csv, written at the end, is the per-animal table behind the panel
## and is the same data as sheet S7 of the supplementary tables.
## ---------------------------------------------------------------------------

.libPaths(c("/home/guozitai/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
suppressMessages({library(dplyr); library(ggplot2)})
outdir <- "."

## ==================== Fig 5b : per-animal cell-type proportions ====================
md <- readRDS("/home/guozitai/20260409_scientific_data/data/meta_25cl.rds")
md$celltype2 <- ifelse(md$celltype=="HSCs","Unassigned",as.character(md$celltype))
md$treatment <- ifelse(grepl("^HS", md$orig.ident), "HS", "PFTN")
lv <- c("T cells","B cells","Monocytes","NK cells","Plasma cells","Dendritic cells","Unassigned")
pd <- md %>% count(orig.ident, treatment, celltype2) %>% group_by(orig.ident) %>%
  mutate(Proportion = 100*n/sum(n)) %>% ungroup() %>%
  mutate(celltype2 = factor(celltype2, levels=lv))
cols <- c("T cells"="#E64B35","B cells"="#4DBBD5","Monocytes"="#00A087","NK cells"="#3C5488",
          "Plasma cells"="#8491B4","Dendritic cells"="#F39B7F","Unassigned"="#BDBDBD")
p <- ggplot(pd, aes(x=orig.ident, y=Proportion, fill=celltype2)) +
  geom_col(width=0.72) +
  facet_grid(~treatment, scales="free_x", space="free_x") +
  scale_fill_manual(values=cols, name="Cell type") +
  scale_y_continuous(expand=expansion(mult=c(0,0.02))) +
  labs(x=NULL, y="Proportion of cells (%)") +
  theme_bw(base_size=9) +
  theme(panel.grid=element_blank(), strip.background=element_rect(fill="grey93",colour=NA),
        axis.text.x=element_text(angle=45,hjust=1), legend.key.size=unit(3.4,"mm"))
ggsave(file.path(outdir,"Fig5b_cellprop_per_animal.pdf"), p, width=4.5, height=3.2, device=cairo_pdf)
ggsave(file.path(outdir,"Fig5b_cellprop_per_animal.png"), p, width=4.5, height=3.2, dpi=600)
write.csv(pd, file.path(outdir,"Fig5b_source_data.csv"), row.names=FALSE)
cat("Fig 5b written\n")
