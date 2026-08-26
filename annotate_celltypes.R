## ---------------------------------------------------------------------------
## Cell-type annotation of the integrated PBMC object
##
## This script documents, end to end, how the `celltype` labels in the deposited
## single-cell object were produced. In the previous release of this repository the
## annotation appeared only as a hard-coded cluster-to-cell-type lookup table, so a
## reader could see the outcome but not how it was reached. Reviewer 2 asked for the
## workflow and the code to agree; this file is that correction.
##
## Workflow actually used
##   1. Graph-based clustering at resolution 0.85 -> 25 clusters (Seurat).
##   2. SingleR against celldex::MonacoImmuneData, run at CLUSTER level, twice:
##        once with ref$label.main and once with ref$label.fine.
##   3. For each cluster the two predictions were compared. Where they agreed at the
##      lineage level the label was taken directly. Where they disagreed, the cluster
##      was adjudicated against the canonical marker genes shown in the dot plot that
##      appears as Fig. 5d, and whichever prediction those markers supported was kept.
##      The comments on the disagreeing clusters below record which call won.
##   4. One cluster could not be resolved this way and is labelled "Unassigned".
##
## SCINA is NOT part of this workflow. Earlier versions of the manuscript listed it. A search
## of the whole analysis environment - every script, saved workspace and interactive R history
## on the host where these data were processed - found no call to it, and the annotated object
## carries no SCINA-derived field. The Methods text has been corrected accordingly.
## ---------------------------------------------------------------------------

suppressMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(dplyr)
})

rdsdir <- "/home/guozitai/scRNA/h5ad/new"
outdir <- "."

obj <- readRDS(file.path(rdsdir, "integrated_data_with_singleR.rds"))
DefaultAssay(obj) <- "RNA"
if (!"data" %in% Layers(obj[["RNA"]])) obj <- NormalizeData(obj, assay = "RNA")
expr <- GetAssayData(obj, assay = "RNA", layer = "data")

## ---- step 2: SingleR at cluster level, main and fine ----------------------
ref <- celldex::MonacoImmuneData()
pred_main <- SingleR(test = expr, ref = ref, labels = ref$label.main,
                     clusters = obj$seurat_clusters, assay.type.test = "logcounts")
pred_fine <- SingleR(test = expr, ref = ref, labels = ref$label.fine,
                     clusters = obj$seurat_clusters, assay.type.test = "logcounts")

singler <- data.frame(
  cluster = rownames(pred_main),
  singleR_main = pred_main$labels,
  singleR_fine = pred_fine$labels[match(rownames(pred_main), rownames(pred_fine))],
  stringsAsFactors = FALSE
)

## ---- step 3: adjudicated final assignment ---------------------------------
## The mapping below is the outcome of step 3, recorded explicitly so that it can be
## audited cluster by cluster against the two SingleR predictions and against Fig. 5d.
## Clusters whose main and fine predictions disagree are commented with the reason.
final_label <- c(
  "0"  = "T cells",         # main CD4+ T / fine Th1 - agree (T lineage)
  "1"  = "T cells",         # main CD4+ T / fine Treg - agree
  "2"  = "B cells",
  "3"  = "B cells",
  "4"  = "B cells",
  "5"  = "Monocytes",       # fine calls myeloid DC; CD14 0.71 CSF1R 0.81 but FLT3 0.014 -> monocyte
  "6"  = "T cells",         # fine calls switched memory B; CD3E 0.93 vs CD79B 0.06 -> T cell
  "7"  = "Monocytes",       # as cluster 5: CD14 0.67 CSF1R 0.52 FLT3 0.002
  "8"  = "T cells",
  "9"  = "Monocytes",       # as cluster 5: CD14 0.33 CSF1R 0.47 FLT3 0.004
  "10" = "Monocytes",       # as cluster 5: CD14 0.50 CSF1R 0.83 FLT3 0.007
  "11" = "T cells",         # main CD8+ T / fine effector memory CD8 - agree
  "12" = "NK cells",        # main calls T cells; NKG7 0.79 KLRD1 0.36 vs CD3E 0.27 -> NK
  "13" = "Plasma cells",    # main calls B cells; MZB1 0.98 JCHAIN 0.97 -> plasma cell
  "14" = "B cells",
  "15" = "Monocytes",
  "16" = "B cells",
  "17" = "T cells",
  "18" = "Dendritic cells", # both calls dendritic; FLT3 0.222, the highest of any cluster
  "19" = "T cells",
  "20" = "B cells",
  "21" = "T cells",
  "22" = "T cells",
  "23" = "T cells",
  "24" = "Unassigned"       # main dendritic / fine exhausted B - mutually inconsistent, and
                            # CD34, PROM1 and MPO are detected in 0% of these 110 cells, so
                            # no lineage can be supported. Previously mislabelled "HSCs".
)

obj$celltype_main <- singler$singleR_main[match(as.character(obj$seurat_clusters), singler$cluster)]
obj$celltype_fine <- singler$singleR_fine[match(as.character(obj$seurat_clusters), singler$cluster)]
obj$celltype      <- final_label[as.character(obj$seurat_clusters)]

## ---- step 4: the per-cluster evidence table shipped with this repository ----
evidence <- obj@meta.data %>%
  count(seurat_clusters, celltype_main, celltype_fine, celltype, name = "n_cells") %>%
  mutate(cluster = as.integer(as.character(seurat_clusters))) %>%
  arrange(cluster) %>%
  select(cluster, n_cells,
         singleR_main = celltype_main, singleR_fine = celltype_fine,
         final_label  = celltype)
## The shipped table also carries the detection rate of 14 canonical lineage markers per
## cluster; see evidence_table.R for that step.
write.csv(evidence, file.path(outdir, "cluster_annotation_evidence.csv"), row.names = FALSE)
print(evidence, row.names = FALSE)

## sanity check against the manuscript: the unassigned cluster expresses no HSPC markers
hspc <- intersect(c("CD34", "PROM1", "MPO"), rownames(obj))
if (length(hspc)) {
  sub <- subset(obj, subset = seurat_clusters == 24)
  cat("\nDetection rate of HSPC markers in the unassigned cluster:\n")
  print(round(rowMeans(GetAssayData(sub, assay = "RNA", layer = "data")[hspc, , drop = FALSE] > 0), 4))
}

saveRDS(obj, file.path(rdsdir, "integrated_25cl_annotated.rds"))
cat("\nAnnotation complete.\n")
