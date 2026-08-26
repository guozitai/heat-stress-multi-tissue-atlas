.libPaths(c("/home/guozitai/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
suppressMessages({library(Seurat); library(dplyr); library(Matrix)})

obj <- readRDS("/home/guozitai/scRNA/h5ad/new/integrated_25cl_annotated.rds")
DefaultAssay(obj) <- "RNA"
if (!"data" %in% Layers(obj[["RNA"]])) obj <- NormalizeData(obj, assay = "RNA")
m  <- GetAssayData(obj, assay = "RNA", layer = "data")
md <- obj@meta.data
md$final_label <- ifelse(md$celltype == "HSCs", "Unassigned", as.character(md$celltype))

## canonical lineage markers, grouped by the lineage they define
panel <- list(
  T        = c("CD3D","CD3E"),
  B        = c("MS4A1","CD79B"),
  Plasma   = c("JCHAIN","MZB1"),
  Monocyte = c("CD14","CSF1R"),
  DC       = c("FLT3"),
  NK       = c("NKG7","KLRD1"),
  HSPC     = c("CD34","PROM1","MPO")
)
genes <- unlist(panel, use.names = FALSE)
genes <- genes[genes %in% rownames(m)]
cat("markers used:", paste(genes, collapse=", "), "\n")
cat("markers absent from the object:",
    paste(setdiff(unlist(panel, use.names=FALSE), genes), collapse=", "), "\n\n")

cl <- as.character(md$seurat_clusters)
ks <- sort(unique(as.integer(cl)))
det <- t(sapply(ks, function(k) round(Matrix::rowMeans(m[genes, which(cl == as.character(k)), drop = FALSE] > 0), 3)))
colnames(det) <- paste0("det_", genes)

base <- do.call(rbind, lapply(ks, function(k){
  x <- md[cl == as.character(k), ]
  data.frame(cluster = k, n_cells = nrow(x),
             singleR_main = names(sort(table(as.character(x$celltype_main)), decreasing = TRUE))[1],
             singleR_fine = names(sort(table(as.character(x$celltype_fine)), decreasing = TRUE))[1],
             final_label  = names(sort(table(as.character(x$final_label)),  decreasing = TRUE))[1],
             stringsAsFactors = FALSE)
}))
## does the label follow SingleR mechanically, or did the two levels disagree?
lineage <- function(s){
  s <- tolower(s)
  if (grepl("^t cells|cd4|cd8|th1|regulatory|effector memory", s)) return("T")
  if (grepl("plasmablast|plasma", s))       return("Plasma")
  if (grepl("b cells|memory b", s))         return("B")
  if (grepl("dendritic", s))                return("DC")
  if (grepl("monocyte", s))                 return("Monocyte")
  if (grepl("natural killer|nk", s))        return("NK")
  s
}
base$agreement <- ifelse(mapply(function(a,b) lineage(a) == lineage(b),
                                base$singleR_main, base$singleR_fine),
                         "concordant", "DISCORDANT - adjudicated on markers")

ev <- cbind(base, as.data.frame(det))
write.csv(ev, "/home/guozitai/20260409_scientific_data/data/cluster_annotation_evidence.csv", row.names = FALSE)
print(ev[, c("cluster","n_cells","singleR_main","singleR_fine","final_label","agreement")], row.names = FALSE)
cat("\n--- detection rates, discordant clusters only ---\n")
print(ev[ev$agreement != "concordant", c("cluster","final_label", colnames(det))], row.names = FALSE)
cat("\n=== EVIDENCE TABLE WRITTEN ===\n")
