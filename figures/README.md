# Figure re-rendering scripts

These scripts produced the **final versions of the panels submitted with the
manuscript**. They supersede the corresponding plotting sections of
`../SciData_Master_RCode_Final.R` for those panels; every other panel comes from
the master script unchanged.

| Script | Produces | What changed relative to the master script |
|---|---|---|
| `fig2c_rerun.R` | `Fig2_panel_c_DEG_jitter.pdf` | Font rendering fix; `set.seed(42)` added so the jitter is reproducible |
| `fig3_rerun.R` | `Fig3_ab_complete.pdf` | Font rendering fix and corrected page box (the earlier export was clipped) |
| `fig5_panels_noun.R` | `Fig5_panel_a_tSNE.pdf`, `_b_cell_proportion`, `_c_DEG_counts`, `_d_marker_dotplot` | Font rendering; unassigned cluster excluded |
| `fig6b_noun.R` | `Fig6_panel_b_module_heatmap.pdf` | Font rendering; unassigned cell type excluded |
| `fig5b_per_animal.R` | `Fig5b_cellprop_per_animal.pdf/.png`, `Fig5b_source_data.csv` | **Revision.** Cell-type proportions per animal (six stacked bars) instead of pooled within treatment group (two bars), as Reviewer 2 asked. Supersedes the Fig. 5b panel of `fig5_panels_noun.R` |
| `fig4b_network.R` | `Fig4b_network.pdf/.png` | **Revision.** miRNA-target network. Two changes: target genes are selected as distinct miRNA-gene pairs (the source table lists each pair once per enriched pathway, so a row-wise selection could draw one gene as several nodes), and every node is labelled - a gene targeted by two miRNAs is drawn, and labelled, once per miRNA. Labels are kept horizontal with automatic collision resolution: every label box is tested against every other and colliding labels are displaced vertically until they clear. The script reports the number of residual overlapping pairs, which is zero |

Panels not listed here (Fig. 2a, the Fig. 2 and Fig. 4 legends, Fig. 4a, Fig. 6a,
Fig. 6c, Fig. 6d) were not re-rendered and come from the master script.
Fig. 1 is a schematic with no R code.

## Two things worth knowing

**1. Text rendering.** `showtext` must be enabled. Without it, `cairo_pdf` lays out
the text itself and rounds every glyph advance to a whole point, which at 5 pt makes
letter spacing visibly uneven (Arial `r` comes out 20% too wide, `t` 28% too narrow,
and `i` and `t` end up the same width). The original panels were produced with
`showtext` on, so the text is outlined. A one-line check:

```bash
pdftotext panel.pdf - | wc -c    # 0 characters = outlined = correct
```

`showtext_opts(dpi=)` must match `ggsave(dpi=)` or font sizes come out wrong.
Fig. 2 and Fig. 3 use `base_family = "Helvetica"`, which resolves to Nimbus Sans
(the URW Helvetica clone) on the analysis machine; `showtext` is pointed at that
same font file to reproduce the original exactly. Fig. 5 and Fig. 6 use Arial.

**2. The unassigned cluster.** One small cluster (110 cells, 0.13% of 82,013) matched
no known cell type — CD34, PROM1 and MPO were all 0% positive, and its transcriptome
was a cross-lineage mixture. It is excluded from the figures, which therefore show
**81,903 cells**, while the manuscript reports the full 82,013 cells and describes
this cluster in the Technical Validation section. Excluding it shifts each cell-type
proportion by at most 0.07 percentage points.

In `fig6b_noun.R` the filtering must be applied to the hME table, not to the Seurat
object: `GetMEs()` returns the module eigengenes stored by hdWGCNA, which still have
82,013 rows, so subsetting the object first causes a length mismatch.

## Running

```bash
export HS_DATA_DIR=/path/to/data      # default: ./data
export HS_OUT_DIR=/path/to/figures    # default: ./figures
Rscript figures/fig5_panels_noun.R
```

These scripts read the same archived intermediate objects as the main script (see
"What this code does and does not run" in the top-level README); the bulk count matrix
is archived at figshare, <https://doi.org/10.6084/m9.figshare.32993768>.
Package versions are in `../sessionInfo.txt`.

## A note on the comments

The original working comments are in Chinese and have been left unchanged. English
translations were added afterwards and are marked `[EN]`.
