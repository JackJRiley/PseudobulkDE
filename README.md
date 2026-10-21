# PseudobulkDE

Pseudobulk differential expression from a Seurat object, using DESeq2 or edgeR.

Single-cell DE tests treat every cell as an independent sample, which
inflates p-values. PseudobulkDE sums counts per biological replicate first,
then tests with standard bulk RNA-seq methods.

## Install

```r
remotes::install_github("JackJRiley/PseudobulkDE")
```

## Usage

```r
library(PseudobulkDE)

res <- PseudobulkDE(seurat_obj,
                    assay = "H3K27ac",          # any assay with raw counts
                    contrast_by = "condition",
                    ident.1 = "stim", 
                    ident.2 = "ctrl",
                    sample_by = "replicate",
                    group_by = "celltype",      # optional
                    covariates = "batch",       # optional
                    method = "DESeq2",          # or "edgeR"
                    normalize = "library"       # or "ncells"
```

Positive `log2FC` means higher in `ident.1`.

`assay` can be any assay holding raw counts.

## sample_by vs group_by

These control two different things:

- **`sample_by` decides what gets summed together.** Cells with the same
  `sample_by` value (and the same `contrast_by` value) are summed into one
  pseudobulk sample. These are the replicates the DE method uses to
  estimate variance, so `sample_by` should identify independent biological
  replicates: mouse, donor, library, etc.
- **`group_by` decides how many DE tests are run.** The data is split by
  `group_by`, and a separate DE test is run for each level. With
  `group_by = "celltype"`, you get one stim vs ctrl test per cell type.
  Results are stacked into one table with a `group` column.

For example, with 2 conditions, 3 replicates and 2 cell types:

```
sample_by = "replicate", group_by = "celltype"

celltype A:  ctrl_r1 ctrl_r2 ctrl_r3 | stim_r1 stim_r2 stim_r3  ->  DE test 1
celltype B:  ctrl_r1 ctrl_r2 ctrl_r3 | stim_r1 stim_r2 stim_r3  ->  DE test 2
```

Without `group_by`, all cell types are summed together and one test is run.
This answers "what changes in the tissue overall", but differences in cell
type proportions between conditions can show up as DE. If you specifically
want to do this, it is recommended to downsample cells to a common n.

Don't put cell type in `sample_by`. That would treat each cell type as an
extra replicate of the same condition within a single test.

Adjusted p-values are corrected within each group's test, not across groups.
Groups with fewer than `min_samples` replicates (default = 2) per condition
are skipped with a warning.

## Normalisation

By default (`normalize = "library"`), DESeq2/edgeR normalise each pseudobulk
by library size, which assumes most genes don't change. A global shift, such
maybe due to a knockout of a core factor, is therefore normalised away, and
only relative changes remain.

Use `normalize = "ncells"` to scale each pseudobulk by its number of cells
instead. The test then asks whether the average cell has more or less signal,
which keeps global effects. This is only appropriate when all samples were
processed and sequenced together (e.g. same batch in split-pool BC) 
so reads per cell are comparable across samples.