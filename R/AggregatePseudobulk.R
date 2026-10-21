#' Aggregate a Seurat object into pseudobulk samples
#'
#' Sums raw counts across cells sharing the same `sample_by` values. If
#' `group_by` is set, cells are first split by it and each level is
#' aggregated separately, giving one pseudobulk matrix per group.
#'
#' @param object A Seurat object.
#' @param assay Assay to aggregate, e.g. `"RNA"` or `"peaks"`.
#' @param sample_by Metadata column(s) defining a pseudobulk sample,
#'   e.g. `c("condition", "replicate")`. Each unique combination becomes one
#'   column of the pseudobulk matrix.
#' @param group_by Optional metadata column to split by first, e.g.
#'   `"celltype"`. Each level gets its own pseudobulk matrix.
#' @param layer Layer holding raw counts. Default `"counts"`.
#' @param min_cells Minimum cells for a pseudobulk sample to be kept.
#'
#' @return A list with one entry per group (`"all"` if no `group_by`), each
#'   holding `counts` (genes x samples) and `samples` (sample metadata with
#'   `n_cells`).
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom stats complete.cases
#' @export
AggregatePseudobulk <- function(object,
                                assay,
                                sample_by,
                                group_by = NULL,
                                layer = "counts",
                                min_cells = 10) {
  counts <- get_counts(object, assay, layer)
  meta <- object[[]][colnames(counts), , drop = FALSE]
  check_columns(meta, c(sample_by, group_by))
  
  if (is.null(group_by)) {
    groups <- rep("all", nrow(meta))
  } else {
    groups <- as.character(meta[[group_by]])
  }
  
  out <- list()
  for (g in unique(groups[!is.na(groups)])) {
    cells <- which(groups == g & complete.cases(meta[sample_by]))
    
    # one sample per unique combination of sample_by values
    key <- do.call(paste, c(meta[cells, sample_by, drop = FALSE], sep = "_"))
    f <- factor(key, levels = unique(key))
    ind <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1,
                        dims = c(length(f), nlevels(f)),
                        dimnames = list(NULL, levels(f)))
    pb <- counts[, cells, drop = FALSE] %*% ind
    
    samples <- meta[cells[match(levels(f), key)], sample_by, drop = FALSE]
    rownames(samples) <- levels(f)
    samples$n_cells <- as.vector(table(f))
    
    keep <- samples$n_cells >= min_cells
    if (any(keep)) {
      out[[g]] <- list(counts = pb[, keep, drop = FALSE],
                       samples = samples[keep, , drop = FALSE])
    }
  }
  out
}