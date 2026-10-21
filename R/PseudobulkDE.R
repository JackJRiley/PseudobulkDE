#' Pseudobulk differential expression
#'
#' Sums raw counts into pseudobulk samples (one per unique combination of
#' `sample_by` + `contrast_by`), then tests `ident.1` vs `ident.2` with DESeq2
#' or edgeR.
#'
#' If `group_by` is set, a separate DE test is run for each level of
#' `group_by` (e.g. one test per cell type), and the results are stacked
#' with a `group` column. Without `group_by`, all cells go into one test.
#'
#' @param object A Seurat object.
#' @param assay Assay to run DE on, e.g. `"RNA"` or `"peaks"`.
#' @param contrast_by Metadata column to test, e.g. `"condition"`.
#' @param ident.1 Level of `contrast_by` to test.
#' @param ident.2 Reference level of `contrast_by`.
#' @param sample_by Metadata column(s) identifying biological replicates,
#'   e.g. `"replicate"`. Cells sharing these values (and the same
#'   `contrast_by` value) are summed into one pseudobulk sample. These are
#'   the samples DESeq2/edgeR use to estimate variance.
#' @param group_by Optional metadata column, e.g. `"celltype"`. The data is
#'   split by this column and a separate DE test is run for each level.
#' @param covariates Optional sample-level metadata columns to include in the
#'   design, e.g. `"batch"`. Design is `~ covariates + contrast_by`.
#' @param layer Layer holding raw counts. Default `"counts"`.
#' @param method `"DESeq2"` (default) or `"edgeR"`.
#' @param normalize How pseudobulk samples are scaled before testing.
#'   `"library"` (default) uses the method's own library-size normalisation,
#'   which assumes most genes don't change, so a global shift is normalised
#'   away. `"ncells"` scales each sample by its number of cells instead, so
#'   the test asks whether the average cell has more or less signal. Use
#'   `"ncells"` when you expect a global change (e.g. loss of a histone mark)
#'   and all samples were sequenced together.
#' @param min_cells Minimum cells for a pseudobulk sample to be kept.
#' @param min_samples Minimum pseudobulk samples needed in each of `ident.1`
#'   and `ident.2`. Groups below this are skipped.
#' @param verbose Print progress messages.
#'
#' @return A data.frame with columns `gene`, `log2FC` (ident.1 over ident.2),
#'   `p_val` and `p_val_adj`. If `group_by` was used, a `group` column says
#'   which test each row came from. Adjusted p-values are corrected within
#'   each group's test, not across groups.
#'
#' @importFrom stats model.matrix reformulate
#' @importFrom edgeR filterByExpr
#' @export
#' @examples
#' \dontrun{
#' res <- PseudobulkDE(seu, assay = "RNA", contrast_by = "condition",
#'                     ident.1 = "stim", ident.2 = "ctrl",
#'                     sample_by = "replicate", group_by = "celltype")
#' }
PseudobulkDE <- function(object,
                         assay,
                         contrast_by,
                         ident.1,
                         ident.2,
                         sample_by,
                         group_by = NULL,
                         covariates = NULL,
                         layer = "counts",
                         method = c("DESeq2", "edgeR"),
                         normalize = c("library", "ncells"),
                         min_cells = 10,
                         min_samples = 2,
                         verbose = TRUE) {
  method <- match.arg(method)
  normalize <- match.arg(normalize)
  meta <- object[[]]
  check_columns(meta, c(contrast_by, sample_by, group_by, covariates))
  if (!all(c(ident.1, ident.2) %in% meta[[contrast_by]])) {
    stop("ident.1 and ident.2 must both be values in '", contrast_by, "'.")
  }
  
  pb <- AggregatePseudobulk(object,
                            assay = assay,
                            sample_by = unique(c(sample_by, contrast_by, covariates)),
                            group_by = group_by,
                            layer = layer,
                            min_cells = min_cells)
  
  results <- list()
  for (g in names(pb)) {
    # keep only the two conditions being compared
    keep <- pb[[g]]$samples[[contrast_by]] %in% c(ident.1, ident.2)
    counts <- as.matrix(pb[[g]]$counts[, keep, drop = FALSE])
    samples <- pb[[g]]$samples[keep, , drop = FALSE]
    
    n1 <- sum(samples[[contrast_by]] == ident.1)
    n2 <- sum(samples[[contrast_by]] == ident.2)
    if (n1 < min_samples || n2 < min_samples) {
      warning("Skipping ", g, ": ", n1, " ", ident.1, " vs ", n2, " ", ident.2,
              " samples.", call. = FALSE)
      next
    }
    if (verbose) message(g, ": ", n1, " ", ident.1, " vs ", n2, " ", ident.2)
    
    # design data: covariates + contrast, with ident.2 as reference
    coldata <- samples[, covariates, drop = FALSE]
    colnames(coldata) <- make.names(colnames(coldata))
    coldata[] <- lapply(coldata, function(x) if (is.character(x)) factor(x) else x)
    coldata$contrast <- factor(ifelse(samples[[contrast_by]] == ident.1, "test", "ref"),
                               levels = c("ref", "test"))
    design <- model.matrix(reformulate(colnames(coldata)), coldata)
    
    # drop lowly expressed genes
    counts <- counts[filterByExpr(counts, design = design), , drop = FALSE]
    
    # NULL = method's own normalisation, otherwise scale by cell number
    n_cells <- if (normalize == "ncells") samples$n_cells else NULL
    
    res <- if (method == "DESeq2") {
      run_deseq2(counts, coldata, n_cells)
    } else {
      run_edger(counts, design, n_cells)
    }
    if (!is.null(group_by)) res <- data.frame(group = g, res)
    results[[g]] <- res
  }
  
  if (length(results) == 0) stop("No groups had enough samples to test.")
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}