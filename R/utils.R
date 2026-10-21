# Pull raw counts, joining split Seurat v5 layers if needed
#' @importFrom SeuratObject Assays Layers JoinLayers LayerData GetAssayData
get_counts <- function(object, assay, layer) {
  if (!assay %in% Assays(object)) {
    stop("Assay '", assay, "' not found. Available: ",
         paste(Assays(object), collapse = ", "))
  }
  a <- object[[assay]]
  if (inherits(a, "Assay5")) {
    if (length(Layers(a, search = layer)) > 1) a <- JoinLayers(a, layers = layer)
    counts <- LayerData(a, layer = layer)
  } else {
    counts <- GetAssayData(a, layer = layer)
  }
  vals <- if (inherits(counts, "dgCMatrix")) counts@x else counts
  if (any(vals %% 1 != 0)) {
    stop("Layer '", layer, "' is not raw counts. Use the counts layer.")
  }
  counts
}

check_columns <- function(meta, cols) {
  missing <- setdiff(cols, colnames(meta))
  if (length(missing) > 0) {
    stop("Metadata column(s) not found: ", paste(missing, collapse = ", "))
  }
}

# DESeq2 Wald test, test vs ref
# If n_cells is given, size factors are proportional to cell number
#' @importFrom DESeq2 DESeqDataSetFromMatrix DESeq results "sizeFactors<-"
#' @importFrom stats reformulate
run_deseq2 <- function(counts, coldata, n_cells = NULL) {
  storage.mode(counts) <- "integer"
  dds <- DESeqDataSetFromMatrix(counts, colData = coldata,
                                design = reformulate(colnames(coldata)))
  if (!is.null(n_cells)) sizeFactors(dds) <- n_cells / exp(mean(log(n_cells)))
  dds <- DESeq(dds, quiet = TRUE)
  r <- as.data.frame(results(dds, contrast = c("contrast", "test", "ref")))
  r <- r[order(r$pvalue), ]
  data.frame(gene = rownames(r),
             log2FC = r$log2FoldChange,
             p_val = r$pvalue,
             p_val_adj = r$padj,
             baseMean = r$baseMean)
}

# edgeR quasi-likelihood F-test, test vs ref
# If n_cells is given, library sizes are proportional to cell number
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
run_edger <- function(counts, design, n_cells = NULL) {
  if (is.null(n_cells)) {
    y <- DGEList(counts)
    y <- calcNormFactors(y)
  } else {
    y <- DGEList(counts, lib.size = n_cells * mean(colSums(counts)) / mean(n_cells))
  }
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  r <- topTags(glmQLFTest(fit, coef = "contrasttest"), n = Inf)$table
  data.frame(gene = rownames(r),
             log2FC = r$logFC,
             p_val = r$PValue,
             p_val_adj = r$FDR,
             logCPM = r$logCPM)
}