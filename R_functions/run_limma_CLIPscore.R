run_limma_interaction_CLIPscore <- function(
    input_df,
    id_col = "gene_name",
    meta_df,
    condition_col = "Condition",
    timepoint_col = "Timepoint",
    batch_col = NULL,                  # <-- NEW
    handling = c("zeros", "filter"),
    voom = TRUE,
    min_reps = 2
) {
  suppressPackageStartupMessages({
    library(limma)
    library(edgeR)
    library(dplyr)
    library(tibble)
    library(genefilter)
  })
  
  handling <- match.arg(handling)
  
  # --- prepare numeric matrix (replicate CLIPscore columns only) ---
  expr_cols <- grep("_rep[0-9]+_CLIPscore$", colnames(input_df), value = TRUE)
  if (length(expr_cols) == 0) stop("No replicate *_rep#_CLIPscore columns found.")
  
  mat <- input_df %>%
    dplyr::select(all_of(id_col), all_of(expr_cols)) %>%
    tibble::column_to_rownames(id_col) %>%
    as.matrix()
  
  # --- handle negatives / filtering ---
  if (handling == "zeros") {
    mat[mat < 0] <- 0
  } else if (handling == "filter") {
    keep <- apply(mat, 1, function(x) any(x > 0, na.rm = TRUE))
    mat <- mat[keep, , drop = FALSE]
    mat[mat < 0] <- 0
  }
  
  # --- ensure metadata matches sample columns ---
  stopifnot(all(colnames(mat) %in% rownames(meta_df)))
  meta_df <- meta_df[colnames(mat), , drop = FALSE]
  
  # --- make sure required cols exist ---
  req <- c(condition_col, timepoint_col)
  if (!is.null(batch_col)) req <- c(req, batch_col)
  missing_cols <- setdiff(req, colnames(meta_df))
  if (length(missing_cols) > 0) {
    stop("meta_df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # --- enforce factors (important for model.matrix) ---
  meta_df[[condition_col]] <- factor(meta_df[[condition_col]])
  meta_df[[timepoint_col]] <- factor(meta_df[[timepoint_col]])
  if (!is.null(batch_col)) meta_df[[batch_col]] <- factor(meta_df[[batch_col]])
  
  # --- design formula ---
  # with batch: ~ 0 + Batch + Condition*Timepoint
  # without:    ~ 0 + Condition*Timepoint
  if (!is.null(batch_col)) {
    fml <- stats::as.formula(paste0("~ 0 + ", batch_col, " + ", condition_col, " * ", timepoint_col))
  } else {
    fml <- stats::as.formula(paste0("~ 0 + ", condition_col, " * ", timepoint_col))
  }
  
  design <- model.matrix(fml, data = meta_df)
  colnames(design) <- make.names(colnames(design))
  message("Design matrix columns: ", paste(colnames(design), collapse = ", "))
  
  # --- voom vs non-voom branch ---
  if (voom) {
    dge <- DGEList(counts = mat)
    dge <- calcNormFactors(dge)
    v <- voom(dge, design, plot = FALSE)
    fit <- lmFit(v, design)
  } else {
    fit <- lmFit(mat, design)
  }
  
  fit <- eBayes(fit)
  
  # --- identify interaction terms ---
  # Make.names introduces "."; model.matrix for ":" becomes "." too.
  # We match "Condition<level>.Timepoint<level>" or similar.
  interaction_terms <- grep(
    paste0("^", make.names(condition_col), ".*\\.", make.names(timepoint_col)),
    colnames(fit$coefficients),
    value = TRUE
  )
  message("Detected interaction terms: ", paste(interaction_terms, collapse = ", "))
  
  if (length(interaction_terms) > 0) {
    
    if (length(interaction_terms) == 1) {
      message("Single interaction term detected — using per-term stats for global interaction.")
      
      single_term <- interaction_terms[1]
      
      term_tab <- limma::topTable(fit, coef = single_term, number = Inf, sort.by = "P")
      term_tab <- tibble::rownames_to_column(term_tab, var = id_col)
      term_tab$contrast <- single_term
      
      global_interaction <- term_tab %>%
        dplyr::select(all_of(id_col), P.Value, adj.P.Val)
      
      per_term <- term_tab
      
    } else {
      message("Multiple interaction terms detected — performing global F-test across all.")
      
      fit_sub <- fit
      fit_sub$coefficients <- fit$coefficients[, interaction_terms, drop = FALSE]
      fit_sub$stdev.unscaled <- fit$stdev.unscaled[, interaction_terms, drop = FALSE]
      
      Fstat <- limma::topTable(fit_sub, coef = NULL, number = Inf, sort.by = "none")
      global_interaction <- tibble::rownames_to_column(Fstat, var = id_col) %>%
        dplyr::select(all_of(id_col), P.Value, adj.P.Val)
      
      term_tables <- lapply(interaction_terms, function(term) {
        tab <- limma::topTable(fit, coef = term, number = Inf, sort.by = "P")
        tab <- tibble::rownames_to_column(tab, var = id_col)
        tab$contrast <- term
        tab
      })
      
      per_term <- dplyr::bind_rows(term_tables)
      
      per_term <- dplyr::left_join(
        per_term,
        global_interaction %>% dplyr::select(all_of(id_col), adj.P.Val),
        by = id_col,
        suffix = c("", ".global")
      )
    }
    
  } else {
    warning("No interaction terms detected — check design matrix.")
    global_interaction <- data.frame()
    per_term <- data.frame()
  }
  
  return(list(
    global_interaction = global_interaction,
    per_term = per_term,
    fit = fit,
    design = design,
    meta = meta_df
  ))
}