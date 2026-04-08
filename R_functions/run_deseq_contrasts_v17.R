run_deseq_contrasts <- function(
    dds,
    my_comparisons,
    Outdirectory,
    tag,
    var1 = NULL,              # main variable (e.g., Condition)
    var2 = NULL,              # optional second variable (e.g., Timepoint)
    ref1 = NULL,              # reference level for var1
    ref2 = NULL,              # reference level for var2
    interaction = FALSE,      # if TRUE, build design with interaction
    pca_label_col = NULL,
    group_col = "Group",
    design_formula = NULL,         # full DESeq2 design (e.g. ~ Condition * Timepoint)
    use_lrt = FALSE,               # run LRT instead of Wald test
    reduced_formula = NULL,        # reduced model for LRT
    anno_df = NULL,
    dds_id_col = NULL,
    anno_id_col = NULL,
    anno_gene_col = "gene_name",
    collapse_method = "none",
    filter_mode = "none",
    filter_threshold = 1,
    one_high_threshold = NULL,
    filter_counts_type = "raw",
    min_reps = "all",
    sig_metric = "padj",
    sig_threshold = 0.05,
    n_label_up = 30,
    n_label_down = 30,
    do_pca = TRUE,
    PCA_ntop = 500,
    do_go = TRUE,
    OrgDb = OrgDb,
    make_venn = TRUE,
    venn_include = "Up and down",
    save_deseq_tables = TRUE,
    save_volcano = TRUE,
    outfile_tag = ""
) {
  
  suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(ggrepel)
    library(dplyr)
    library(tibble)
    library(VennDetail)
    library(clusterProfiler)
    library(org.Mm.eg.db)
  })
  
  # object to store LRT results for return
  lrt_df <- NULL
  # --- Annotation subsetting ---
  if (!is.null(anno_df) && nrow(anno_df) > 0 &&
      !is.null(dds_id_col) && !is.null(anno_id_col)) {
    
    message("Subsetting annotation to IDs present in both anno_df and dds...")
    dds_ids <- rownames(dds)
    anno_df <- anno_df[anno_df[[anno_id_col]] %in% dds_ids, , drop = FALSE]
    dds <- dds[dds_ids %in% anno_df[[anno_id_col]], ]
    message("After subsetting: ", nrow(anno_df), " annotation rows, ",
            nrow(dds), " rows in dds.")
  }
  
  # --- Output directories ---
  if (isTRUE(save_deseq_tables)) {
    dir.create(file.path(Outdirectory, "DESeq"), showWarnings = FALSE, recursive = TRUE)
  }
  
  if (isTRUE(save_volcano) &&
      (isTRUE(use_lrt) || (!is.null(my_comparisons) && nrow(my_comparisons) > 0))) {
    dir.create(file.path(Outdirectory, "Volcano"), showWarnings = FALSE, recursive = TRUE)
  }
  
  if (isTRUE(do_pca)) {
    dir.create(file.path(Outdirectory, "PCA"), showWarnings = FALSE, recursive = TRUE)
  }
  
  if (isTRUE(do_go)) {
    dir.create(file.path(Outdirectory, "GOanalysis"), showWarnings = FALSE, recursive = TRUE)
  }
  
  if (!isTRUE(use_lrt) && isTRUE(make_venn)) {
    dir.create(file.path(Outdirectory, "Venn"), showWarnings = FALSE, recursive = TRUE)
  }
  
  # --- Collapse transcripts ---
  if (collapse_method == "max_mean" &&
      !is.null(anno_df) && nrow(anno_df) > 0 &&
      !is.null(dds_id_col) && !is.null(anno_id_col)) {
    
    message("Collapsing transcripts to highest mean expression per gene...")
    counts_mat <- counts(dds)
    means <- rowMeans(counts_mat)
    anno_merge <- anno_df %>%
      dplyr::select(all_of(anno_id_col), all_of(anno_gene_col))
    merge_df <- data.frame(id = rownames(dds), mean_expr = means) %>%
      left_join(anno_merge, by = setNames(anno_id_col, "id"))
    keep_rows <- merge_df %>%
      group_by(.data[[anno_gene_col]]) %>%
      arrange(desc(mean_expr), id) %>%
      slice_head(n = 1) %>%
      pull(id)
    dds <- dds[rownames(dds) %in% keep_rows, ]
    message("Collapsed to ", nrow(dds), " rows (highest mean per gene).")
    
    # --- Keep annotation synchronized ---
    if (!is.null(anno_df) && nrow(anno_df) > 0 && !is.null(anno_id_col)) {
      anno_df <- anno_df[anno_df[[anno_id_col]] %in% rownames(dds), , drop = FALSE]
      message("Annotation table trimmed to ", nrow(anno_df), " rows after collapsing.")
    }
  }
  
  # --- Filtering ---
  if (filter_mode != "none") {
    message("Filtering genes: mode = ", filter_mode,
            ", threshold = ", filter_threshold,
            ", counts type = ", filter_counts_type)
    
    # Normalize input
    norm_flag <- tolower(filter_counts_type) %in% c("norm", "normalized", "normed")
    C <- if (norm_flag) counts(dds, normalized = TRUE) else counts(dds)
    
    grp <- colData(dds)[[group_col]]
    if (!is.factor(grp)) grp <- factor(grp)
    
    if (is.null(one_high_threshold)) one_high_threshold <- filter_threshold
    
    keep <- switch(
      filter_mode,
      "all_samples_atleast" = rowSums(C >= filter_threshold) == ncol(C),
      
      "any_group_all" = {
        pass <- rep(FALSE, nrow(C))
        for (g in levels(grp)) {
          smp <- rownames(colData(dds))[grp == g]
          if (length(smp) == 0) next
          pass <- pass | (rowSums(C[, smp, drop = FALSE] >= filter_threshold) == length(smp))
        }
        pass
      },
      
      "any_group_atleast" = {
        pass <- rep(FALSE, nrow(C))
        for (g in levels(grp)) {
          smp <- rownames(colData(dds))[grp == g]
          if (length(smp) == 0) next
          mr <- if (identical(min_reps, "all")) length(smp) else as.integer(min_reps)
          pass <- pass | (rowSums(C[, smp, drop = FALSE] >= filter_threshold) >= mr)
        }
        pass
      },
      
      "global_atleast" = {
        mr <- if (identical(min_reps, "all")) ncol(C) else as.integer(min_reps)
        rowSums(C >= filter_threshold) >= mr
      },
      
      "group_all_and_one" = {
        pass <- rep(FALSE, nrow(C))
        for (g in levels(grp)) {
          smp <- rownames(colData(dds))[grp == g]
          if (length(smp) == 0) next
          group_counts <- C[, smp, drop = FALSE]
          all_pass <- rowSums(group_counts >= filter_threshold) == length(smp)
          one_high <- matrixStats::rowMaxs(as.matrix(group_counts)) >= one_high_threshold
          pass <- pass | (all_pass & one_high)
        }
        pass
      },
      
      stop("Unknown filter_mode: ", filter_mode)
    )
    
    message("Keeping ", sum(keep), " of ", nrow(C), " genes (",
            round(100 * sum(keep) / nrow(C), 2), "%) after filtering.")
    dds <- dds[keep, ]
    message("Retained ", nrow(dds), " genes after filtering.")
    
    # --- Keep annotation synchronized ---
    if (!is.null(anno_df) && nrow(anno_df) > 0 && !is.null(anno_id_col)) {
      anno_df <- anno_df[anno_df[[anno_id_col]] %in% rownames(dds), , drop = FALSE]
      message("Annotation table trimmed to ", nrow(anno_df), " rows after filtering.")
    }
  }
  
  # --- Flexible design building ---
  if (!is.null(var1)) {
    message("Building design from variable(s): ", var1,
            if (!is.null(var2)) paste0(" and ", var2) else "")
    
    # Set reference levels if specified
    if (!is.null(ref1) && var1 %in% colnames(colData(dds))) {
      dds[[var1]] <- relevel(dds[[var1]], ref = ref1)
      message("Reference for ", var1, " set to: ", ref1)
    }
    if (!is.null(var2) && !is.null(ref2) && var2 %in% colnames(colData(dds))) {
      dds[[var2]] <- relevel(dds[[var2]], ref = ref2)
      message("Reference for ", var2, " set to: ", ref2)
    }
    
    if (!is.null(var2)) {
      if (isTRUE(interaction)) {
        design_formula <- as.formula(paste0("~ ", var1, " * ", var2))
        message("Using interaction design: ", deparse(design_formula))
      } else {
        design_formula <- as.formula(paste0("~ ", var1, " + ", var2))
        message("Using additive design: ", deparse(design_formula))
      }
    } else {
      design_formula <- as.formula(paste0("~ ", var1))
      message("Using single-factor design: ", deparse(design_formula))
    }
    
    design(dds) <- design_formula
  } else if (!is.null(design_formula)) {
    message("Using provided design formula: ", deparse(design_formula))
    design(dds) <- design_formula
  } else {
    message("var1/design_formula not provided: keeping existing design: ", deparse(design(dds)))
    # do not overwrite design(dds)
  }
  
  # --- Run DESeq ---
  message("Running DESeq2...")
  if (is.null(mcols(dds)$dispersion) || all(is.na(mcols(dds)$dispersion))) {
    if (isTRUE(use_lrt)) {
      if (is.null(reduced_formula)) stop("You must supply a reduced_formula when use_lrt = TRUE")
      message("Performing LRT (Likelihood Ratio Test)...")
      dds <- DESeq(dds, test = "LRT", reduced = reduced_formula)
    } else {
      dds <- DESeq(dds)
    }
  }
  
  # ============================
  #  LRT GLOBAL INTERACTION MODE
  # ============================
  if (isTRUE(use_lrt)) {
    message("Summarizing global interaction (LRT) results...")
    lrt_res <- results(dds)
    lrt_df  <- as.data.frame(lrt_res)
    if (!is.null(dds_id_col)) {
      lrt_df <- tibble::rownames_to_column(lrt_df, var = dds_id_col)
    }
    
    if (!is.null(anno_df) && nrow(anno_df) > 0 &&
        !is.null(dds_id_col) && !is.null(anno_id_col)) {
      lrt_df <- lrt_df %>%
        left_join(anno_df, by = setNames(anno_id_col, dds_id_col))
    }
    
    # order by adjusted p-value
    if ("padj" %in% colnames(lrt_df)) {
      lrt_df <- lrt_df[order(lrt_df$padj, na.last = TRUE), ]
    }
    top_lrt <- head(lrt_df, 100)
    
    if (isTRUE(save_deseq_tables)) {
      write.csv(
        lrt_df,
        file = file.path(Outdirectory, "DESeq",
                         paste0(tag, "_LRT_full_results_", outfile_tag, ".csv")),
        row.names = FALSE
      )
      }
    
    
    message("Saved LRT global interaction results (", nrow(lrt_df), " genes).")
    
    # --- Volcano + GO for LRT interaction ---
    if (!("log2FoldChange" %in% colnames(lrt_df))) {
      warning("LRT results have no 'log2FoldChange' column; skipping volcano/GO.")
    } else if (!(sig_metric %in% colnames(lrt_df))) {
      warning("LRT results have no column '", sig_metric,
              "'; please use sig_metric = 'padj' or 'pvalue'. Skipping volcano/GO.")
    } else {
      # remove NAs in significance metric
      lrt_df_clean <- lrt_df %>% filter(!is.na(.data[[sig_metric]])) %>% ungroup()
      
      contrast_name <- if (nzchar(outfile_tag)) {
        paste0("LRT_interaction_", outfile_tag)
      } else {
        "LRT_interaction"
      }
      
      up_genes <- lrt_df_clean %>%
        filter(log2FoldChange > 0, .data[[sig_metric]] < sig_threshold)
      down_genes <- lrt_df_clean %>%
        filter(log2FoldChange < 0, .data[[sig_metric]] < sig_threshold)
      
      # Base unlabeled volcano
      vol_lrt_unlabeled <- ggplot(lrt_df_clean,
                                  aes(x = log2FoldChange,
                                      y = -log10(.data[[sig_metric]]))) +
        geom_point(color = "grey", alpha = 0.5) +
        geom_point(data = up_genes,   color = "red3") +
        geom_point(data = down_genes, color = "dodgerblue") +
        theme_classic() +
        geom_hline(yintercept = -log10(sig_threshold),
                   linetype = 'dotted', colour = "darkgrey") +
        labs(
          title = contrast_name,
          y     = paste0("-log10(", sig_metric, ")"),
          x     = "Log2 Fold Change (LRT full model coefficient)"
        ) +
        annotate("text",
                 x = Inf, y = Inf,
                 label = paste0("Decreased (",
                                nrow(down_genes), ")"),
                 color = "dodgerblue", hjust = 2, vjust = 2) +
        annotate("text",
                 x = Inf, y = Inf,
                 label = paste0("Increased (",
                                nrow(up_genes), ")"),
                 color = "red3", hjust = 2, vjust = 4)
      print(vol_lrt_unlabeled)
      if (isTRUE(save_volcano)) {
        ggsave(
          file.path(Outdirectory, "Volcano",
                    paste0(tag, "_", contrast_name, "_volcano_",
                           sig_metric, "_unlabeled_", outfile_tag, ".pdf")),
          plot = vol_lrt_unlabeled, width = 7, height = 5
        )
      }
      
      # Labeled volcano (if gene names present)
      if (anno_gene_col %in% colnames(lrt_df_clean)) {
        vol_lrt_labeled <-  ggplot(lrt_df_clean,
                                   aes(x = log2FoldChange,
                                       y = -log10(.data[[sig_metric]]))) +
          geom_point(color = "grey", alpha = 0.5) +
          geom_point(data = up_genes,   color = "red3") +
          geom_point(data = down_genes, color = "dodgerblue") +
          theme_classic() +
          geom_hline(yintercept = -log10(sig_threshold),
                     linetype = 'dotted', colour = "darkgrey") +
          labs(
            title = contrast_name,
            y     = paste0("-log10(", sig_metric, ")"),
            x     = "Log2 Fold Change (LRT full model coefficient)"
          ) +
          # Upregulated labels
          ggrepel::geom_label_repel(
            data = up_genes %>% slice_min(order_by = .data[[sig_metric]], n = n_label_up),
            aes(label = .data[[anno_gene_col]]),
            size = 3.5,
            color = "red3",
            max.overlaps = Inf,nudge_x = 2,nudge_y=2
          ) +
          
          # Downregulated labels
          ggrepel::geom_label_repel(
            data = down_genes %>% slice_min(order_by = .data[[sig_metric]], n = n_label_down),
            aes(label = .data[[anno_gene_col]]),
            size = 3.5,
            color = "dodgerblue",
            max.overlaps = Inf,nudge_x = -2,nudge_y=2
          )
        print(vol_lrt_labeled)
        if (isTRUE(save_volcano)) {
          ggsave(
            file.path(Outdirectory, "Volcano",
                      paste0(tag, "_", contrast_name, "_volcano_",
                             sig_metric, "_labeled_", outfile_tag, ".pdf")),
            plot = vol_lrt_labeled, width = 7, height = 5
          )
        }
      } else {
        message("No column '", anno_gene_col,
                "' in LRT results → skipping labeled volcano.")
      }
      
      # --- GO for LRT up/down ---
      if (isTRUE(do_go) && (anno_gene_col %in% colnames(lrt_df_clean))) {
        # UP
        if (nrow(up_genes) >= 25) {
          message("Running GO for LRT UP genes (", nrow(up_genes), ")...")
          go_up <- tryCatch(
            enrichGO(
              gene          = unique(up_genes[[anno_gene_col]]),
              OrgDb         = OrgDb,
              keyType       = "SYMBOL",
              readable      = TRUE,
              ont           = "all",
              pvalueCutoff  = 0.05,
              qvalueCutoff  = 0.10
            ),
            error = function(e) NULL
          )
          if (!is.null(go_up) && inherits(go_up, "enrichResult") &&
              nrow(go_up@result) > 0) {
            write.csv(
              go_up@result,
              file.path(Outdirectory, "GOanalysis",
                        paste0(tag, "_", contrast_name,
                               "_GO_Up_results_", sig_metric, "_", outfile_tag, ".csv")),
              row.names = FALSE
            )
            GO_Up_plot <- ggplot(go_up@result %>%
                                   arrange(p.adjust) %>%
                                   slice_head(n = 20),
                                 aes(Count, Description, fill = p.adjust)) +
              geom_bar(stat = "identity") +
              theme_classic(base_size = 18) +
              labs(title = paste0("GO Upregulated (LRT): ", contrast_name),
                   y = NULL)
            ggsave(
              file.path(Outdirectory, "GOanalysis",
                        paste0(tag, "_", contrast_name,
                               "_GO_Up_plot_", sig_metric, "_", outfile_tag, ".pdf")),
              plot = GO_Up_plot, width = 12, height = 8, units = "in"
            )
          } else {
            message("⚠️ No enriched GO terms found for LRT UP genes.")
          }
        } else message("Skipping LRT UP GO: fewer than 25 up genes.")
        
        # DOWN
        if (nrow(down_genes) >= 25) {
          message("Running GO for LRT DOWN genes (", nrow(down_genes), ")...")
          go_down <- tryCatch(
            enrichGO(
              gene          = unique(down_genes[[anno_gene_col]]),
              OrgDb         = OrgDb,
              keyType       = "SYMBOL",
              readable      = TRUE,
              ont           = "all",
              pvalueCutoff  = 0.05,
              qvalueCutoff  = 0.10
            ),
            error = function(e) NULL
          )
          if (!is.null(go_down) && inherits(go_down, "enrichResult") &&
              nrow(go_down@result) > 0) {
            write.csv(
              go_down@result,
              file.path(Outdirectory, "GOanalysis",
                        paste0(tag, "_", contrast_name,
                               "_GO_Down_results_", sig_metric, "_", outfile_tag, ".csv")),
              row.names = FALSE
            )
            GO_Down_plot <- ggplot(go_down@result %>%
                                     arrange(p.adjust) %>%
                                     slice_head(n = 20),
                                   aes(Count, Description, fill = p.adjust)) +
              geom_bar(stat = "identity") +
              theme_classic(base_size = 18) +
              labs(title = paste0("GO Downregulated (LRT): ", contrast_name),
                   y = NULL)
            ggsave(
              file.path(Outdirectory, "GOanalysis",
                        paste0(tag, "_", contrast_name,
                               "_GO_Down_plot_", sig_metric, "_", outfile_tag, ".pdf")),
              plot = GO_Down_plot, width = 12, height = 8, units = "in"
            )
          } else {
            message("⚠️ No enriched GO terms found for LRT DOWN genes.")
          }
        } else message("Skipping LRT DOWN GO: fewer than 25 down genes.")
      } else if (isTRUE(do_go)) {
        message("No column '", anno_gene_col,
                "' in LRT results → skipping GO for LRT.")
      }
    } # end volcano/GO-if
  } # end use_lrt block
  
  # --- PCA ---
  # --- PCA ---
  if (isTRUE(do_pca) && ncol(dds) >= 2) {
    message("Generating PCA plot...")
    tx <- tryCatch(
      { vst(dds, blind = FALSE) },
      error = function(e) {
        message("vst() failed: ", e$message)
        message("→ Falling back to varianceStabilizingTransformation().")
        varianceStabilizingTransformation(dds, blind = FALSE)
      }
    )
    
    pca_df <- plotPCA(tx, intgroup = c(group_col), returnData = TRUE, ntop = PCA_ntop)
    percentVar <- round(100 * attr(pca_df, "percentVar"))
    
    if (!is.null(pca_label_col) && pca_label_col %in% colnames(colData(dds))) {
      pca_df$Label <- colData(dds)[[pca_label_col]][match(rownames(pca_df), rownames(colData(dds)))]
    } else {
      pca_df$Label <- rownames(pca_df)
    }
    
    # Optional: wrap very long labels (keeps them from being one long line)
    # You can tune width= e.g. 18/22 depending on how long your strings are.
    pca_df$Label_wrapped <- stringr::str_wrap(as.character(pca_df$Label), width = 22)
    
    pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = .data[[group_col]])) +
      geom_point(size = 4) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      coord_fixed() +
      theme_classic(base_size = 14)
    
    ggsave(
      file.path(Outdirectory, "PCA",
                paste0(tag, "_PCA_unlabeled_", group_col, "_", outfile_tag, ".pdf")),
      plot = pca_plot, width = 7, height = 5
    )
    
    # Expand plotting window so labels have space, and don't clip them
    x_rng <- range(pca_df$PC1, na.rm = TRUE)
    y_rng <- range(pca_df$PC2, na.rm = TRUE)
    x_pad <- diff(x_rng) * 0.15
    y_pad <- diff(y_rng) * 0.15
    if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
    if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1
    
    pca_plot_labeled <- pca_plot +
      coord_cartesian(
        xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
        ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad),
        clip = "off"
      ) +
      ggrepel::geom_text_repel(
        aes(label = Label_wrapped),
        size = 2.5,
        color = "black",
        box.padding = 0.6,
        point.padding = 0.5,
        min.segment.length = 0,   # always draw a segment if needed
        segment.alpha = 0.6,
        max.overlaps = Inf,
        force = 10,
        force_pull = 0.5,
        max.iter = 20000,
        seed = 1
      )
    
    ggsave(
      file.path(Outdirectory, "PCA",
                paste0(tag, "_PCA_labeled_", group_col, "_", outfile_tag, ".pdf")),
      plot = pca_plot_labeled, width = 7, height = 5
    )
  }
  
  # --- Initialize ---
  venn_sets_up <- list()
  venn_sets_down <- list()
  results_list <- list()
  
  # --- Skip pairwise contrast loop if LRT or no comparisons ---
  if (isTRUE(use_lrt) || is.null(my_comparisons) || nrow(my_comparisons) == 0) {
    message("use_lrt=TRUE or no contrasts provided → skipping pairwise contrasts, volcano, GO, and Venn steps for individual contrasts.")
  } else {
    # --- Pairwise contrast loop (unchanged) ---
    for (i in seq_len(nrow(my_comparisons))) {
      treat   <- as.character(my_comparisons$treat[i])
      control <- as.character(my_comparisons$control[i])
      contrast_name <- paste0(treat, "_vs_", control)
      message("Processing contrast: ", contrast_name)
      res <- results(dds, contrast = c(group_col, treat, control))
      res_df <- as.data.frame(res)
      if (!is.null(dds_id_col)) res_df <- tibble::rownames_to_column(res_df, var = dds_id_col)
      if (!is.null(anno_df) && nrow(anno_df) > 0 &&
          !is.null(dds_id_col) && !is.null(anno_id_col)) {
        res_df <- res_df %>% left_join(anno_df, by = setNames(anno_id_col, dds_id_col))
      }
      raw_counts <- counts(dds, normalized = FALSE)
      counts_df <- as.data.frame(raw_counts)
      if (!is.null(dds_id_col)) counts_df <- tibble::rownames_to_column(counts_df, var = dds_id_col)
      res_df <- res_df %>% left_join(counts_df, by = dds_id_col)
      res_df <- res_df[order(res_df[[sig_metric]], na.last = TRUE), ]
      res_df <- res_df %>% ungroup()
      if (isTRUE(save_deseq_tables)) {
        write.csv(
          res_df,
          file = file.path(Outdirectory, "DESeq",
                           paste0(tag, "_DESeq_results_", contrast_name, "_",
                                  sig_metric, "_", outfile_tag, ".csv")),
          row.names = FALSE
        )
      }
     
      results_list[[contrast_name]] <- res_df
      
      # Volcano + GO for pairwise (unchanged)
      up_genes <- res_df %>% filter(log2FoldChange > 0, .data[[sig_metric]] < sig_threshold)
      down_genes <- res_df %>% filter(log2FoldChange < 0, .data[[sig_metric]] < sig_threshold)
      print(nrow(up_genes))
      print(nrow(down_genes))
      vol_unlabeled <- ggplot(res_df, aes(x = log2FoldChange,
                                          y = -log10(.data[[sig_metric]]))) +
        geom_point(color = "grey", alpha = 0.5) +
        geom_point(data = up_genes, color = "red3") +
        geom_point(data = down_genes, color = "dodgerblue") +
        theme_classic() +
        geom_hline(yintercept = 1.3, linetype = 'dotted', colour = "darkgrey") +
        scale_y_continuous(trans = "log1p") +
        labs(title = contrast_name,
             y = paste0("-log10(", sig_metric, ")"),
             x = paste0("Log2 Fold Change: ", treat, " vs ", control)) +
        annotate("text",
                 x = Inf, y = Inf,
                 label = paste0("Decreased (",
                                nrow(res_df %>% filter(log2FoldChange < 0,
                                                       .data[[sig_metric]] < sig_threshold)), ")"),
                 color = "dodgerblue", hjust = 2, vjust = 2) +
        annotate("text",
                 x = Inf, y = Inf,
                 label = paste0("Increased (",
                                nrow(res_df %>%
                                       filter(log2FoldChange > 0,
                                              .data[[sig_metric]] < sig_threshold)), ")"),
                 color = "red3", hjust = 2, vjust = 4)
      print(vol_unlabeled)
      if (isTRUE(save_volcano)) {
        ggsave(
          file.path(Outdirectory, "Volcano",
                    paste0(tag, "_", contrast_name, "_volcano_", sig_metric,
                           "_unlabeled_", outfile_tag, ".pdf")),
          plot = vol_unlabeled, width = 7, height = 5
        )
      }
      vol_labeled <-  ggplot(res_df, aes(x = log2FoldChange,
                                         y = -log10(.data[[sig_metric]]))) +
        geom_point(color = "grey", alpha = 0.5) +
        geom_point(data = up_genes, color = "red3") +
        geom_point(data = down_genes, color = "dodgerblue") +
        theme_classic() +
        geom_hline(yintercept = 1.3, linetype = 'dotted', colour = "darkgrey") +
        scale_y_continuous(trans = "log1p") +
        labs(title = contrast_name,
             y = paste0("-log10(", sig_metric, ")"),
             x = paste0("Log2 Fold Change: ", treat, " vs ", control)) +
        # Upregulated labels
        ggrepel::geom_label_repel(
          data = up_genes %>% slice_min(order_by = .data[[sig_metric]], n = n_label_up),
          aes(label = .data[[anno_gene_col]]),
          size = 3.5,
          color = "red3",
          max.overlaps = Inf,nudge_x = 2,nudge_y=2
        ) +
        
        # Downregulated labels
        ggrepel::geom_label_repel(
          data = down_genes %>% slice_min(order_by = .data[[sig_metric]], n = n_label_down),
          aes(label = .data[[anno_gene_col]]),
          size = 3.5,
          color = "dodgerblue",
          max.overlaps = Inf,nudge_x = -2,nudge_y=2
        )
      if (isTRUE(save_volcano)) {
        ggsave(
          file.path(Outdirectory, "Volcano",
                    paste0(tag, "_", contrast_name, "_volcano_", sig_metric,
                           "_labeled_", outfile_tag, ".pdf")),
          plot = vol_labeled, width = 7, height = 5
        )
        }
      print(vol_labeled)
      # GO for pairwise (unchanged)
      if (isTRUE(do_go) && (anno_gene_col %in% names(res_df))) {
        # UP
        if (nrow(up_genes) >= 25) {
          message("Running GO for UP genes (", nrow(up_genes), ")...")
          go_up <- tryCatch(
            enrichGO(
              gene          = up_genes[[anno_gene_col]],
              OrgDb         = OrgDb,
              keyType       = "SYMBOL",
              readable      = TRUE,
              ont           = "all",
              pvalueCutoff  = 0.05,
              qvalueCutoff  = 0.10
            ),
            error = function(e) NULL
          )
          if (!is.null(go_up) && inherits(go_up, "enrichResult") &&
              nrow(go_up@result) > 0) {
            write.csv(
              go_up@result,
              file.path(Outdirectory, "GOanalysis",
                        paste0(tag,"_",contrast_name,"_GO_Up_results_",
                               sig_metric,"_",outfile_tag,".csv")),
              row.names = FALSE
            )
            GO_Up_plot <- ggplot(go_up@result %>%
                                   arrange(p.adjust) %>%
                                   slice_head(n = 20),
                                 aes(Count, Description, fill = p.adjust)) +
              geom_bar(stat = "identity") +
              theme_classic(base_size = 18) +
              labs(title = paste0("GO Upregulated: ", contrast_name), y = NULL)
            ggsave(file.path(Outdirectory, "GOanalysis",
                             paste0(tag,"_",contrast_name,"_GO_Up_plot_",
                                    sig_metric,"_",outfile_tag,".pdf")),
                   plot = GO_Up_plot, width = 12, height = 8, units = "in")
          } else {
            message("⚠️ No enriched GO terms found for UP genes (", contrast_name, ").")
          }
        } else message("Skipping GO analysis: fewer than 25 up genes")
        
        # DOWN
        if (nrow(down_genes) >= 25) {
          message("Running GO for DOWN genes (", nrow(down_genes), ")...")
          go_down <- tryCatch(
            enrichGO(
              gene          = down_genes[[anno_gene_col]],
              OrgDb         = OrgDb,
              keyType       = "SYMBOL",
              readable      = TRUE,
              ont           = "all",
              pvalueCutoff  = 0.05,
              qvalueCutoff  = 0.10
            ),
            error = function(e) NULL
          )
          if (!is.null(go_down) && inherits(go_down, "enrichResult") &&
              nrow(go_down@result) > 0) {
            write.csv(
              go_down@result,
              file.path(Outdirectory, "GOanalysis",
                        paste0(tag,"_",contrast_name,"_GO_Down_results_",
                               sig_metric,"_",outfile_tag,".csv")),
              row.names = FALSE
            )
            GO_Down_plot <- ggplot(go_down@result %>%
                                     arrange(p.adjust) %>%
                                     slice_head(n = 20),
                                   aes(Count, Description, fill = p.adjust)) +
              geom_bar(stat = "identity") +
              theme_classic(base_size = 18) +
              labs(title = paste0("GO Downregulated: ", contrast_name), y = NULL)
            ggsave(file.path(Outdirectory, "GOanalysis",
                             paste0(tag,"_",contrast_name,"_GO_Down_plot_",
                                    sig_metric,"_",outfile_tag,".pdf")),
                   plot = GO_Down_plot, width = 12, height = 8, units = "in")
          } else {
            message("⚠️ No enriched GO terms found for DOWN genes (", contrast_name, ").")
          }
        } else message("Skipping GO analysis: fewer than 25 down genes")
      }
      
      # Store for Venn
      if (anno_gene_col %in% names(up_genes)) {
        venn_sets_up[[contrast_name]] <- unique(up_genes[[anno_gene_col]])
      } else if ("gene" %in% names(up_genes)) {
        venn_sets_up[[contrast_name]] <- unique(up_genes$gene)
      } else {
        venn_sets_up[[contrast_name]] <- rownames(up_genes)
      }
      
      if (anno_gene_col %in% names(down_genes)) {
        venn_sets_down[[contrast_name]] <- unique(down_genes[[anno_gene_col]])
      } else if ("gene" %in% names(down_genes)) {
        venn_sets_down[[contrast_name]] <- unique(down_genes$gene)
      } else {
        venn_sets_down[[contrast_name]] <- rownames(down_genes)
      }
    } # end for
  } # end else (non-LRT pairwise)
  
  # --- Venn plots ---
  if (!isTRUE(use_lrt) && isTRUE(make_venn)) {
    set_list <- list()
    if (venn_include %in% c("Up", "Up and down")) {
      for (nm in names(venn_sets_up)) set_list[[paste0("Up_", nm)]] <- venn_sets_up[[nm]]
    }
    if (venn_include %in% c("Down", "Up and down")) {
      for (nm in names(venn_sets_down)) set_list[[paste0("Down_", nm)]] <- venn_sets_down[[nm]]
    }
    set_list <- set_list[vapply(set_list, function(x) length(x) > 0, logical(1))]
    
    if (length(set_list) >= 2) {
      ven_obj <- venndetail(set_list)
      png(file.path(Outdirectory, "Venn",
                    paste0(tag, "_UpSet_", outfile_tag, ".png")),
          width = 1200, height = 800)
      print(plot(ven_obj, type = "upset", cex = 8))
      dev.off()
      write.csv(
        result(ven_obj),
        file = file.path(Outdirectory, "Venn",
                         paste0(tag, "_UpSet_results_", outfile_tag, ".csv")),
        row.names = FALSE
      )
    }
  }
  
  message("=== DESeq2 pipeline complete for: ", tag, " ===")
  return(list(
    dds = dds,
    lrt_df = lrt_df,      # NULL unless use_lrt = TRUE
    results_list = results_list   # empty list if no pairwise contrasts
  ))
  
}
