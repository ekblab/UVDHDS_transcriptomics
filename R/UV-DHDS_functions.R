## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************
## collection of functions used for the DHDS project

##*********************************************************************************************************
## plot_transcript_dist
## This function creates a bar plot showing the distribution of significantly altered transcripts 
## by transcript/genetype annotation (e.g. "protein-coding", "lncRNA", "miRNA").

## Arguments:
## data      – a data.frame containing transcript-level results (e.g. DESeq2 output with annotation columns)
## anno.col  – character string specifying the column name to group transcripts by (e.g. "GENETYPE_biomaRt")
plot_transcript_dist <- function(data, anno.col){
  n_cols <- length(unique(data[[anno.col]]))
  table(data[[anno.col]]) %>%
    as.data.frame() %>%
    mutate(Var1 = reorder(Var1, Freq, decreasing = TRUE)) %>%
    ggplot(aes(Var1, Freq, fill = Var1)) +
    geom_bar(stat = "identity", position = position_dodge(), color = "black", alpha = 0.7) +
    geom_text(aes(label = Freq), nudge_y = 200) +
    ylab("Number of significantly altered transcripts \n between Melanocytes and DSCs") +
    scale_fill_manual(values = rev(colorRampPalette(c("#F7FBFF","#C6DBEF", "#6BAED6", "#2171B5", "#08306B"))(n_cols))) +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_blank(),
          legend.title = element_blank()) +
    ggtitle(paste0(anno.col, " annotation"))
}

##*********************************************************************************************************
## extract_and_plot_venn
## This function identifies significantly altered transcripts between two cell types 
## (e.g. DSCs vs. Melanocytes), filters them by expression and optionally by biotype, 
## and visualizes overlaps via a Venn diagram. It also saves CSV files of exclusive and shared gene sets.

## Arguments:
## dds                – DESeqDataSet object with raw counts and metadata (used to extract sample groups)
## results_object     – DESeq2 results table with ENSEMBL IDs, p-values, and optional biotype columns (annotated with AnnoDBI or biomaRt)
## pThres             – adjusted p-value threshold for significance filtering (default: 0.05)
## expression_threshold – numeric count threshold for average expression to consider a gene expressed (default: 1)
## biotype_filter     – optional character vector to filter by gene biotypes (e.g. "lncRNA")
## output_dir         – path to save the result CSV files (default: "Results/cell_line")
## plot_title         – optional character string for the Venn diagram title

## Returns:
## A ggvenn object displaying overlap of significantly expressed genes across cell types.
extract_and_plot_venn <- function(dds, results_object, 
                                  pThres = 0.05, 
                                  expression_threshold = 1,
                                  biotype_filter = NULL,
                                  output_dir = "Results/cell_line", 
                                  plot_title = NULL) {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Step 1: Filter for significantly differentially expressed genes
  sig_counts <- assay(dds)[!is.na(results_object$padj) & results_object$padj < pThres, ]
  metadata <- colData(dds)
  
  # Step 2: Optional gene type (biotype) filtering using embedded annotation
  if (!is.null(biotype_filter)) {
    if ("GENETYPE_biomaRt" %in% colnames(results_object)) {
      gene_types <- results_object$GENETYPE_biomaRt
    } else if ("GENETYPE_AnnoDBI" %in% colnames(results_object)) {
      gene_types <- results_object$GENETYPE_AnnoDBI
    } else {
      stop("No gene type annotation column (GENETYPE_biomaRt or GENETYPE_AnnoDBI) found in results_object.")
    }
    
    valid_genes <- results_object$ENSEMBL[!is.na(gene_types) & gene_types %in% biotype_filter]
    sig_counts <- sig_counts[rownames(sig_counts) %in% valid_genes, ]
  }
  
  # Step 3: Prepare count subsets
  dsc_counts <- data.frame(sig_counts[, metadata$cell == "DSC"])
  melanocyte_counts <- data.frame(sig_counts[, metadata$cell == "Melanocytes"])
  
  dsc_counts$sum_counts <- rowMeans(dsc_counts)
  melanocyte_counts$sum_counts <- rowMeans(melanocyte_counts)
  
  # Step 4: Remove genes with low expression
  dsc_counts <- dsc_counts %>% filter(sum_counts > expression_threshold)
  melanocyte_counts <- melanocyte_counts %>% filter(sum_counts > expression_threshold)
  
  # Step 5: Prepare Venn inputs
  venn_input <- list(
    DSCs = rownames(dsc_counts),
    Melanocytes = rownames(melanocyte_counts)
  )
  
  # Step 6: Identify exclusive and common genes
  ex_DSC <- setdiff(venn_input$DSCs, venn_input$Melanocytes)
  ex_Mel <- setdiff(venn_input$Melanocytes, venn_input$DSCs)
  common_genes <- intersect(venn_input$DSCs, venn_input$Melanocytes)
  
  # Step 7: Save gene lists to CSV
  if (!is.null(biotype_filter)) {
    biotype_suffix <- paste(biotype_filter, collapse = "_")
    filename_suffix <- paste0("_", biotype_suffix)
  } else {
    filename_suffix <- ""
  }
  
  write.csv(data.frame(Gene = ex_DSC, ENSEMBL_input = paste0(ex_DSC, ",")), 
            file.path(output_dir, paste0("DSC_exclusive_transcripts", filename_suffix, ".csv")), 
            row.names = FALSE)
  
  write.csv(data.frame(Gene = ex_Mel, ENSEMBL_input = paste0(ex_Mel, ",")), 
            file.path(output_dir, paste0("Melanocytes_exclusive_transcripts", filename_suffix, ".csv")), 
            row.names = FALSE)
  
  write.csv(data.frame(Gene = common_genes, ENSEMBL_input = paste0(common_genes, ",")), 
            file.path(output_dir, paste0("DSC_Melanocytes_overlap_transcripts", filename_suffix, ".csv")), 
            row.names = FALSE)
  
  # Step 8: Plot Venn
  n_diff <- nrow(sig_counts)
  n_venn <- length(ex_DSC) + length(ex_Mel) + length(common_genes)
  
  p_venn <- ggvenn(venn_input, 
                   fill_color = c("#CD534CFF", "#0073C2FF"), 
                   stroke_size = 0.5, set_name_size = 7, text_size = 7) +
    theme(plot.caption = element_text(size = 12)) +
    labs(
      title = plot_title,
      caption = paste0(n_diff, " differentially expressed genes.\n",
                       "Exclusive genes defined as average expression ≤ ", expression_threshold, 
                       " counts in one cell type and > ", expression_threshold, " in the other.\n",
                       n_diff - n_venn, " genes were excluded due to low expression in both.")
    )
  
  return(p_venn)
}

##*********************************************************************************************************
## plot_biotype_heatmap
## This function generates a heatmap  of differentially expressed genes 
## from variance-stabilized RNA-seq data, optionally filtered by gene biotype 
## (e.g. "lncRNA", "miRNA", "protein-coding"). It includes annotations for cell type, time point, and donor.

## Arguments:
## dds              – DESeqDataSet object (used to extract sample metadata)
## vsd              – variance-stabilized DESeq2 object (assay(vsd) must contain transformed counts)
## results_object   – DESeq2 results table with ENSEMBL IDs and gene biotype annotations (columns: ENSEMBL, GENETYPE_biomaRt or GENETYPE_AnnoDBI)
## biotype_filter   – character vector specifying which biotypes to include (e.g. "lncRNA"); if NULL, all genes are used
## heatmap_fontsize – numeric, font size for heatmap annotations and legends (default: 18)
## plot_title       – optional character string for custom plot title
## show_row_dend    – logical indicating whether row dendrogram should be plotted. Default is FALSE
## show_row_names   – logical indicating whether row names should be plotted. Default is FALSE
## plot heatmap based on biotype (default is all transcripts)
plot_biotype_heatmap <- function(dds, vsd, results_object, 
                                 biotype_filter = NULL,
                                 heatmap_fontsize = 18,
                                 show_row_dend = FALSE,
                                 show_row_names = FALSE, 
                                 plot_title = NULL, ...) {
  
  # Step 1: Biotype filtering using ENSEMBL IDs
  if (!is.null(biotype_filter)) {
    if ("GENETYPE_biomaRt" %in% colnames(results_object)) {
      gene_types <- results_object$GENETYPE_biomaRt
    } else if ("GENETYPE_AnnoDBI" %in% colnames(results_object)) {
      gene_types <- results_object$GENETYPE_AnnoDBI
    } else {
      stop("Biotype filter requested, but no GENETYPE_biomaRt or GENETYPE_AnnoDBI column found.")
    }
    
    valid_genes <- results_object$ENSEMBL[!is.na(gene_types) & gene_types %in% biotype_filter]
    sig_genes <- intersect(results_object$ENSEMBL, valid_genes)
  } else {
    sig_genes <- results_object$ENSEMBL
  }
  
  if (length(sig_genes) == 0) {
    stop("No genes found for specified biotype filter.")
  }

  # Step 2: Prepare expression matrix from VST-transformed object (Control only)
  mat <- assay(vsd)[sig_genes, ]
  
  # Step 3: Z-score scaling across rows

  datScaled <- mat %>% t() %>% scale() %>% t()
  datScaled <- datScaled[complete.cases(datScaled) & rowSums(is.infinite(datScaled)) == 0, ]

  # Step 4: Annotation bar
  anno_df <- as.data.frame(colData(vsd)[, c("cell", "donor")])
  colnames(anno_df) <- c("Cell_type", "Donor")
  donors <- colorRampPalette(rev(brewer.pal(n = 7, name = "Set1")))(length(unique(anno_df$Donor)))
  # donors <- pal_npg("nrc")(length(unique(anno_df$Donor)))
  names(donors) <- unique(anno_df$Donor)
  
  annColors <- HeatmapAnnotation(
    df = anno_df,
    col = list(
      Cell_type = c("DSC" = ggsci::pal_npg("nrc")(4)[1],
                    "Melanocytes" = ggsci::pal_npg("nrc")(4)[2]),
      Donor = donors
    ),
    annotation_legend_param = list(
      Cell_type = list(nrow = 1), 
      Time = list(nrow = 1),
      title_gp = gpar(fontsize = heatmap_fontsize),
      labels_gp = gpar(fontsize = heatmap_fontsize)
    ),
    annotation_name_gp = gpar(fontsize = heatmap_fontsize, fontface = "bold")
  )
  
  # Step 5: Color function
  colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  col_fun <- colorRamp2(c(-2, 0, 2), c(colors[1], colors[51], colors[100]))
  
  # Step 6: Title and file suffix
  gene_count <- length(sig_genes)
  title_text <- ifelse(is.null(plot_title), paste0("Differentially Expressed Genes (", gene_count, ")"), plot_title)
  
  # Replace ENSEMBL rownames with gene symbols if available
  if ("SYMBOL" %in% colnames(results_object)) {
    ensembl_to_symbol <- results_object$SYMBOL
    names(ensembl_to_symbol) <- results_object$ENSEMBL
    
    # Ensure only genes present in both the matrix and the annotation
    common_genes <- intersect(rownames(datScaled), names(ensembl_to_symbol))
    
    # Replace ENSEMBL IDs with SYMBOLs for rownames
    rownames(datScaled)[rownames(datScaled) %in% common_genes] <- ensembl_to_symbol[common_genes]
  }

  # Step 7: Generate and save the heatmap
  Ht <- Heatmap(datScaled,
                col = col_fun,
                top_annotation = annColors,
                clustering_method_row = "average",
                clustering_method_columns = "average",
                clustering_distance_row = "pearson",
                clustering_distance_column = "euclidean",
                show_row_dend = show_row_dend,
                show_row_names = show_row_names,
                show_column_names = FALSE,
                use_raster = TRUE,
                column_names_gp = gpar(fontsize = 10),
                column_title = title_text,
                ...,
                heatmap_legend_param = list(
                  title = "row Z-score",
                  at = seq(-2, 2, by = 1),
                  color_bar = "continuous",
                  title_position = "topcenter",
                  legend_direction = "horizontal",
                  legend_width = unit(4, "cm")
                ))
  
  return(Ht)
}

##*********************************************************************************************************
## get_significance_annotations_from_list
## This function extracts significance labels (e.g. "*", "**", "***") for selected genes 
## from a list of DESeq2 result tables. It returns a data frame suitable for ggplot annotation 
## using `stat_pvalue_manual()`.

## Arguments:
## results_list – named list of DESeq2 results tables (with SYMBOL and padj columns)
## genes        – character vector of gene symbols for which significance annotations should be extracted

## Returns:
## A data frame with gene names, comparison groups, and significance labels
get_significance_annotations_from_list <- function(results_list, genes) {
  sig_df_list <- list()
  
  for (name in names(results_list)) {
    res <- results_list[[name]]
    
    # Try to extract: cell type, condition, timepoint
    # E.g., name = "Melanocytes_1x_UVB900_24h_vs_Melanocytes_control_24h"
    parts <- strsplit(name, "_vs_")[[1]]
    group1 <- parts[1]   # e.g., "Melanocytes_1x_UVB900_24h"
    group2 <- parts[2] # e.g., "Melanocytes_control_24h"
    
    # For each gene
    for (gene in genes) {
      row <- res[!is.na(res$SYMBOL) & res$SYMBOL == gene, ]
      if (nrow(row) == 0) next  # Skip if gene not found
      
      padj <- row$padj[1]
      label <- case_when(
        is.na(padj)      ~ "",
        padj < 0.001     ~ "***",
        padj < 0.01      ~ "**",
        padj < 0.05      ~ "*",
        TRUE             ~ ""
      )
      
      sig_df_list[[length(sig_df_list) + 1]] <- data.frame(
        gene = gene,
        group1 = group1,
        group2 = group2,
        label = label,
        y.position = NA  # Will be added later
      )
    }
  }
  
  sig_df <- do.call(rbind, sig_df_list)
  sig_df
}

##*********************************************************************************************************
## plot_gene_expression
## This function generates boxplots of normalized gene expression values 
## from a DESeq2 VST object, faceted by gene, and annotated with significance 
## asterisks based on a list of DESeq2 results.

## Arguments:
## vsd_obj   – DESeq2 VST (variance-stabilized transformation) object containing expression data
## anno_obj  – annotation data frame with SYMBOL and ENSEMBL columns
## res       – named list of DESeq2 results tables (same as used in get_significance_annotations_from_list)
## genes     – character vector of gene symbols to plot
## ncol      – number of columns in the facet plot (default: 1)

## Returns:
## A ggplot object showing boxplots of normalized expression with significance annotations
plot_gene_expression <- function(vsd_obj, anno_obj, res, genes, ncol = 1) {
  
  # Prepare expression matrix
  ind <- which(!is.na(anno_obj$SYMBOL) & anno_obj$SYMBOL %in% genes)
  vsd_counts <- assay(vsd_obj)
  df <- t(vsd_counts[ind, ])
  colnames(df) <- anno_obj$SYMBOL[match(colnames(df), anno_obj$ENSEMBL)]
  
  df <- cbind(colData(vsd_obj), df) %>%
    as.data.frame() %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "vsd")
  
  df <- df %>%
    mutate(
      cell = factor(cell),
      irradiation = factor(irradiation),
      condition = factor(condition),
      gene = factor(gene, levels = genes),
      cell_irradiation = factor(interaction(cell, irradiation, sep = "_", drop = TRUE))
    )
  
  sig_df <- get_significance_annotations_from_list(res, genes = genes)
  
  max_y <- df %>%
    group_by(gene, condition, time, cell) %>%
    summarise(y = max(vsd, na.rm = TRUE), .groups = "drop") %>% 
    ungroup() %>%
    group_by(gene, time) %>%
    mutate(y = max(y))
  
  sig_df <- sig_df %>%
    left_join(max_y, by = c("gene", "group2" = "condition")) %>%
    mutate(y.position = y + 0.3)
  
  cell_irrad_colors <- setNames(
    c(
      colorRampPalette(c("#CD534CFF", "#892822"))(4),  # Shades for DSC
      colorRampPalette(c("#0073C2FF", "#003366"))(4)   # Shades for Melanocytes
    )[1:nlevels(df$cell_irradiation)],
    levels(df$cell_irradiation)
  )
  
  ggplot(df, aes(x = condition, y = vsd)) +
    geom_boxplot(
      aes(fill = cell_irradiation),
      alpha = 0.8, width = 0.6, outlier.shape = NA, linewidth = 0.6, color = "black"
    ) +
    geom_point(
      aes(fill = cell_irradiation),
      position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.6),
      shape = 21, size = 1.6, stroke = 0.3, color = "black"
    ) +
    facet_wrap(~gene, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = cell_irrad_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) + 
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(size = 0.6, color = "black"),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(0.6, "cm"),
      strip.text = element_text(face = "bold", size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
    ) +
    ylab("Normalized expression") + 
    stat_pvalue_manual(
      data = sig_df,
      label = "label",
      x = "group1",     # Only show star above the irradiated condition
      y.position = "y.position",
      tip.length = 0,   # No vertical line
      size = 3.5,
      bracket.size = 0  # No bracket
    )
}

##*********************************************************************************************************
## extract_gene_list
## This function extracts lists of ENSEMBL gene IDs from a named list of DESeq2 results tables,
## based on combinations of specified conditions, cell type, and time point. It assumes a consistent
## naming convention in the list: "Cell_Condition_Time_vs_Cell_control_Time".

## Arguments:
## ls         – named list of DESeq2 results tables, where each table has an ENSEMBL column
## conditions – character vector of conditions to extract (e.g. "1x_UVB900")
## cell       – character string specifying the cell type (e.g. "Melanocytes")
## time       – character string specifying the time point (e.g. "24h")

## Returns:
## A list of ENSEMBL gene vectors corresponding to the specified comparisons
extract_gene_list <- function(ls, conditions, cell, time){
  vec <- unlist(sapply(conditions, function(condition){
    paste0(cell, "_",condition, "_", time, "_vs_", cell, "_control_", time)
  }))
  
  lapply(vec, function(x){
    ls[[x]]$ENSEMBL
  })
}

##*********************************************************************************************************
## plot_venn
## This function creates a Venn diagram visualizing the overlap of ENSEMBL gene sets
## extracted for different conditions within a specific cell type and time point.
## It uses the `extract_gene_list()` function to retrieve gene sets from a named list
## of DESeq2 results and plots the overlap with `ggvenn`.

## Arguments:
## ls         – named list of DESeq2 results tables, each containing an ENSEMBL column
## vsd_obj    – variance-stabilized DESeq2 object (not used inside the function but kept for compatibility)
## cell       – character string specifying the cell type (e.g., "Melanocytes")
## conditions – character vector of conditions to include in the Venn diagram (e.g., c("1x_UVB900", "3x_UVB300"))
## time       – character string specifying the time point (e.g., "24h")

## Returns:
## A ggplot object showing the Venn diagram of gene overlaps between specified conditions
plot_venn <- function(ls, vsd_obj, cell, conditions, time){
  
  venn_input <- extract_gene_list(ls = ls, conditions = conditions, cell = cell, time = time)
  names(venn_input) <- stringr::str_replace_all(conditions, "_", " ")
  
  ggvenn(venn_input, fill_color = c( "#CD534CFF", "#0073C2FF", "springgreen"),  stroke_size = 0.5, set_name_size = 4, text_size = 4) +
    ggtitle(paste(cell, time)) 
}


##*********************************************************************************************************
## parse_irradiation_input
## This function parses a character vector of irradiation condition strings formatted like 
## "1x_UVB300" or "3x_IR" and extracts components to construct a human-readable description.
## It separates the dose (e.g., "1x"), irradiation type (e.g., "UVB", "UVA", "IR"), and energy (e.g., "300") 
## to produce strings like "1x 300J/m² UVB" or "3x IR".

## Arguments:
## input_vec – character vector of irradiation condition strings to parse

## Returns:
## A character vector of parsed irradiation descriptions
parse_irradiation_input <- function(input_vec) {
  # Extract components using regex
  matches <- str_match(input_vec, "(\\d+x)_(UV[AB]|IR)(\\d*)")
  
  dose   <- matches[, 2]  # "1x", "3x", etc.
  type   <- matches[, 3]  # "UVB", "UVA", "IR"
  energy <- matches[, 4]  # "300", "500", or ""
  
  # Construct final strings
  result <- ifelse(
    energy != "",
    paste0(dose, " ", energy, "J/m² ", type),
    paste0(dose, " ", type)
  )
  return(result)
}

##*********************************************************************************************************
## plot_dose_effect
## This function plots the dose-response effect of overlapping significantly differentially expressed genes 
## across specified irradiation conditions within a given cell type and time point.
## It extracts gene lists using `extract_gene_list()`, finds the common genes across all conditions,
## and visualizes their normalized expression values (variance stabilized) as violin plots with jittered points
## and boxplots. Statistical comparisons between conditions are shown using pairwise tests (wilcoxon test).

## Arguments:
## ls          – named list of DESeq2 results tables containing gene sets for different conditions
## vsd_obj     – variance-stabilized DESeq2 object with expression data and metadata
## conditions  – character vector of irradiation conditions to compare (e.g., c("1x_UVB300", "3x_UVB300", "1x_UVB900"))
## cell        – character string specifying the cell type (e.g., "Melanocytes")
## time        – character string specifying the time point (e.g., "24h")
## point.alpha – numeric controlling the transparency of individual points in the plot (default: 0.1)

## Returns:
## A ggplot object showing violin plots of expression distributions for the intersecting genes 
## across irradiation conditions, with pairwise statistical comparisons.
plot_dose_effect <- function(ls, vsd_obj, conditions, cell, time, point.alpha = 0.1){
  ls <- extract_gene_list(ls = ls, conditions = conditions, cell = cell, time = time)
  
  ## calculate intersecting genes
  overlap <- Reduce(intersect,ls)
  ind <- which(rownames(dds) %in% overlap)
  
  ## convert input into correct format
  cond <- parse_irradiation_input(conditions)
  
  ## change vsd_obj to plot gene counts
  metadata <- colData(vsd_obj)
  vsd_obj <- vsd_obj[ind, metadata$time == time & metadata$cell == cell & metadata$irradiation %in% cond]
  df <-  assay(vsd_obj) %>% 
    t() %>% 
    data.frame() %>%
    rownames_to_column("sample") %>% 
    pivot_longer(cols = -sample, names_to = "gene", values_to = "vsd") %>%
    left_join(data.frame(colData(vsd_obj)), by = "sample") %>% 
    filter(donor != "03_24")
  
  max_y <- max(df$vsd)
  
  ggplot(df, aes(x = irradiation, y = vsd)) +
    geom_violin(
      aes(fill = irradiation), color = "black", width = 0.9, alpha = 0.5, trim = TRUE
    ) +
    geom_jitter(
      position = position_jitter(width = 0.15, height = 0),
      size = 0.9, shape = 21, alpha = point.alpha
    ) +
    geom_boxplot(
      width = 0.1,
      outlier.shape = NA,
      alpha = 0.5,
      color = "black"
    ) +
    scale_fill_manual(values = c("#CD534CFF", "#0073C2FF", "springgreen")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    theme_classic(base_size = 12) +
    labs(
      y = "Normalized expression (VSD)", 
      x = NULL, 
      title = paste0("UVB dose response of \n overlapping genes in ", cell, " (", time, ")")
    ) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold"),
      strip.background = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_text()
    ) +
    stat_compare_means(comparisons = list(c(cond[1], cond[2]),
                                          c(cond[1], cond[3]), 
                                          c(cond[2], cond[3])
    ), label.y = c(max_y*1.05, max_y*1.1, max_y*1.15), paired = FALSE, size = 4, tip.length = 0)
}


##*********************************************************************************************************
## Function to plot control gene expression (DSC vs. Melanocytes) for selected genes with DESeq2 significance
# plot_control_expression_comparison <- function(vsd_obj, results_object, genes, nrow = NULL, ncol = NULL) {
#   
#   ## Extract expression matrix
#   vsd_counts <- assay(vsd_obj)
#   metadata <- colData(vsd_obj)
#   
#   ## Map ENSEMBL IDs to SYMBOLs
#   ensembl_to_symbol <- results_object %>% 
#     filter(SYMBOL %in% genes) %>% 
#     dplyr::select(ENSEMBL, SYMBOL) %>%
#     distinct()
#   
#   if(any(colnames(vsd_counts) != metadata$sample) & all(toupper(colnames(vsd_counts)) == toupper(metadata$ID))){
#     colnames(vsd_counts) <- metadata$sample
#   } else if(any(colnames(vsd_counts) != metadata$sample) & any(toupper(colnames(vsd_counts)) != toupper(metadata$ID))){
#     stop("The column names of the vsd count matrix are not identical to the metadata IDs.")
#   }
#   
#   ## Extract relevant VSD rows
#   gene_rows <- which(rownames(vsd_counts) %in% ensembl_to_symbol$ENSEMBL)
#   df_expr <- t(vsd_counts[gene_rows, ]) %>%
#     as.data.frame() %>%
#     rownames_to_column("sample") 
#   
#   ## Rename columns with SYMBOL
#   colnames(df_expr)[-1] <- ensembl_to_symbol$SYMBOL[match(colnames(df_expr)[-1], ensembl_to_symbol$ENSEMBL)]
#   
#   ## Long format
#   df_long <- df_expr %>%
#     pivot_longer(cols = -sample, names_to = "gene", values_to = "vsd") %>% 
#     left_join(as.data.frame(metadata), by = "sample")
#   
#   ## Add significance annotations
#   sig_annotations <- results_object %>%
#     filter(SYMBOL %in% genes) %>%
#     mutate(label = case_when(
#       is.na(padj)        ~ "",
#       padj < 0.001       ~ "***",
#       padj < 0.01        ~ "**",
#       padj < 0.05        ~ "*",
#       TRUE               ~ ""
#     )) %>%
#     dplyr::select(SYMBOL, label) %>%
#     distinct()
#   
#   df_long <- df_long %>%
#     left_join(sig_annotations, by = c("gene" = "SYMBOL"))
#   
#   ## Plot
#   p <- ggplot(df_long, aes(x = cell, y = vsd, fill = cell)) +
#     geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5, color = "black") +
#     geom_jitter(width = 0.15, size = 1.5, shape = 21, stroke = 0.3, alpha = 0.8, color = "black") +
#     facet_wrap(~ gene, scales = "free_y", nrow = nrow, ncol = ncol) +
#     geom_text(data = distinct(df_long, gene, label),
#               aes(label = label, x = 2, y = Inf),
#               vjust = 1.2, inherit.aes = FALSE, size = 5) +
#     scale_fill_manual(values = c("DSC" = "#CD534CFF", "Melanocytes" = "#0073C2FF")) +
#     ylab("Normalized expression (VSD)") +
#     theme_classic(base_size = 13) +
#     theme(
#       axis.line = element_line(size = 0.6, color = "black"),
#       axis.text = element_text(color = "black"),
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       axis.title.x = element_blank(),
#       legend.position = "bottom",
#       legend.title = element_blank(),
#       legend.key.size = unit(0.6, "cm"),
#       strip.text = element_text(face = "bold", size = 10),
#       panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
#       plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
#     )
#   
#   return(p)
# }

plot_control_expression_comparison <- function(
    vsd.obj, 
    results.object, 
    genes, 
    nrow = NULL, 
    ncol = NULL,
    merge.plots = FALSE
) {
  ## Extract expression matrix
  vsd_counts <- assay(vsd.obj)
  metadata <- colData(vsd.obj)
  
  ## Map ENSEMBL IDs to SYMBOLs
  ensembl_to_symbol <- results.object %>% 
    filter(SYMBOL %in% genes) %>% 
    dplyr::select(ENSEMBL, SYMBOL) %>%
    distinct()
  
  if(any(colnames(vsd_counts) != metadata$sample) & all(toupper(colnames(vsd_counts)) == toupper(metadata$ID))){
    colnames(vsd_counts) <- metadata$sample
  } else if(any(colnames(vsd_counts) != metadata$sample) & any(toupper(colnames(vsd_counts)) != toupper(metadata$ID))){
    stop("The column names of the vsd count matrix are not identical to the metadata IDs.")
  }
  
  ## Extract relevant VSD rows
  gene_rows <- which(rownames(vsd_counts) %in% ensembl_to_symbol$ENSEMBL)
  df_expr <- t(vsd_counts[gene_rows, ]) %>%
    as.data.frame() %>%
    rownames_to_column("sample") 
  
  ## Rename columns with SYMBOL
  colnames(df_expr)[-1] <- ensembl_to_symbol$SYMBOL[match(colnames(df_expr)[-1], ensembl_to_symbol$ENSEMBL)]
  
  ## Long format
  df_long <- df_expr %>%
    pivot_longer(cols = -sample, names_to = "gene", values_to = "vsd") %>% 
    left_join(as.data.frame(metadata), by = "sample")
  
  ## Add significance annotations
  sig_annotations <- results.object %>%
    filter(SYMBOL %in% genes) %>%
    mutate(label = case_when(
      is.na(padj)        ~ "",
      padj < 0.001       ~ "***",
      padj < 0.01        ~ "**",
      padj < 0.05        ~ "*",
      TRUE               ~ ""
    )) %>%
    dplyr::select(SYMBOL, label) %>%
    distinct()
  
  df_long <- df_long %>%
    left_join(sig_annotations, by = c("gene" = "SYMBOL"))
  
  ## Universal ggplot components
  fill_vals <- c("DSC" = "#CD534CFF", "Melanocytes" = "#0073C2FF")
  text_anno <- geom_text(
    data = distinct(df_long, gene, label),
    aes(label = label, x = if(merge.plots) gene else 2, y = Inf),
    vjust = 1.2, inherit.aes = FALSE, size = 5
  )
  
  ## Construct plot depending on 'merge'
  if (!merge.plots) {
    p <- ggplot(df_long, aes(x = cell, y = vsd, fill = cell)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5, color = "black") +
      geom_jitter(width = 0.15, size = 1.5, shape = 21, stroke = 0.3, alpha = 0.8, color = "black") +
      facet_wrap(~ gene, scales = "free_y", nrow = nrow, ncol = ncol)
  } else {
    p <- ggplot(df_long, aes(x = gene, y = vsd, fill = cell)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5, color = "black", position = position_dodge(width = 0.7)) +
      geom_jitter(shape = 21, stroke = 0.3, alpha = 0.8, color = "black", 
                  position = position_dodge(width = 0.7)) +
      geom_vline(xintercept = seq(1.5, length(unique(df_long$gene)) - 0.5, by = 1), linetype = "dashed", color = "gray40") +
      xlab("Gene")
  }
  
  ## Add universal ggplot components
  p <- p +
    text_anno +
    scale_fill_manual(values = fill_vals) +
    ylab("Normalized expression (VSD)") +
    theme_classic(base_size = 13) +
    theme(
      axis.line = element_line(size = 0.6, color = "black"),
      axis.text = element_text(color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.size = unit(0.6, "cm"),
      strip.text = element_text(face = "bold", size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
    )

  return(p)
}

##*********************************************************************************************************
## read edgeR counts for miRNA data
read_edgeR_counts <- function(counts.file){
  counts_edgeR <- t(read.csv(file = counts.file)) 
  colnames(counts_edgeR) <- counts_edgeR[1,]
  counts_edgeR <- counts_edgeR[-1,] 
  names <- rownames(counts_edgeR)
  counts_edgeR <- apply(counts_edgeR, 2, as.numeric)
  rownames(counts_edgeR) <- names 
  as.data.frame(counts_edgeR)
}

