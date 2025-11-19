## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##*********************************************************************************************************
## collection of functions used for the UVDHDS_transcriptomics project

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
## plot_control_expression_comparison
## This function generates boxplots of normalized gene expression values (DESeq2 VSD) 
## for a given list of genes in control (non-irradiated) samples, comparing DSCs and Melanocytes.
## Optionally, the plots can be merged to show both cell types for all genes in a single panel, or faceted by gene.
## Significance is indicated with asterisks above the plots, based on DESeq2 padj values.

## Arguments:
## vsd.obj        – DESeq2 VSD (variance-stabilized transformation) object containing normalized counts
## results.object – data frame or DESeq2 results table containing gene annotations (SYMBOL, ENSEMBL, padj columns)
## genes          – character vector of gene symbols to plot
## nrow           – number of rows in the facet plot (only if merge.plots = FALSE)
## ncol           – number of columns in the facet plot (only if merge.plots = FALSE)
## merge.plots    – logical; if TRUE, boxplots for all genes are shown in one panel with cell types grouped, 
##                  if FALSE, each gene gets a separate facet

## Returns:
## A ggplot object showing boxplots of normalized expression for selected genes in DSCs and Melanocytes,
## with significance annotations ("*", "**", "***") above the plots (based on DESeq2 results).
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

