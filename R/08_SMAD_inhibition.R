## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<HEAD>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

##*********************************************************************************************************

## load packages
library(tidyverse)
library(patchwork)
library(ggprism)
library(rstatix)
library(readxl)
library(ekbSeq)

##*********************************************************************************************************
## Load data

dat <- read_excel("Data/SMAD_inhibition_treatment.xlsx") %>%
  mutate(Gene = str_replace_all(.$Gene, "Nestin", "NES"))

#------------------------------------------------------------
# 1. Collapse technical replicates within biological replicate
#------------------------------------------------------------

dat_bio <- dat %>%
  group_by(Donor, Time, Gene) %>%
  summarise(
    ddCT = mean(ddCT, na.rm = TRUE),
    .groups = "drop"
  ) 

#------------------------------------------------------------
# 2. Summarise biological replicates: mean ± SD
#------------------------------------------------------------

dat_sum <- dat_bio %>%
  group_by(Time, Gene) %>%
  summarise(
    mean_ddCT = mean(ddCT, na.rm = TRUE),
    sd_ddCT   = sd(ddCT, na.rm = TRUE),
    n         = n(),
    .groups = "drop"
  )

#------------------------------------------------------------
# 2. One-sample t-test vs 0, adjusted within each gene
#------------------------------------------------------------

stat_ddct <- dat_bio %>%
  group_by(Gene, Time) %>%
  t_test(ddCT ~ 1, mu = 0) %>%
  ungroup() %>%
  group_by(Gene) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  ungroup()

#------------------------------------------------------------
# 3. Prepare significance labels
#------------------------------------------------------------

sig_labels <- dat_sum %>%
  left_join(
    stat_ddct %>%
      select(Gene, Time, p, p.adj, p.adj.signif),
    by = c("Gene", "Time")
  ) %>%
  group_by(Gene) %>%
  mutate(
    y_range = max(mean_ddCT + sd_ddCT, na.rm = TRUE) -
      min(mean_ddCT - sd_ddCT, na.rm = TRUE),
    y_offset = if_else(y_range == 0, 0.1, 0.08 * y_range),
    y.position = mean_ddCT + sd_ddCT + y_offset
  ) %>%
  ungroup() %>%
  filter(!is.na(p.adj.signif), p.adj.signif != "ns")

#------------------------------------------------------------
# 3. Plotting function
#------------------------------------------------------------

plot_gene_ddct <- function(gene_name,
                           data_sum,
                           data_bio,
                           sig_data,
                           show_y_title = FALSE) {
  
  plot_sum <- data_sum %>%
    filter(Gene == gene_name)
  
  plot_bio <- data_bio %>%
    filter(Gene == gene_name)
  
  plot_sig <- sig_data %>%
    filter(Gene == gene_name)
  
  p <- ggplot(plot_sum, aes(x = Time, y = mean_ddCT)) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      linetype = "dashed",
      color = "grey60"
    )
  
  p <- p +
    geom_line(
      aes(group = 1),
      linewidth = 0.8, 
      color = "black"
    ) +
    geom_errorbar(
      aes(
        ymin = mean_ddCT - sd_ddCT,
        ymax = mean_ddCT + sd_ddCT
      ),
      width = 3,
      linewidth = 0.6,
      color = "black"
    ) +
    geom_point(
      size = 2.4,
      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.8
    )
  
  if (nrow(plot_sig) > 0) {
    p <- p +
      geom_text(
        data = plot_sig,
        aes(
          x = Time,
          y = y.position,
          label = p.adj.signif
        ),
        inherit.aes = FALSE,
        size = 5,
        color = "black",
        vjust = 0
      )
  }
  
  p <- p +
    labs(
      title = gene_name,
      x = "Time (h)",
      y = if (show_y_title) expression(Delta * Delta * CT) else NULL
    ) +
    theme_prism(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text = element_text(size = 9),
      panel.grid.major.y = element_line(linewidth = 0.25, color = "grey88"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(5.5, 5.5, 5.5, 5.5)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.10, 0.25))
    ) +
    scale_x_continuous(
      limits = c(0, 75), 
      breaks = c(0, 24, 48, 72)
    )
  return(p)
}

#------------------------------------------------------------
# 4. Plot results
#------------------------------------------------------------

genes_to_plot <- c("COL1A1", "SMAD7", "SOX2", "NES", "MITF")

plots <- map2(
  genes_to_plot,
  seq_along(genes_to_plot),
  ~ plot_gene_ddct(
    gene_name = .x,
    data_sum = dat_sum,
    data_bio = dat_bio,
    sig_data = sig_labels,
    show_y_title = .y == 1
  )
)

combined_plot <- wrap_plots(plots, nrow = 1) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 12)
  )

export_plot_dual("Results/smad_inhibition_results", combined_plot, width = 15, height = 2.5)



