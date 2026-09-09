rm(list = ls())

########## Load required packages ##########
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)

########## Read data ##########
coi_data <- read.csv("COI_data.csv", check.names = FALSE, stringsAsFactors = FALSE)
rrna_data <- read.csv("18S_data.csv", check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv("metadata.csv", check.names = FALSE, stringsAsFactors = FALSE)

########## Analysis settings ##########
nperm <- 9999
analysis_seed <- 20260721
elevation_levels <- c("300", "600", "900", "1200", "1600", "1800", "2100")

metadata <- metadata %>%
  transmute(
    Sample = as.character(`OTU ID`),
    Group = factor(Group, levels = c("T", "S"), labels = c("Topsoil", "Subsoil")),
    Elevation = factor(as.character(Elevation), levels = elevation_levels),
    PlotID = factor(paste0("P", sub("^[TtSs][_-]?", "", as.character(`OTU ID`))))
  )

########## Prepare marker-specific complete pairs ##########
prepare_marker <- function(data, metadata, marker_name) {
  meta <- metadata %>%
    filter(Sample %in% intersect(metadata$Sample, colnames(data))) %>%
    group_by(PlotID) %>%
    filter(n() == 2, n_distinct(Group) == 2, n_distinct(Elevation) == 1) %>%
    ungroup() %>%
    arrange(PlotID, Group) %>%
    droplevels()

  community <- t(as.matrix(data[, meta$Sample, drop = FALSE]))
  storage.mode(community) <- "numeric"
  colnames(community) <- as.character(data$ASV)
  community <- community[, colSums(community) > 0, drop = FALSE]

  taxonomy <- data %>%
    transmute(
      ASV_ID = as.character(ASV),
      Genus = trimws(sub("^[a-z]__", "", as.character(Genus)))
    )

  cat(
    "\n", marker_name,
    "\nComplete plots: ", n_distinct(meta$PlotID),
    "\nSamples: ", nrow(meta), "\n",
    sep = ""
  )

  list(community = community, metadata = meta, taxonomy = taxonomy)
}

########## Restricted permutation matrices ##########
make_whole_plot_permutations <- function(metadata, nset, seed) {
  set.seed(seed)
  plot_ids <- unique(as.character(metadata$PlotID))
  layer_ids <- levels(metadata$Group)
  permutations <- matrix(NA_integer_, nrow = nset, ncol = nrow(metadata))

  for (b in seq_len(nset)) {
    source_plots <- sample(plot_ids)
    perm_b <- integer(nrow(metadata))

    for (j in seq_along(plot_ids)) {
      for (layer in layer_ids) {
        destination <- which(
          as.character(metadata$PlotID) == plot_ids[j] &
            as.character(metadata$Group) == layer
        )
        source <- which(
          as.character(metadata$PlotID) == source_plots[j] &
            as.character(metadata$Group) == layer
        )
        perm_b[destination] <- source
      }
    }

    permutations[b, ] <- perm_b
  }

  permutations
}

make_within_plot_permutations <- function(metadata, nset, seed) {
  set.seed(seed)
  plot_ids <- unique(as.character(metadata$PlotID))
  permutations <- matrix(NA_integer_, nrow = nset, ncol = nrow(metadata))

  for (b in seq_len(nset)) {
    perm_b <- seq_len(nrow(metadata))

    for (plot_id in plot_ids) {
      rows <- which(as.character(metadata$PlotID) == plot_id)
      perm_b[rows] <- sample(rows, 2)
    }

    permutations[b, ] <- perm_b
  }

  permutations
}

########## PERMANOVA helper ##########
extract_permanova_term <- function(result, terms, marker_name) {
  result_df <- as.data.frame(result, check.names = FALSE)
  result_df$Term <- rownames(result_df)
  selected <- result_df %>% filter(Term %in% terms)

  tibble(
    Marker = marker_name,
    Term = selected$Term,
    Df = selected$Df,
    Pseudo_F = selected$F,
    R2 = selected$R2,
    P_value = selected$`Pr(>F)`
  )
}

########## PCoA, SIMPER, and PERMANOVA ##########
analyse_marker <- function(dataset, marker_name, seed_offset = 0) {
  community <- dataset$community
  meta <- dataset$metadata

  # Same Bray-Curtis distance for PCoA and PERMANOVA
  bray <- vegdist(community, method = "bray")

  whole_perm <- make_whole_plot_permutations(
    meta, nperm, analysis_seed + seed_offset
  )
  within_perm <- make_within_plot_permutations(
    meta, nperm, analysis_seed + seed_offset + 1
  )

  permanova_whole <- adonis2(
    bray ~ Elevation * Group,
    data = meta,
    permutations = whole_perm,
    by = "terms"
  )
  permanova_within <- adonis2(
    bray ~ Elevation * Group,
    data = meta,
    permutations = within_perm,
    by = "terms"
  )

  permanova <- bind_rows(
    extract_permanova_term(permanova_whole, "Elevation", marker_name),
    extract_permanova_term(permanova_within, "Group", marker_name),
    extract_permanova_term(
      permanova_within,
      c("Elevation:Group", "Group:Elevation"),
      marker_name
    )
  ) %>%
    mutate(
      Term = case_when(
        Term == "Group" ~ "Soil layer",
        Term %in% c("Elevation:Group", "Group:Elevation") ~ "Elevation × Soil layer",
        TRUE ~ Term
      ),
      Significance = case_when(
        P_value < 0.001 ~ "***",
        P_value < 0.01 ~ "**",
        P_value < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )

  # PCoA
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  positive_eigenvalues <- pcoa$eig[pcoa$eig > 0]
  variance_percent <- 100 * pcoa$eig[1:2] / sum(positive_eigenvalues)

  pcoa_scores <- as.data.frame(pcoa$points[, 1:2, drop = FALSE])
  colnames(pcoa_scores) <- c("PCo1", "PCo2")
  pcoa_scores$Sample <- rownames(pcoa_scores)

  pcoa_scores <- pcoa_scores %>%
    left_join(meta %>% select(Sample, PlotID, Elevation, Group), by = "Sample") %>%
    mutate(
      SoilLayer = factor(
        Group,
        levels = c("Topsoil", "Subsoil"),
        labels = c("0-10 cm", "10-20 cm")
      )
    )

  # SIMPER using all retained ASVs
  simper_result <- simper(community, group = meta$Group, permutations = 999)[[1]]
  asv_ids <- rownames(simper_result)
  if (is.null(asv_ids)) asv_ids <- names(simper_result$average)

  simper_top5 <- tibble(
    ASV_ID = as.character(asv_ids),
    Average_contribution = as.numeric(simper_result$average)
  ) %>%
    mutate(
      Contribution_percent = Average_contribution / sum(Average_contribution) * 100
    ) %>%
    left_join(dataset$taxonomy, by = "ASV_ID") %>%
    mutate(
      Genus = case_when(
        is.na(Genus) ~ "Unknown",
        trimws(Genus) == "" ~ "Unknown",
        tolower(trimws(Genus)) %in% c(
          "unidentified", "unclassified", "uncultured", "unknown", "na"
        ) ~ "Unknown",
        TRUE ~ Genus
      )
    ) %>%
    group_by(Genus) %>%
    summarise(TotalContribution = sum(Contribution_percent), .groups = "drop") %>%
    arrange(desc(TotalContribution)) %>%
    slice_head(n = 5) %>%
    mutate(DisplayName = Genus)

  list(
    permanova = permanova,
    pcoa_scores = pcoa_scores,
    pcoa_variance = variance_percent,
    simper = simper_top5
  )
}

########## Run analyses ##########
coi_dataset <- prepare_marker(coi_data, metadata, "COI")
rrna_dataset <- prepare_marker(rrna_data, metadata, "18S rRNA")

coi_results <- analyse_marker(coi_dataset, "COI", seed_offset = 0)
rrna_results <- analyse_marker(rrna_dataset, "18S rRNA", seed_offset = 100)

combined_permanova <- bind_rows(coi_results$permanova, rrna_results$permanova)

########## Show results ##########
cat("\nPERMANOVA results:\n")
print(combined_permanova)

cat("\nPCoA explained variance (%):\n")
print(
  tibble(
    Marker = c("COI", "18S rRNA"),
    PCo1 = c(coi_results$pcoa_variance[1], rrna_results$pcoa_variance[1]),
    PCo2 = c(coi_results$pcoa_variance[2], rrna_results$pcoa_variance[2])
  )
)

cat("\nCOI SIMPER top five genera:\n")
print(coi_results$simper)

cat("\n18S rRNA SIMPER top five genera:\n")
print(rrna_results$simper)

########## PCoA plot ##########
plot_pcoa <- function(scores, variance_percent, permanova_table, marker_name) {
  layer_colors <- c("0-10 cm" = "#fd79cd", "10-20 cm" = "#875540")

  annotation <- permanova_table %>%
    filter(Marker == marker_name, Term == "Soil layer")

  p_text <- ifelse(
    annotation$P_value < 0.001,
    "P < 0.001",
    paste0("P = ", format(round(annotation$P_value, 3), nsmall = 3))
  )

  label_text <- paste0(
    "PERMANOVA\nR² = ", round(annotation$R2, 3), "\n", p_text
  )

  x_range <- range(scores$PCo1)
  y_range <- range(scores$PCo2)

  ggplot(
    scores,
    aes(x = PCo1, y = PCo2, color = SoilLayer, fill = SoilLayer)
  ) +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(
      geom = "polygon",
      alpha = 0.2,
      level = 0.95,
      show.legend = FALSE
    ) +
    scale_color_manual(values = layer_colors, name = "Soil layer") +
    scale_fill_manual(values = layer_colors, name = "Soil layer") +
    labs(
      x = paste0("PCo1 (", round(variance_percent[1], 2), "%)"),
      y = paste0("PCo2 (", round(variance_percent[2], 2), "%)")
    ) +
    annotate(
      "text",
      x = x_range[1] + 0.02 * diff(x_range),
      y = y_range[2] - 0.05 * diff(y_range),
      label = label_text,
      size = 4,
      hjust = 0,
      vjust = 1,
      color = "black",
      fontface = "bold"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
}

########## SIMPER plot ##########
plot_simper <- function(data, marker_name) {
  fill_col <- ifelse(marker_name == "COI", "#F64007", "#5EAFCF")

  plot_data <- data %>%
    arrange(TotalContribution) %>%
    mutate(DisplayName = factor(DisplayName, levels = DisplayName))

  y_max <- ceiling(max(plot_data$TotalContribution) / 5) * 5

  ggplot(plot_data, aes(x = DisplayName, y = TotalContribution)) +
    geom_col(fill = fill_col, alpha = 0.8) +
    coord_flip() +
    labs(x = NULL, y = "Contribution to dissimilarity (%)") +
    scale_y_continuous(limits = c(0, y_max), expand = c(0, 0)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black")
    )
}

########## PERMANOVA plot ##########
effect_colors <- c(
  "Elevation" = "#27B3B8",
  "Soil layer" = "#F17870",
  "Elevation × Soil layer" = "#717A83"
)

plot_permanova <- function(data, marker_name) {
  plot_data <- data %>%
    filter(Marker == marker_name) %>%
    mutate(
      Effect = factor(
        Term,
        levels = c("Elevation × Soil layer", "Soil layer", "Elevation")
      ),
      Variance_explained = R2 * 100,
      Percent_label = paste0(format(round(Variance_explained, 1), nsmall = 1), "%"),
      Percent_x = pmax(Variance_explained - 0.8, Variance_explained * 0.55),
      Star_x = Variance_explained + 0.6,
      Plot_significance = ifelse(P_value < 0.05, Significance, "")
    )

  ggplot(
    plot_data,
    aes(x = Variance_explained, y = Effect, fill = Effect)
  ) +
    geom_col(width = 0.72) +
    geom_text(
      aes(x = Percent_x, label = Percent_label),
      hjust = 1,
      size = 4.8,
      fontface = "bold",
      family = "serif",
      color = "black"
    ) +
    geom_text(
      data = plot_data %>% filter(Plot_significance != ""),
      aes(x = Star_x, label = Plot_significance),
      hjust = 0,
      size = 6,
      fontface = "bold",
      family = "serif",
      color = "black"
    ) +
    scale_fill_manual(values = effect_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, 30),
      breaks = c(0, 10, 20, 30),
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(title = marker_name, x = "Variance explained (%)", y = NULL) +
    theme_bw(base_size = 14, base_family = "serif") +
    theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.title.x = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 13, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = margin(t = 8, r = 15, b = 8, l = 8)
    )
}

########## Generate Figure 3A and 3B ##########
p_coi_pcoa <- plot_pcoa(
  coi_results$pcoa_scores,
  coi_results$pcoa_variance,
  combined_permanova,
  "COI"
)
p_rrna_pcoa <- plot_pcoa(
  rrna_results$pcoa_scores,
  rrna_results$pcoa_variance,
  combined_permanova,
  "18S rRNA"
)

p_coi_simper <- plot_simper(coi_results$simper, "COI")
p_rrna_simper <- plot_simper(rrna_results$simper, "18S rRNA")

p_coi_permanova <- plot_permanova(combined_permanova, "COI")
p_rrna_permanova <- plot_permanova(combined_permanova, "18S rRNA")

# Relative widths follow the original plot dimensions.
p_fig3a <- p_coi_pcoa + p_coi_simper + p_coi_permanova +
  plot_layout(ncol = 3, widths = c(6.3, 3.45, 6))

p_fig3b <- p_rrna_pcoa + p_rrna_simper + p_rrna_permanova +
  plot_layout(ncol = 3, widths = c(6.3, 3.45, 6))

print(p_fig3a)
print(p_fig3b)
