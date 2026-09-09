rm(list = ls())

########## Load required packages ##########
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)

########## Read data ##########
coi_data <- read.csv(
  "COI_data.csv",
  check.names = FALSE
)

rrna_data <- read.csv(
  "18S_data.csv",
  check.names = FALSE
)

metadata <- read.csv(
  "metadata.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

########## Factor levels and plot settings ##########
marker_levels <- c(
  "COI",
  "18S rRNA"
)

layer_levels <- c(
  "Topsoil",
  "Subsoil"
)

elevation_levels <- c(
  "300",
  "600",
  "900",
  "1200",
  "1600",
  "1800",
  "2100"
)

marker_offset <- c(
  "COI" = -0.20,
  "18S rRNA" = 0.20
)

fill_colours <- c(
  "COI" = "#F87A52",
  "18S rRNA" = "#188CBA"
)

point_colours <- c(
  "COI" = "#F64007",
  "18S rRNA" = "#188CBA"
)

tax_cols <- c(
  "ASV",
  "Domain",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

########## Calculate observed ASV richness ##########
calculate_richness <- function(
    data,
    marker_name
) {

  sample_cols <- setdiff(
    colnames(data),
    tax_cols
  )

  abundance <- as.matrix(
    data[
      ,
      sample_cols,
      drop = FALSE
    ]
  )

  storage.mode(abundance) <- "numeric"

  tibble(
    Sample = sample_cols,
    Richness = colSums(
      abundance > 0,
      na.rm = TRUE
    ),
    Marker = marker_name
  )
}

richness_data <- bind_rows(
  calculate_richness(
    coi_data,
    "COI"
  ),
  calculate_richness(
    rrna_data,
    "18S rRNA"
  )
)

metadata_clean <- metadata %>%
  transmute(
    Sample = as.character(`OTU ID`),
    Group = factor(
      as.character(Group),
      levels = c(
        "T",
        "S"
      ),
      labels = layer_levels
    ),
    Elevation = factor(
      as.character(Elevation),
      levels = elevation_levels
    ),
    PlotID = sub(
      "^[TtSs][_-]?",
      "",
      as.character(`OTU ID`)
    )
  )

richness_data <- richness_data %>%
  left_join(
    metadata_clean,
    by = "Sample"
  ) %>%
  mutate(
    Marker = factor(
      Marker,
      levels = marker_levels
    ),
    PlotID = factor(PlotID)
  ) %>%
  filter(
    !is.na(Richness),
    !is.na(Group),
    !is.na(Elevation),
    !is.na(PlotID)
  ) %>%
  droplevels()

cat(
  "\nASV richness data summary:\n",
  "Rows: ",
  nrow(richness_data),
  "\nSamples: ",
  n_distinct(richness_data$Sample),
  "\nPlots: ",
  n_distinct(richness_data$PlotID),
  "\n",
  sep = ""
)

########## Helper functions ##########
custom_signif <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

paired_wilcox_result <- function(
    data,
    column_1,
    column_2
) {

  complete_data <- data %>%
    filter(
      !is.na(.data[[column_1]]),
      !is.na(.data[[column_2]])
    )

  if (nrow(complete_data) < 2) {
    return(
      tibble(
        n_pairs = nrow(complete_data),
        statistic = NA_real_,
        p = NA_real_
      )
    )
  }

  test <- wilcox.test(
    complete_data[[column_1]],
    complete_data[[column_2]],
    paired = TRUE,
    exact = FALSE
  )

  tibble(
    n_pairs = nrow(complete_data),
    statistic = unname(
      test$statistic
    ),
    p = test$p.value
  )
}

y_max <- max(
  richness_data$Richness,
  na.rm = TRUE
)

y_min <- min(
  richness_data$Richness,
  na.rm = TRUE
)

y_range <- y_max - y_min

if (y_range == 0) {
  y_range <- max(
    abs(y_max),
    1
  )
}

y_step <- y_range * 0.10

########## Fig. S5A: overall ASV richness ##########
overall_wide <- richness_data %>%
  select(
    Sample,
    Marker,
    Richness
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = Marker,
    values_from = Richness
  )

overall_stats <- paired_wilcox_result(
  overall_wide,
  "COI",
  "18S rRNA"
) %>%
  mutate(
    comparison = "Overall: COI vs 18S rRNA",
    group1 = "COI",
    group2 = "18S rRNA",
    p.adj = p,
    signif = custom_signif(p.adj),
    xmin = 1,
    xmax = 2,
    y.position = y_max + y_step
  )

p_overall <- ggplot(
  richness_data,
  aes(
    x = Marker,
    y = Richness,
    fill = Marker
  )
) +
  geom_boxplot(
    alpha = 0.9,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_point(
    aes(
      color = Marker
    ),
    position = position_jitter(
      width = 0.2
    ),
    alpha = 0.9,
    size = 2.25
  )

if (
  !is.na(overall_stats$p.adj) &&
  overall_stats$p.adj < 0.05
) {
  p_overall <- p_overall +
    stat_pvalue_manual(
      overall_stats,
      label = "signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

p_overall <- p_overall +
  scale_fill_manual(
    values = fill_colours
  ) +
  scale_color_manual(
    values = point_colours
  ) +
  scale_x_discrete(
    limits = marker_levels
  ) +
  labs(
    x = "Genetic marker",
    y = "ASV richness"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 11,
      color = "black"
    )
  )

########## Fig. S5B: ASV richness by soil layer ##########
marker_layer_wide <- richness_data %>%
  select(
    Group,
    Sample,
    Marker,
    Richness
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = Marker,
    values_from = Richness
  )

marker_layer_stats <- marker_layer_wide %>%
  group_by(Group) %>%
  group_modify(
    ~ paired_wilcox_result(
      .x,
      "COI",
      "18S rRNA"
    )
  ) %>%
  ungroup() %>%
  mutate(
    group1 = "COI",
    group2 = "18S rRNA",
    p.adj = p.adjust(
      p,
      method = "BH"
    ),
    signif = custom_signif(p.adj),
    group_index = match(
      as.character(Group),
      layer_levels
    ),
    xmin = group_index +
      marker_offset["COI"],
    xmax = group_index +
      marker_offset["18S rRNA"],
    y.position = y_max + y_step
  )

layer_marker_wide <- richness_data %>%
  select(
    PlotID,
    Marker,
    Group,
    Richness
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = Group,
    values_from = Richness
  )

layer_marker_stats <- layer_marker_wide %>%
  group_by(Marker) %>%
  group_modify(
    ~ paired_wilcox_result(
      .x,
      "Topsoil",
      "Subsoil"
    )
  ) %>%
  ungroup() %>%
  mutate(
    group1 = "Topsoil",
    group2 = "Subsoil",
    p.adj = p.adjust(
      p,
      method = "BH"
    ),
    signif = custom_signif(p.adj),
    xmin = 1 +
      marker_offset[
        as.character(Marker)
      ],
    xmax = 2 +
      marker_offset[
        as.character(Marker)
      ],
    y.position =
      y_max +
      y_step * 2 +
      ifelse(
        Marker == "COI",
        0,
        y_step
      )
  )

p_layer <- ggplot(
  richness_data,
  aes(
    x = Group,
    y = Richness,
    fill = Marker
  )
) +
  geom_boxplot(
    alpha = 0.9,
    outlier.shape = NA,
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  geom_point(
    aes(
      color = Marker
    ),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.8
    ),
    alpha = 0.9,
    size = 2.25
  )

marker_layer_sig <- marker_layer_stats %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

layer_marker_sig <- layer_marker_stats %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

if (nrow(marker_layer_sig) > 0) {
  p_layer <- p_layer +
    stat_pvalue_manual(
      marker_layer_sig,
      label = "signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

if (nrow(layer_marker_sig) > 0) {
  p_layer <- p_layer +
    stat_pvalue_manual(
      layer_marker_sig,
      label = "signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

p_layer <- p_layer +
  scale_fill_manual(
    values = fill_colours
  ) +
  scale_color_manual(
    values = point_colours
  ) +
  labs(
    x = "Soil layer",
    y = "ASV richness"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 11,
      color = "black"
    ),
    legend.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      color = "black"
    )
  )

########## Fig. S5C: ASV richness across elevations ##########
elevation_index <- setNames(
  seq_along(
    elevation_levels
  ),
  elevation_levels
)

marker_elevation_wide <- richness_data %>%
  select(
    Elevation,
    Sample,
    Marker,
    Richness
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = Marker,
    values_from = Richness
  )

marker_elevation_stats <- marker_elevation_wide %>%
  group_by(Elevation) %>%
  group_modify(
    ~ paired_wilcox_result(
      .x,
      "COI",
      "18S rRNA"
    )
  ) %>%
  ungroup() %>%
  mutate(
    group1 = "COI",
    group2 = "18S rRNA",
    p.adj = p.adjust(
      p,
      method = "BH"
    ),
    signif = custom_signif(p.adj),
    x_center = elevation_index[
      as.character(Elevation)
    ],
    xmin = x_center +
      marker_offset["COI"],
    xmax = x_center +
      marker_offset["18S rRNA"],
    y.position = y_max + y_step
  )

kruskal_results <- list()
dunn_results <- list()

for (marker_name in marker_levels) {

  marker_data <- richness_data %>%
    filter(
      Marker == marker_name
    ) %>%
    droplevels()

  kw <- rstatix::kruskal_test(
    marker_data,
    Richness ~ Elevation
  ) %>%
    mutate(
      Marker = marker_name
    )

  kruskal_results[[marker_name]] <- kw

  if (
    !is.na(kw$p[1]) &&
    kw$p[1] < 0.05
  ) {

    dunn_results[[marker_name]] <-
      rstatix::dunn_test(
        marker_data,
        Richness ~ Elevation,
        p.adjust.method = "bonferroni"
      ) %>%
      filter(
        p.adj < 0.05
      ) %>%
      mutate(
        Marker = marker_name,
        signif = custom_signif(p.adj),
        xmin =
          elevation_index[
            as.character(group1)
          ] +
          marker_offset[marker_name],
        xmax =
          elevation_index[
            as.character(group2)
          ] +
          marker_offset[marker_name]
      )
  }
}

kruskal_results <- bind_rows(
  kruskal_results
)

dunn_results <- bind_rows(
  dunn_results
)

if (nrow(dunn_results) > 0) {
  dunn_results <- dunn_results %>%
    arrange(
      Marker,
      xmin,
      xmax
    ) %>%
    mutate(
      line_number = row_number(),
      y.position =
        y_max +
        y_step * 2 +
        y_step * line_number
    )
}

p_elevation <- ggplot(
  richness_data,
  aes(
    x = Elevation,
    y = Richness,
    fill = Marker
  )
) +
  geom_boxplot(
    alpha = 0.9,
    outlier.shape = NA,
    position = position_dodge(
      width = 0.8
    ),
    width = 0.7
  ) +
  geom_point(
    aes(
      color = Marker
    ),
    position = position_jitterdodge(
      jitter.width = 0.2,
      dodge.width = 0.8
    ),
    alpha = 0.9,
    size = 2.25
  )

marker_elevation_sig <- marker_elevation_stats %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

if (nrow(marker_elevation_sig) > 0) {
  p_elevation <- p_elevation +
    stat_pvalue_manual(
      marker_elevation_sig,
      label = "signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

if (nrow(dunn_results) > 0) {
  p_elevation <- p_elevation +
    stat_pvalue_manual(
      dunn_results,
      label = "signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

p_elevation <- p_elevation +
  scale_fill_manual(
    values = fill_colours
  ) +
  scale_color_manual(
    values = point_colours
  ) +
  labs(
    x = "Elevation (m)",
    y = "ASV richness"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text = element_text(
      size = 10,
      color = "black"
    ),
    axis.title = element_text(
      size = 11,
      color = "black"
    ),
    legend.text = element_text(
      color = "black"
    ),
    legend.title = element_text(
      color = "black"
    )
  )

########## Show statistical results ##########
cat(
  "\nOverall marker comparison:\n"
)
print(
  overall_stats
)

cat(
  "\nMarker comparisons within soil layers:\n"
)
print(
  marker_layer_stats
)

cat(
  "\nSoil-layer comparisons within markers:\n"
)
print(
  layer_marker_stats
)

cat(
  "\nMarker comparisons within elevations:\n"
)
print(
  marker_elevation_stats
)

cat(
  "\nKruskal-Wallis tests among elevations:\n"
)
print(
  kruskal_results
)

if (nrow(dunn_results) > 0) {
  cat(
    "\nSignificant Dunn-Bonferroni comparisons:\n"
  )
  print(
    dunn_results
  )
}

########## Display Fig. S5A-C ##########
print(
  p_overall
)

print(
  p_layer
)

print(
  p_elevation
)
