rm(list = ls())

########## Load required packages ##########
library(vegan)
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

########## Set factor order ##########
marker_levels <- c(
  "COI",
  "18S rRNA"
)

elevation_levels <- c(
  "300", "600", "900", "1200",
  "1600", "1800", "2100"
)

layer_levels <- c(
  "Topsoil",
  "Subsoil"
)

marker_offset <- c(
  "COI" = -0.20,
  "18S rRNA" = 0.20
)

########## Prepare metadata ##########
metadata <- metadata %>%
  transmute(
    Sample = as.character(`OTU ID`),
    Group = factor(
      Group,
      levels = c("T", "S"),
      labels = layer_levels
    ),
    Elevation = factor(
      as.character(Elevation),
      levels = elevation_levels
    ),
    PlotID = paste0(
      "P",
      sub(
        "^[TtSs][_-]?",
        "",
        as.character(`OTU ID`)
      )
    )
  )

########## Calculate Shannon diversity ##########
calculate_shannon <- function(data, metadata, marker_name) {
  
  sample_cols <- intersect(
    metadata$Sample,
    colnames(data)
  )
  
  abundance <- as.matrix(
    data[, sample_cols]
  )
  
  tibble(
    Sample = sample_cols,
    Shannon = diversity(
      t(abundance),
      index = "shannon"
    ),
    Marker = marker_name
  )
}

coi_alpha <- calculate_shannon(
  coi_data,
  metadata,
  "COI"
)

rrna_alpha <- calculate_shannon(
  rrna_data,
  metadata,
  "18S rRNA"
)

alpha_data <- bind_rows(
  coi_alpha,
  rrna_alpha
) %>%
  left_join(
    metadata,
    by = "Sample"
  ) %>%
  mutate(
    Marker = factor(
      Marker,
      levels = marker_levels
    ),
    Group = factor(
      Group,
      levels = layer_levels
    ),
    Elevation = factor(
      Elevation,
      levels = elevation_levels
    ),
    PlotID = factor(PlotID)
  )

########## Helper functions ##########
get_significance <- function(p) {
  case_when(
    is.na(p) ~ "ns",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

paired_wilcox <- function(data, var1, var2) {
  
  complete_data <- data %>%
    filter(
      !is.na(.data[[var1]]),
      !is.na(.data[[var2]])
    )
  
  test <- wilcox.test(
    complete_data[[var1]],
    complete_data[[var2]],
    paired = TRUE,
    exact = FALSE
  )
  
  tibble(
    n_pairs = nrow(complete_data),
    statistic = unname(test$statistic),
    p = test$p.value
  )
}

get_y_parameters <- function(data) {
  
  y_max <- max(
    data$Shannon,
    na.rm = TRUE
  )
  
  y_min <- min(
    data$Shannon,
    na.rm = TRUE
  )
  
  y_range <- y_max - y_min
  
  if (
    !is.finite(y_range) ||
    y_range == 0
  ) {
    y_range <- max(
      abs(y_max),
      1
    )
  }
  
  list(
    y_max = y_max,
    step = y_range * 0.10
  )
}

########## Figure 2C: overall marker comparison ##########
overall_wide <- alpha_data %>%
  select(
    Sample,
    Marker,
    Shannon
  ) %>%
  pivot_wider(
    names_from = Marker,
    values_from = Shannon
  )

stats_2c <- paired_wilcox(
  overall_wide,
  "COI",
  "18S rRNA"
) %>%
  mutate(
    group1 = "COI",
    group2 = "18S rRNA",
    p.adj = p,
    significance = get_significance(p.adj)
  )

print(stats_2c)

y_2c <- get_y_parameters(
  alpha_data
)

annotation_2c <- stats_2c %>%
  mutate(
    xmin = 1,
    xmax = 2,
    y.position =
      y_2c$y_max +
      y_2c$step
  )

p2c <- ggplot(
  alpha_data,
  aes(
    x = Marker,
    y = Shannon,
    fill = Marker
  )
) +
  geom_boxplot(
    alpha = 0.9,
    outlier.shape = NA,
    width = 0.7
  ) +
  geom_point(
    aes(color = Marker),
    position = position_jitter(
      width = 0.2
    ),
    alpha = 0.9,
    size = 2.25
  )

if (
  !is.na(stats_2c$p.adj) &&
  stats_2c$p.adj < 0.05
) {
  p2c <- p2c +
    stat_pvalue_manual(
      annotation_2c,
      label = "significance",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

p2c <- p2c +
  labs(
    x = "Genetic marker",
    y = "Shannon diversity"
  ) +
  scale_fill_manual(
    values = c(
      "COI" = "#F87A52",
      "18S rRNA" = "#188CBA"
    )
  ) +
  scale_color_manual(
    values = c(
      "COI" = "#F64007",
      "18S rRNA" = "#188CBA"
    )
  ) +
  scale_x_discrete(
    limits = marker_levels
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

print(p2c)

########## Figure 2D: elevation comparison ##########
elevation_wide <- alpha_data %>%
  select(
    Elevation,
    Sample,
    Marker,
    Shannon
  ) %>%
  pivot_wider(
    names_from = Marker,
    values_from = Shannon
  )

stats_marker_elevation <- elevation_wide %>%
  group_by(Elevation) %>%
  group_modify(
    ~ paired_wilcox(
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
    significance = get_significance(
      p.adj
    )
  )

kw_results <- list()
dunn_results <- list()

for (marker in marker_levels) {
  
  marker_data <- alpha_data %>%
    filter(
      Marker == marker
    ) %>%
    droplevels()
  
  kw <- kruskal_test(
    marker_data,
    Shannon ~ Elevation
  ) %>%
    mutate(
      Marker = marker
    )
  
  kw_results[[marker]] <- kw
  
  if (
    !is.na(kw$p[1]) &&
    kw$p[1] < 0.05
  ) {
    
    dunn_results[[marker]] <-
      dunn_test(
        marker_data,
        Shannon ~ Elevation,
        p.adjust.method = "bonferroni"
      ) %>%
      filter(
        p.adj < 0.05
      ) %>%
      mutate(
        Marker = marker,
        significance =
          get_significance(p.adj)
      )
  }
}

kw_results <- bind_rows(
  kw_results
)

dunn_results <- bind_rows(
  dunn_results
)

print(stats_marker_elevation)
print(kw_results)
print(dunn_results)

y_2d <- get_y_parameters(
  alpha_data
)

elevation_index <- setNames(
  seq_along(
    elevation_levels
  ),
  elevation_levels
)

annotation_marker_elevation <-
  stats_marker_elevation %>%
  mutate(
    x_center =
      elevation_index[
        as.character(Elevation)
      ],
    xmin =
      x_center +
      marker_offset["COI"],
    xmax =
      x_center +
      marker_offset["18S rRNA"],
    y.position =
      y_2d$y_max +
      y_2d$step
  ) %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

if (
  nrow(dunn_results) > 0
) {
  
  annotation_dunn <- dunn_results %>%
    mutate(
      xmin =
        elevation_index[
          as.character(group1)
        ] +
        marker_offset[
          as.character(Marker)
        ],
      xmax =
        elevation_index[
          as.character(group2)
        ] +
        marker_offset[
          as.character(Marker)
        ]
    ) %>%
    arrange(
      Marker,
      xmin,
      xmax
    ) %>%
    mutate(
      y.position =
        y_2d$y_max +
        y_2d$step * 2 +
        y_2d$step *
        row_number()
    )
}

p2d <- ggplot(
  alpha_data,
  aes(
    x = Elevation,
    y = Shannon,
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
    aes(color = Marker),
    position =
      position_jitterdodge(
        jitter.width = 0.2,
        dodge.width = 0.8
      ),
    alpha = 0.9,
    size = 2.25
  )

if (
  nrow(annotation_marker_elevation) > 0
) {
  p2d <- p2d +
    stat_pvalue_manual(
      annotation_marker_elevation,
      label = "significance",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

if (
  nrow(dunn_results) > 0
) {
  p2d <- p2d +
    stat_pvalue_manual(
      annotation_dunn,
      label = "significance",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

p2d <- p2d +
  labs(
    x = "Elevation (m)",
    y = "Shannon diversity"
  ) +
  scale_fill_manual(
    values = c(
      "COI" = "#F87A52",
      "18S rRNA" = "#188CBA"
    )
  ) +
  scale_color_manual(
    values = c(
      "COI" = "#F64007",
      "18S rRNA" = "#188CBA"
    )
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

print(p2d)

########## Figure 2E: soil-layer comparison ##########
marker_within_layer_wide <- alpha_data %>%
  select(
    Group,
    Sample,
    Marker,
    Shannon
  ) %>%
  pivot_wider(
    names_from = Marker,
    values_from = Shannon
  )

stats_marker_layer <-
  marker_within_layer_wide %>%
  group_by(Group) %>%
  group_modify(
    ~ paired_wilcox(
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
    significance =
      get_significance(p.adj)
  )

layer_within_marker_wide <-
  alpha_data %>%
  select(
    PlotID,
    Marker,
    Group,
    Shannon
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = Shannon
  )

stats_layer_marker <-
  layer_within_marker_wide %>%
  group_by(Marker) %>%
  group_modify(
    ~ paired_wilcox(
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
    significance =
      get_significance(p.adj)
  )

print(stats_marker_layer)
print(stats_layer_marker)

y_2e <- get_y_parameters(
  alpha_data
)

annotation_marker_layer <-
  stats_marker_layer %>%
  mutate(
    group_index = match(
      as.character(Group),
      layer_levels
    ),
    xmin =
      group_index +
      marker_offset["COI"],
    xmax =
      group_index +
      marker_offset["18S rRNA"],
    y.position =
      y_2e$y_max +
      y_2e$step
  ) %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

annotation_layer_marker <-
  stats_layer_marker %>%
  mutate(
    xmin =
      1 +
      marker_offset[
        as.character(Marker)
      ],
    xmax =
      2 +
      marker_offset[
        as.character(Marker)
      ],
    y.position =
      y_2e$y_max +
      y_2e$step * 2 +
      ifelse(
        Marker == "COI",
        0,
        y_2e$step
      )
  ) %>%
  filter(
    !is.na(p.adj),
    p.adj < 0.05
  )

p2e <- ggplot(
  alpha_data,
  aes(
    x = Group,
    y = Shannon,
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
    aes(color = Marker),
    position =
      position_jitterdodge(
        jitter.width = 0.2,
        dodge.width = 0.8
      ),
    alpha = 0.9,
    size = 2.25
  )

if (
  nrow(annotation_marker_layer) > 0
) {
  p2e <- p2e +
    stat_pvalue_manual(
      annotation_marker_layer,
      label = "significance",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

if (
  nrow(annotation_layer_marker) > 0
) {
  p2e <- p2e +
    stat_pvalue_manual(
      annotation_layer_marker,
      label = "significance",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.01,
      hide.ns = TRUE
    )
}

p2e <- p2e +
  labs(
    x = "Soil layer",
    y = "Shannon diversity"
  ) +
  scale_fill_manual(
    values = c(
      "COI" = "#F87A52",
      "18S rRNA" = "#188CBA"
    )
  ) +
  scale_color_manual(
    values = c(
      "COI" = "#F64007",
      "18S rRNA" = "#188CBA"
    )
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

print(p2e)