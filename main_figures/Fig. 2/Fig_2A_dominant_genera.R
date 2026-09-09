rm(list = ls())

########## Load required packages ##########
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggtext)
library(cowplot)

########## Read data ##########
coi_data <- read.csv("COI_data.csv", check.names = FALSE)
rrna_data <- read.csv("18S_data.csv", check.names = FALSE)
metadata <- read.csv("metadata.csv", check.names = FALSE)

########## Prepare metadata ##########
metadata <- metadata %>%
  transmute(
    SampleID = as.character(`OTU ID`),
    Layer = factor(
      Group,
      levels = c("T", "S"),
      labels = c("Topsoil", "Subsoil")
    ),
    Elevation = as.numeric(Elevation)
  )

########## Set colors and factor order ##########
marker_cols <- c(
  "COI" = "#FC7854",
  "18S rRNA" = "#5AB1CD"
)

marker_order <- c("COI", "18S rRNA")
layer_order <- c("Topsoil", "Subsoil")
group_order <- c("Bdelloidea", "Monogononta", "Other")

########## Prepare genus-level abundance data ##########
prepare_genus_data <- function(data, metadata, marker_name) {
  
  sample_cols <- intersect(
    metadata$SampleID,
    colnames(data)
  )
  
  data %>%
    mutate(
      Class_clean = sub("^[a-z]__", "", trimws(as.character(Class))),
      Order_clean = sub("^[a-z]__", "", trimws(as.character(Order))),
      Genus = sub("^[a-z]__", "", trimws(as.character(Genus))),
      RotiferGroup = case_when(
        Order_clean %in% c(
          "Adinetida",
          "Philodinida",
          "Philodinavida"
        ) ~ "Bdelloidea",
        Order_clean %in% c(
          "Ploima",
          "Flosculariida",
          "Collothecida",
          "Collothecaceae",
          "Testudinellida"
        ) ~ "Monogononta",
        Class_clean == "Bdelloidea" ~ "Bdelloidea",
        Class_clean == "Monogononta" ~ "Monogononta",
        TRUE ~ "Other"
      )
    ) %>%
    select(
      Genus,
      RotiferGroup,
      all_of(sample_cols)
    ) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "SampleID",
      values_to = "Abundance"
    ) %>%
    mutate(
      SampleID = as.character(SampleID),
      Abundance = as.numeric(Abundance)
    ) %>%
    left_join(
      metadata,
      by = "SampleID"
    ) %>%
    filter(
      !is.na(Genus),
      Genus != "",
      !is.na(Layer)
    ) %>%
    group_by(
      SampleID,
      Layer,
      Elevation,
      Genus,
      RotiferGroup
    ) %>%
    summarise(
      Abundance = sum(Abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Marker = marker_name
    )
}

coi_genus <- prepare_genus_data(
  coi_data,
  metadata,
  "COI"
)

rrna_genus <- prepare_genus_data(
  rrna_data,
  metadata,
  "18S rRNA"
)

genus_sample <- bind_rows(
  coi_genus,
  rrna_genus
) %>%
  mutate(
    Marker = as.character(Marker),
    Layer = as.character(Layer),
    RotiferGroup = ifelse(
      RotiferGroup %in% group_order,
      RotiferGroup,
      "Other"
    )
  )

########## Calculate mean relative abundance ##########
sample_marker_layer <- genus_sample %>%
  distinct(
    Marker,
    SampleID,
    Layer
  )

marker_genus_info <- genus_sample %>%
  filter(Abundance > 0) %>%
  count(
    Marker,
    Genus,
    RotiferGroup,
    sort = TRUE
  ) %>%
  group_by(
    Marker,
    Genus
  ) %>%
  slice_max(
    n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    Marker,
    Genus,
    RotiferGroup
  )

genus_rel_observed <- genus_sample %>%
  group_by(
    Marker,
    SampleID,
    Layer
  ) %>%
  mutate(
    SampleTotal = sum(Abundance, na.rm = TRUE),
    RelAbund = ifelse(
      SampleTotal > 0,
      Abundance / SampleTotal * 100,
      0
    )
  ) %>%
  ungroup() %>%
  group_by(
    Marker,
    SampleID,
    Layer,
    Genus,
    RotiferGroup
  ) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    RelAbund = sum(RelAbund, na.rm = TRUE),
    .groups = "drop"
  )

genus_rel <- sample_marker_layer %>%
  inner_join(
    marker_genus_info,
    by = "Marker"
  ) %>%
  left_join(
    genus_rel_observed %>%
      select(
        Marker,
        SampleID,
        Layer,
        Genus,
        Abundance,
        RelAbund
      ),
    by = c(
      "Marker",
      "SampleID",
      "Layer",
      "Genus"
    )
  ) %>%
  mutate(
    Abundance = replace_na(Abundance, 0),
    RelAbund = replace_na(RelAbund, 0),
    Marker = factor(
      Marker,
      levels = marker_order
    ),
    Layer = factor(
      Layer,
      levels = layer_order
    ),
    RotiferGroup = factor(
      RotiferGroup,
      levels = group_order
    )
  )

########## Select top 10 genera for each marker ##########
mean_marker_genus <- genus_rel %>%
  group_by(
    Marker,
    Genus
  ) %>%
  summarise(
    MeanRel = mean(RelAbund, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Marker = as.character(Marker)
  )

top10_by_marker <- mean_marker_genus %>%
  group_by(Marker) %>%
  slice_max(
    order_by = MeanRel,
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    Marker,
    desc(MeanRel)
  )

top10_union <- unique(
  top10_by_marker$Genus
)

########## Determine genus groups and plotting order ##########
genus_group_info <- marker_genus_info %>%
  filter(
    Genus %in% top10_union
  ) %>%
  count(
    Genus,
    RotiferGroup,
    sort = TRUE
  ) %>%
  group_by(Genus) %>%
  slice_max(
    n,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  select(
    Genus,
    RotiferGroup
  )

genus_membership <- mean_marker_genus %>%
  filter(
    Genus %in% top10_union,
    MeanRel > 0
  ) %>%
  group_by(Genus) %>%
  summarise(
    in_COI = any(Marker == "COI"),
    in_18S = any(Marker == "18S rRNA"),
    TotalRel = sum(MeanRel, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Membership = case_when(
      in_COI & in_18S ~ "Shared",
      in_COI & !in_18S ~ "COI only",
      !in_COI & in_18S ~ "18S rRNA only",
      TRUE ~ "Other"
    ),
    Membership = factor(
      Membership,
      levels = c(
        "Shared",
        "COI only",
        "18S rRNA only",
        "Other"
      )
    )
  )

genus_info <- genus_membership %>%
  left_join(
    genus_group_info,
    by = "Genus"
  ) %>%
  mutate(
    RotiferGroup = ifelse(
      is.na(RotiferGroup),
      "Other",
      as.character(RotiferGroup)
    ),
    RotiferGroup = factor(
      RotiferGroup,
      levels = group_order
    )
  ) %>%
  arrange(
    RotiferGroup,
    Membership,
    desc(TotalRel)
  ) %>%
  mutate(
    x_id = row_number()
  )

genus_order <- genus_info$Genus
genus_labels <- paste0(
  "<i>",
  genus_info$Genus,
  "</i>"
)

########## Prepare abundance-bar data ##########
bar_data <- mean_marker_genus %>%
  filter(
    Genus %in% top10_union
  ) %>%
  complete(
    Genus = genus_order,
    Marker = marker_order,
    fill = list(MeanRel = 0)
  ) %>%
  left_join(
    genus_info %>%
      select(
        Genus,
        x_id,
        RotiferGroup,
        Membership
      ),
    by = "Genus"
  ) %>%
  mutate(
    Marker = factor(
      Marker,
      levels = marker_order
    )
  )

########## Prepare detection-matrix data ##########
row_position <- tibble(
  Row = c(
    "COI | Topsoil",
    "COI | Subsoil",
    "18S rRNA | Topsoil",
    "18S rRNA | Subsoil"
  ),
  Row_id = c(4, 3, 2, 1),
  row_order_top = 1:4
)

row_bg_grey <- row_position %>%
  mutate(
    ymin = Row_id - 0.5,
    ymax = Row_id + 0.5
  ) %>%
  filter(
    row_order_top %% 2 == 0
  )

detect_data <- genus_sample %>%
  filter(
    Genus %in% top10_union
  ) %>%
  group_by(
    Marker,
    Layer,
    Genus
  ) %>%
  summarise(
    Detected = sum(Abundance, na.rm = TRUE) > 0,
    .groups = "drop"
  ) %>%
  complete(
    Marker = marker_order,
    Layer = layer_order,
    Genus = genus_order,
    fill = list(Detected = FALSE)
  ) %>%
  mutate(
    Row = paste(
      Marker,
      Layer,
      sep = " | "
    )
  ) %>%
  left_join(
    row_position,
    by = "Row"
  ) %>%
  left_join(
    genus_info %>%
      select(
        Genus,
        x_id,
        RotiferGroup,
        Membership
      ),
    by = "Genus"
  ) %>%
  mutate(
    Marker = factor(
      Marker,
      levels = marker_order
    ),
    Layer = factor(
      Layer,
      levels = layer_order
    )
  )

line_data <- detect_data %>%
  filter(Detected) %>%
  arrange(
    Genus,
    Row_id
  )

########## Calculate displayed genus numbers ##########
set_size_data <- detect_data %>%
  group_by(
    Row,
    Marker,
    Layer,
    Row_id,
    row_order_top
  ) %>%
  summarise(
    Genus_number = sum(Detected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Marker = factor(
      Marker,
      levels = marker_order
    ),
    Layer = factor(
      Layer,
      levels = layer_order
    ),
    Genus_number_left = -Genus_number
  )

max_set_size <- max(
  set_size_data$Genus_number,
  na.rm = TRUE
)

########## Prepare group separators ##########
separator_data <- genus_info %>%
  group_by(RotiferGroup) %>%
  summarise(
    max_x = max(x_id),
    .groups = "drop"
  ) %>%
  arrange(max_x)

if (nrow(separator_data) > 1) {
  separator_data <- separator_data[
    1:(nrow(separator_data) - 1),
  ] %>%
    mutate(
      xintercept = max_x + 0.5
    )
} else {
  separator_data <- tibble(
    xintercept = numeric()
  )
}

group_label_data <- genus_info %>%
  group_by(RotiferGroup) %>%
  summarise(
    xmin = min(x_id),
    xmax = max(x_id),
    xmid = mean(c(xmin, xmax)),
    .groups = "drop"
  )

########## Show results ##########
print(top10_by_marker)

print(
  genus_info %>%
    select(
      Genus,
      RotiferGroup,
      Membership,
      TotalRel
    )
)

print(
  set_size_data %>%
    select(
      Marker,
      Layer,
      Genus_number
    )
)

########## Plot upper abundance bars ##########
max_y <- max(
  bar_data$MeanRel,
  na.rm = TRUE
)

p_bar <- ggplot(
  bar_data,
  aes(
    x = x_id,
    y = MeanRel,
    fill = Marker
  )
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.68,
    color = "black",
    linewidth = 0.25
  ) +
  geom_vline(
    data = separator_data,
    aes(xintercept = xintercept),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey60"
  ) +
  geom_text(
    data = group_label_data,
    aes(
      x = xmid,
      y = max_y * 1.10,
      label = RotiferGroup
    ),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = marker_cols,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = genus_info$x_id,
    labels = genus_labels,
    expand = expansion(
      add = c(0.5, 0.5)
    )
  ) +
  scale_y_continuous(
    limits = c(0, max_y * 1.18),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    x = NULL,
    y = "Mean relative abundance (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      color = "black",
      size = 10
    ),
    axis.title.y = element_text(
      color = "black",
      size = 11
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(
      size = 10
    ),
    plot.margin = margin(
      10,
      5,
      8,
      5
    )
  )

########## Plot left genus-number bars ##########
x_breaks_pos <- pretty(
  c(0, max_set_size),
  n = 4
)

x_breaks_pos <- x_breaks_pos[
  x_breaks_pos >= 0
]

x_breaks <- -rev(
  x_breaks_pos
)

x_min <- min(x_breaks)
label_offset <- max_set_size * 0.05
bar_half_height <- 0.28

p_setsize <- ggplot(set_size_data) +
  geom_rect(
    aes(
      xmin = Genus_number_left,
      xmax = 0,
      ymin = Row_id - bar_half_height,
      ymax = Row_id + bar_half_height,
      fill = Marker
    ),
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(
      x = Genus_number_left - label_offset,
      y = Row_id,
      label = Genus_number
    ),
    hjust = 1,
    size = 3.2,
    color = "black"
  ) +
  scale_fill_manual(
    values = marker_cols,
    guide = "none"
  ) +
  scale_y_continuous(
    position = "right",
    breaks = row_position$Row_id,
    labels = NULL,
    limits = c(0.5, 4.5),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(
      x_min - label_offset * 3,
      0
    ),
    breaks = x_breaks,
    labels = abs(x_breaks),
    expand = c(0, 0)
  ) +
  labs(
    x = "Genus number",
    y = NULL
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_blank(),
    axis.line.x.bottom = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.line.y.right = element_line(
      color = "black",
      linewidth = 0.5
    ),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      color = "black",
      size = 9
    ),
    axis.ticks.x = element_line(
      color = "black",
      linewidth = 0.4
    ),
    axis.title.x = element_text(
      color = "black",
      size = 10
    ),
    panel.grid = element_blank(),
    plot.margin = margin(
      0,
      2,
      5,
      5
    )
  )

########## Plot detection matrix ##########
row_label_map <- tibble(
  Row_id = c(4, 3, 2, 1),
  label = c(
    "<span style='color:#FC7854'>COI Topsoil</span>",
    "<span style='color:#FC7854'>COI Subsoil</span>",
    "<span style='color:#5AB1CD'>18S rRNA Topsoil</span>",
    "<span style='color:#5AB1CD'>18S rRNA Subsoil</span>"
  )
)

p_matrix <- ggplot(
  detect_data,
  aes(
    x = x_id,
    y = Row_id
  )
) +
  geom_rect(
    data = row_bg_grey,
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = "grey96",
    color = NA
  ) +
  geom_vline(
    data = separator_data,
    aes(xintercept = xintercept),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey60"
  ) +
  geom_line(
    data = line_data,
    aes(group = Genus),
    color = "grey45",
    linewidth = 0.65
  ) +
  geom_point(
    data = detect_data %>%
      filter(!Detected),
    shape = 21,
    fill = "grey85",
    color = "grey85",
    size = 2.7,
    stroke = 0.2
  ) +
  geom_point(
    data = detect_data %>%
      filter(Detected),
    aes(fill = Marker),
    shape = 21,
    color = "black",
    size = 3.4,
    stroke = 0.25
  ) +
  scale_fill_manual(
    values = marker_cols,
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = genus_info$x_id,
    labels = genus_labels,
    expand = expansion(
      add = c(0.5, 0.5)
    )
  ) +
  scale_y_continuous(
    breaks = row_label_map$Row_id,
    labels = row_label_map$label,
    limits = c(0.5, 4.5),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_blank(),
    axis.text.x = element_markdown(
      size = 9,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_markdown(
      size = 10
    ),
    axis.ticks = element_blank(),
    plot.margin = margin(
      0,
      5,
      5,
      2
    )
  )

########## Combine Figure 2A ##########
aligned_right <- align_plots(
  p_bar,
  p_matrix,
  align = "v",
  axis = "lr"
)

top_row <- plot_grid(
  NULL,
  aligned_right[[1]],
  ncol = 2,
  rel_widths = c(1, 6),
  align = "h"
)

bottom_row <- plot_grid(
  p_setsize,
  aligned_right[[2]],
  ncol = 2,
  rel_widths = c(1, 6),
  align = "h"
)

p2a <- plot_grid(
  top_row,
  bottom_row,
  ncol = 1,
  rel_heights = c(2.25, 1.35),
  align = "v"
)

print(p2a)