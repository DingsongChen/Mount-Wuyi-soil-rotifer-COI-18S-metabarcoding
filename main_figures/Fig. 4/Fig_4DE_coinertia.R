rm(list = ls())

########## Load required packages ##########
library(ade4)
library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)

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

########## Analysis settings ##########
nperm <- 9999
analysis_seed <- 20260804
number_of_axes <- 2
panel_aspect_ratio <- 309 / 331

elevation_levels <- c(
  300, 600, 900, 1200,
  1600, 1800, 2100
)

elevation_colours <- c(
  "300" = "#2C6AA0",
  "600" = "#4B91C3",
  "900" = "#2A9D8F",
  "1200" = "#6D9F38",
  "1600" = "#C9821C",
  "1800" = "#D95F2D",
  "2100" = "#A6294B"
)

environment_columns <- c(
  "pH" = "pH",
  "SWC" = "SWC",
  "NO3_N" = "NO3--N",
  "NH4_N" = "NH4+-N",
  "TC" = "TC",
  "TP" = "TP",
  "AP" = "AP",
  "MBN" = "MBN",
  "AcP" = "AcP",
  "LAP" = "LAP"
)

########## Prepare metadata ##########
metadata <- metadata %>%
  transmute(
    Sample = as.character(`OTU ID`),
    Group = factor(
      Group,
      levels = c("T", "S"),
      labels = c("Topsoil", "Subsoil")
    ),
    Elevation = as.numeric(Elevation),
    PlotID = sub(
      "^[TtSs][_-]?",
      "",
      as.character(`OTU ID`)
    ),
    pH = as.numeric(pH),
    SWC = as.numeric(SWC),
    NO3_N = as.numeric(`NO3--N`),
    NH4_N = as.numeric(`NH4+-N`),
    TC = as.numeric(TC),
    TP = as.numeric(TP),
    AP = as.numeric(AP),
    MBN = as.numeric(MBN),
    AcP = as.numeric(AcP),
    LAP = as.numeric(LAP)
  )

########## Retain common complete plot pairs ##########
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

coi_samples <- setdiff(
  colnames(coi_data),
  tax_cols
)

rrna_samples <- setdiff(
  colnames(rrna_data),
  tax_cols
)

common_samples <- Reduce(
  intersect,
  list(
    coi_samples,
    rrna_samples,
    metadata$Sample
  )
)

metadata_final <- metadata %>%
  filter(
    Sample %in% common_samples
  ) %>%
  filter(
    complete.cases(
      across(
        c(
          Elevation,
          all_of(
            names(environment_columns)
          )
        )
      )
    )
  )

complete_plot_ids <- metadata_final %>%
  group_by(PlotID) %>%
  summarise(
    n_samples = n(),
    n_layers = n_distinct(Group),
    n_elevations = n_distinct(Elevation),
    .groups = "drop"
  ) %>%
  filter(
    n_samples == 2,
    n_layers == 2,
    n_elevations == 1
  ) %>%
  pull(PlotID)

metadata_final <- metadata_final %>%
  filter(
    PlotID %in% complete_plot_ids
  ) %>%
  arrange(
    PlotID,
    Group
  ) %>%
  droplevels()

final_samples <- metadata_final$Sample

cat(
  "Common complete plots:",
  n_distinct(metadata_final$PlotID),
  "\nSamples per marker:",
  nrow(metadata_final),
  "\n"
)

########## Prepare environmental matrix ##########
environment_raw <- metadata_final %>%
  select(
    all_of(
      names(environment_columns)
    )
  ) %>%
  as.data.frame()

rownames(environment_raw) <- final_samples

########## Prepare community matrices ##########
prepare_community <- function(
    data,
    samples
) {

  community <- t(
    as.matrix(
      data[
        ,
        samples,
        drop = FALSE
      ]
    )
  )

  storage.mode(community) <- "numeric"

  community <- community[
    ,
    colSums(
      community,
      na.rm = TRUE
    ) > 0,
    drop = FALSE
  ]

  community
}

coi_community <- prepare_community(
  coi_data,
  final_samples
)

rrna_community <- prepare_community(
  rrna_data,
  final_samples
)

coi_hellinger <- decostand(
  coi_community,
  method = "hellinger"
)

rrna_hellinger <- decostand(
  rrna_community,
  method = "hellinger"
)

########## Co-inertia helper functions ##########
rv_from_gram <- function(
    gram_x,
    gram_y
) {

  numerator <- sum(
    gram_x * gram_y
  )

  denominator <- sqrt(
    sum(gram_x^2) *
      sum(gram_y^2)
  )

  numerator / denominator
}

community_row_space_scores <- function(
    hellinger_matrix,
    tolerance = 1e-10
) {

  centered <- scale(
    hellinger_matrix,
    center = TRUE,
    scale = FALSE
  )

  gram <- tcrossprod(
    centered
  )

  gram <- (
    gram +
      t(gram)
  ) / 2

  eigen_result <- eigen(
    gram,
    symmetric = TRUE
  )

  keep_axes <- which(
    eigen_result$values >
      max(eigen_result$values) *
      tolerance
  )

  scores <- sweep(
    eigen_result$vectors[
      ,
      keep_axes,
      drop = FALSE
    ],
    2,
    sqrt(
      eigen_result$values[
        keep_axes
      ]
    ),
    FUN = "*"
  )

  scores <- scale(
    scores,
    center = TRUE,
    scale = FALSE
  )

  rownames(scores) <- rownames(
    hellinger_matrix
  )

  scores
}

fit_coinertia <- function(
    environment_raw,
    community_hellinger,
    marker_name
) {

  community_scores <- community_row_space_scores(
    community_hellinger
  )

  environment_pca <- dudi.pca(
    environment_raw,
    center = TRUE,
    scale = TRUE,
    scannf = FALSE,
    nf = min(
      ncol(environment_raw),
      nrow(environment_raw) - 1
    )
  )

  community_pca <- dudi.pca(
    as.data.frame(
      community_scores
    ),
    center = FALSE,
    scale = FALSE,
    scannf = FALSE,
    nf = min(
      ncol(community_scores),
      nrow(community_scores) - 1
    )
  )

  coinertia_model <- coinertia(
    environment_pca,
    community_pca,
    scannf = FALSE,
    nf = number_of_axes
  )

  environment_gram <- tcrossprod(
    as.matrix(
      environment_pca$tab
    )
  )

  community_gram <- tcrossprod(
    as.matrix(
      community_pca$tab
    )
  )

  direct_rv <- rv_from_gram(
    environment_gram,
    community_gram
  )

  if (
    abs(
      coinertia_model$RV -
        direct_rv
    ) > 1e-7
  ) {
    stop(
      paste0(
        marker_name,
        ": direct and ade4 RV coefficients do not match."
      )
    )
  }

  list(
    marker = marker_name,
    model = coinertia_model,
    environment_gram = environment_gram,
    community_gram = community_gram
  )
}

########## Plot-block permutations ##########
make_plot_block_permutations <- function(
    metadata,
    nperm,
    seed
) {

  plot_blocks <- split(
    seq_len(
      nrow(metadata)
    ),
    as.character(
      metadata$PlotID
    )
  )

  set.seed(seed)

  permutation_matrix <- matrix(
    NA_integer_,
    nrow = nperm,
    ncol = nrow(metadata)
  )

  for (i in seq_len(nperm)) {

    permuted_plots <- sample(
      names(plot_blocks),
      length(plot_blocks),
      replace = FALSE
    )

    permutation_matrix[i, ] <- unlist(
      plot_blocks[
        permuted_plots
      ],
      use.names = FALSE
    )
  }

  permutation_matrix
}

run_rv_permutation <- function(
    environment_gram,
    community_gram,
    permutation_matrix
) {

  observed_rv <- rv_from_gram(
    environment_gram,
    community_gram
  )

  null_rv <- apply(
    permutation_matrix,
    1,
    function(index) {

      permuted_community <- community_gram[
        index,
        index,
        drop = FALSE
      ]

      rv_from_gram(
        environment_gram,
        permuted_community
      )
    }
  )

  p_value <- (
    1 +
      sum(
        null_rv >= observed_rv
      )
  ) /
    (
      length(null_rv) +
        1
    )

  list(
    RV = observed_rv,
    P = p_value
  )
}

########## Fit co-inertia models ##########
coi_model <- fit_coinertia(
  environment_raw,
  coi_hellinger,
  "COI"
)

rrna_model <- fit_coinertia(
  environment_raw,
  rrna_hellinger,
  "18S rRNA"
)

permutation_matrix <- make_plot_block_permutations(
  metadata_final,
  nperm,
  analysis_seed
)

coi_test <- run_rv_permutation(
  coi_model$environment_gram,
  coi_model$community_gram,
  permutation_matrix
)

rrna_test <- run_rv_permutation(
  rrna_model$environment_gram,
  rrna_model$community_gram,
  permutation_matrix
)

########## Extract plotting data ##########
extract_coinertia_scores <- function(
    model_object,
    metadata,
    permutation_test
) {

  model <- model_object$model

  eigenvalues <- as.numeric(
    model$eig
  )

  axis_percent <- 100 *
    eigenvalues /
    sum(eigenvalues)

  environment_scores <- as.data.frame(
    model$mX[
      ,
      1:2,
      drop = FALSE
    ]
  )

  community_scores <- as.data.frame(
    model$mY[
      ,
      1:2,
      drop = FALSE
    ]
  )

  sample_scores <- tibble(
    Sample = rownames(
      environment_scores
    ),
    Environment_axis1 =
      environment_scores[, 1],
    Environment_axis2 =
      environment_scores[, 2],
    Community_axis1 =
      community_scores[, 1],
    Community_axis2 =
      community_scores[, 2]
  ) %>%
    left_join(
      metadata %>%
        select(
          Sample,
          Elevation
        ),
      by = "Sample"
    )

  list(
    scores = sample_scores,
    summary = tibble(
      Marker = model_object$marker,
      RV = permutation_test$RV,
      P = permutation_test$P,
      Axis1_percent = axis_percent[1],
      Axis2_percent = axis_percent[2]
    )
  )
}

coi_result <- extract_coinertia_scores(
  coi_model,
  metadata_final,
  coi_test
)

rrna_result <- extract_coinertia_scores(
  rrna_model,
  metadata_final,
  rrna_test
)

model_summary <- bind_rows(
  coi_result$summary,
  rrna_result$summary
)

cat("\nCo-inertia model summary:\n")
print(model_summary)

########## Shared axis limits ##########
all_axis1 <- c(
  coi_result$scores$Environment_axis1,
  coi_result$scores$Community_axis1,
  rrna_result$scores$Environment_axis1,
  rrna_result$scores$Community_axis1
)

all_axis2 <- c(
  coi_result$scores$Environment_axis2,
  coi_result$scores$Community_axis2,
  rrna_result$scores$Environment_axis2,
  rrna_result$scores$Community_axis2
)

axis1_range <- range(
  all_axis1,
  finite = TRUE
)

axis2_range <- range(
  all_axis2,
  finite = TRUE
)

axis1_padding <- max(
  diff(axis1_range) * 0.08,
  0.10
)

axis2_padding <- max(
  diff(axis2_range) * 0.08,
  0.10
)

shared_x_limits <- c(
  axis1_range[1] - axis1_padding,
  axis1_range[2] + axis1_padding
)

shared_y_limits <- c(
  axis2_range[1] - axis2_padding,
  axis2_range[2] + axis2_padding
)

########## Plotting function ##########
format_p_value <- function(p) {
  if (p < 0.001) {
    "P < 0.001"
  } else {
    paste0(
      "P = ",
      format(
        round(
          p,
          3
        ),
        nsmall = 3
      )
    )
  }
}

plot_coinertia <- function(
    result,
    marker_name
) {

  scores <- result$scores
  summary <- result$summary

  scores$Elevation_factor <- factor(
    as.character(
      scores$Elevation
    ),
    levels = as.character(
      elevation_levels
    ),
    ordered = TRUE
  )

  point_data <- bind_rows(
    tibble(
      Axis1 = scores$Environment_axis1,
      Axis2 = scores$Environment_axis2,
      Position = "Environment",
      Elevation_factor =
        scores$Elevation_factor
    ),
    tibble(
      Axis1 = scores$Community_axis1,
      Axis2 = scores$Community_axis2,
      Position = "Community",
      Elevation_factor =
        scores$Elevation_factor
    )
  )

  point_data$Position <- factor(
    point_data$Position,
    levels = c(
      "Environment",
      "Community"
    )
  )

  annotation_text <- paste0(
    "RV = ",
    format(
      round(
        summary$RV,
        3
      ),
      nsmall = 3
    ),
    "\n",
    format_p_value(
      summary$P
    )
  )

  ggplot() +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.40,
      colour = "grey60"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.40,
      colour = "grey60"
    ) +
    geom_segment(
      data = scores,
      aes(
        x = Environment_axis1,
        y = Environment_axis2,
        xend = Community_axis1,
        yend = Community_axis2,
        colour = Elevation_factor
      ),
      linewidth = 0.60,
      alpha = 0.62,
      arrow = grid::arrow(
        length = grid::unit(
          0.075,
          "inches"
        ),
        type = "closed"
      )
    ) +
    geom_point(
      data = point_data,
      aes(
        x = Axis1,
        y = Axis2,
        colour = Elevation_factor,
        shape = Position
      ),
      size = 2.8,
      stroke = 0.80,
      alpha = 0.95
    ) +
    scale_shape_manual(
      values = c(
        "Environment" = 1,
        "Community" = 16
      )
    ) +
    scale_colour_manual(
      values = elevation_colours,
      breaks = names(
        elevation_colours
      ),
      name = "Elevation (m)",
      drop = FALSE
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = annotation_text,
      hjust = 1.08,
      vjust = 1.12,
      size = 4.0,
      family = "serif",
      colour = "black"
    ) +
    labs(
      x = paste0(
        "Co-inertia axis 1 (",
        format(
          round(
            summary$Axis1_percent,
            1
          ),
          nsmall = 1
        ),
        "%)"
      ),
      y = paste0(
        "Co-inertia axis 2 (",
        format(
          round(
            summary$Axis2_percent,
            1
          ),
          nsmall = 1
        ),
        "%)"
      ),
      title = marker_name,
      shape = NULL
    ) +
    coord_cartesian(
      xlim = shared_x_limits,
      ylim = shared_y_limits,
      expand = FALSE,
      clip = "off"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      ),
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text = element_text(
        colour = "black"
      ),
      axis.title = element_text(
        colour = "black"
      ),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key.width = grid::unit(
        0.75,
        "lines"
      ),
      aspect.ratio = panel_aspect_ratio,
      plot.margin = margin(
        8,
        10,
        8,
        8
      )
    )
}

########## Figure 4D: COI co-inertia ##########
p_coi <- plot_coinertia(
  coi_result,
  "COI"
)

########## Figure 4E: 18S rRNA co-inertia ##########
p_rrna <- plot_coinertia(
  rrna_result,
  "18S rRNA"
)

########## Display plots ##########
fig4de <- wrap_plots(
  p_coi,
  p_rrna,
  ncol = 2,
  guides = "collect"
) &
  theme(
    legend.position = "bottom"
  )

print(p_coi)
print(p_rrna)
print(fig4de)
