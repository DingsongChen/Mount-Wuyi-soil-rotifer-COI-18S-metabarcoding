rm(list = ls())

########## Load required packages ##########
library(vegan)
library(permute)
library(ggplot2)
library(ggrepel)
library(rdacca.hp)
library(dplyr)

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

metadata$Sample <- as.character(
  metadata[["OTU ID"]]
)

########## Analysis settings ##########
nperm <- 9999
analysis_seed <- 20260720

marker_colours <- c(
  "COI" = "#F87A52",
  "18S rRNA" = "#5EAFCF"
)

panel_aspect_ratio <- 309 / 331
bar_outline_1pt <- 25.4 / 72

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

display_names <- c(
  "pH" = "pH",
  "SWC" = "SWC",
  "NO3_N" = "NO3-N",
  "NH4_N" = "NH4-N",
  "TC" = "TC",
  "TP" = "TP",
  "AP" = "AP",
  "MBN" = "MBN",
  "AcP" = "AcP",
  "LAP" = "LAP"
)

########## Helper functions ##########
significance_symbol <- function(p) {
  case_when(
    is.na(p) ~ "",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ ""
  )
}

extract_hp_contribution <- function(hp_table) {

  hp_df <- as.data.frame(
    hp_table,
    check.names = FALSE
  )

  hp_df$Variable <- rownames(hp_df)
  rownames(hp_df) <- NULL

  if (!"Individual" %in% colnames(hp_df)) {
    stop(
      "The 'Individual' column was not found in rdacca.hp output."
    )
  }

  tibble(
    Variable = hp_df$Variable,
    Individual_R2 = as.numeric(
      hp_df[["Individual"]]
    ),
    Contribution_percent =
      as.numeric(
        hp_df[["Individual"]]
      ) * 100
  )
}

prepare_marker_data <- function(
    data,
    metadata
) {

  sample_names <- intersect(
    metadata$Sample,
    colnames(data)
  )

  metadata_marker <- metadata[
    match(
      sample_names,
      metadata$Sample
    ),
    ,
    drop = FALSE
  ]

  env_raw <- data.frame(
    row.names = sample_names,
    check.names = FALSE
  )

  for (internal_name in names(environment_columns)) {

    source_name <- environment_columns[[internal_name]]

    env_raw[[internal_name]] <- as.numeric(
      metadata_marker[[source_name]]
    )
  }

  community <- t(
    as.matrix(
      data[
        ,
        sample_names,
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

  keep <- complete.cases(env_raw) &
    rowSums(
      community,
      na.rm = TRUE
    ) > 0

  env_raw <- env_raw[
    keep,
    ,
    drop = FALSE
  ]

  community <- community[
    keep,
    ,
    drop = FALSE
  ]

  community <- community[
    ,
    colSums(
      community,
      na.rm = TRUE
    ) > 0,
    drop = FALSE
  ]

  if (
    !identical(
      rownames(env_raw),
      rownames(community)
    )
  ) {
    stop(
      "Sample order does not match between environment and community matrices."
    )
  }

  zero_variance <- names(
    which(
      vapply(
        env_raw,
        function(x) sd(
          x,
          na.rm = TRUE
        ) == 0,
        logical(1)
      )
    )
  )

  if (length(zero_variance) > 0) {
    stop(
      paste(
        "Zero-variance environmental variables:",
        paste(
          zero_variance,
          collapse = ", "
        )
      )
    )
  }

  env_scaled <- as.data.frame(
    scale(env_raw),
    check.names = FALSE
  )

  rownames(env_scaled) <- rownames(env_raw)

  list(
    community = community,
    environment = env_scaled
  )
}

########## Run dbRDA, envfit and hierarchical partitioning ##########
run_marker_analysis <- function(
    data,
    metadata,
    marker_name
) {

  prepared <- prepare_marker_data(
    data,
    metadata
  )

  community <- prepared$community
  env_scaled <- prepared$environment

  bray_distance <- vegdist(
    community,
    method = "bray"
  )

  db_rda_model <- capscale(
    bray_distance ~ .,
    data = env_scaled
  )

  set.seed(analysis_seed)

  permutation_matrix <- shuffleSet(
    n = nrow(env_scaled),
    nset = nperm,
    control = how(
      nperm = nperm
    )
  )

  overall_test <- anova(
    db_rda_model,
    permutations = permutation_matrix
  )

  environmental_fit <- envfit(
    db_rda_model,
    env_scaled,
    permutations = permutation_matrix
  )

  environmental_vectors <- as.data.frame(
    scores(
      environmental_fit,
      display = "vectors",
      choices = 1:2
    )
  )

  colnames(
    environmental_vectors
  )[1:2] <- c(
    "dbRDA1",
    "dbRDA2"
  )

  environmental_vectors$Variable <- rownames(
    environmental_vectors
  )

  environmental_vectors$r2 <- as.numeric(
    environmental_fit$vectors$r[
      environmental_vectors$Variable
    ]
  )

  environmental_vectors$P_value <- as.numeric(
    environmental_fit$vectors$pvals[
      environmental_vectors$Variable
    ]
  )

  environmental_vectors$Significance <- vapply(
    environmental_vectors$P_value,
    significance_symbol,
    character(1)
  )

  environmental_vectors$Display_name <- unname(
    display_names[
      environmental_vectors$Variable
    ]
  )

  environmental_vectors$Label <- paste0(
    environmental_vectors$Display_name,
    environmental_vectors$Significance
  )

  significant_vectors <- environmental_vectors %>%
    filter(
      !is.na(P_value),
      P_value < 0.05
    )

  site_scores <- as.data.frame(
    scores(
      db_rda_model,
      display = "sites",
      choices = 1:2
    )
  )

  colnames(
    site_scores
  )[1:2] <- c(
    "dbRDA1",
    "dbRDA2"
  )

  constrained_eigenvalues <- db_rda_model$CCA$eig

  axis1_percent <- round(
    constrained_eigenvalues[1] /
      sum(constrained_eigenvalues) *
      100,
    2
  )

  axis2_percent <- round(
    constrained_eigenvalues[2] /
      sum(constrained_eigenvalues) *
      100,
    2
  )

  total_constrained_percent <- round(
    sum(constrained_eigenvalues) /
      db_rda_model$tot.chi *
      100,
    2
  )

  if (nrow(significant_vectors) > 0) {

    site_x_limit <- max(
      abs(site_scores$dbRDA1),
      na.rm = TRUE
    )

    site_y_limit <- max(
      abs(site_scores$dbRDA2),
      na.rm = TRUE
    )

    vector_x_limit <- max(
      abs(significant_vectors$dbRDA1),
      na.rm = TRUE
    )

    vector_y_limit <- max(
      abs(significant_vectors$dbRDA2),
      na.rm = TRUE
    )

    possible_scales <- c(
      if (vector_x_limit > 0) {
        site_x_limit / vector_x_limit
      } else {
        NA_real_
      },
      if (vector_y_limit > 0) {
        site_y_limit / vector_y_limit
      } else {
        NA_real_
      }
    )

    possible_scales <- possible_scales[
      is.finite(possible_scales) &
        possible_scales > 0
    ]

    arrow_multiplier <- if (
      length(possible_scales) == 0
    ) {
      1
    } else {
      0.70 * min(possible_scales)
    }

    significant_vectors <- significant_vectors %>%
      mutate(
        Plot_x = dbRDA1 * arrow_multiplier,
        Plot_y = dbRDA2 * arrow_multiplier,
        Label_x = Plot_x * 1.15,
        Label_y = Plot_y * 1.15
      )
  }

  db_rda_plot <- ggplot(
    site_scores,
    aes(
      x = dbRDA1,
      y = dbRDA2
    )
  ) +
    geom_point(
      size = 5,
      alpha = 0.8,
      shape = 16,
      colour = marker_colours[[marker_name]]
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      colour = "black"
    ) +
    labs(
      x = paste0(
        "dbRDA1 (",
        axis1_percent,
        "%)"
      ),
      y = paste0(
        "dbRDA2 (",
        axis2_percent,
        "%)"
      ),
      title = paste0(
        marker_name,
        " rotifer communities"
      ),
      subtitle = paste0(
        "Total constrained variance: ",
        total_constrained_percent,
        "%"
      )
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
      aspect.ratio = panel_aspect_ratio,
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      )
    ) +
    coord_cartesian(
      clip = "off"
    )

  if (nrow(significant_vectors) > 0) {

    db_rda_plot <- db_rda_plot +
      geom_segment(
        data = significant_vectors,
        aes(
          x = 0,
          y = 0,
          xend = Plot_x,
          yend = Plot_y
        ),
        inherit.aes = FALSE,
        arrow = arrow(
          length = grid::unit(
            0.02,
            "npc"
          )
        ),
        colour = "black",
        linewidth = 0.8
      ) +
      geom_text_repel(
        data = significant_vectors,
        aes(
          x = Label_x,
          y = Label_y,
          label = Label
        ),
        inherit.aes = FALSE,
        size = 3.5,
        colour = "black",
        bg.color = "white",
        bg.r = 0.15,
        max.overlaps = 50,
        box.padding = 0.5,
        point.padding = 0.5,
        min.segment.length = 0.2
      )
  }

  hierarchical_result <- rdacca.hp(
    bray_distance,
    env_scaled,
    method = "dbRDA",
    type = "R2",
    scale = FALSE
  )

  contribution_data <- extract_hp_contribution(
    hierarchical_result$Hier.part
  ) %>%
    mutate(
      Marker = marker_name,
      Factor = unname(
        display_names[Variable]
      )
    )

  adjusted_r2 <- RsquareAdj(
    db_rda_model
  )

  overall_p <- as.data.frame(
    overall_test,
    check.names = FALSE
  )[1, "Pr(>F)"]

  model_summary <- tibble(
    Marker = marker_name,
    Samples = nrow(community),
    ASVs = ncol(community),
    R2 = unname(
      adjusted_r2$r.squared
    ),
    Adjusted_R2 = unname(
      adjusted_r2$adj.r.squared
    ),
    Total_constrained_percent =
      total_constrained_percent,
    Axis1_percent_of_constrained =
      axis1_percent,
    Axis2_percent_of_constrained =
      axis2_percent,
    Overall_P = overall_p
  )

  list(
    plot = db_rda_plot,
    envfit = environmental_vectors,
    significant_vectors = significant_vectors,
    contribution = contribution_data,
    model_summary = model_summary
  )
}

########## Figure 4A: COI dbRDA ##########
coi_result <- run_marker_analysis(
  coi_data,
  metadata,
  "COI"
)

########## Figure 4B: 18S rRNA dbRDA ##########
rrna_result <- run_marker_analysis(
  rrna_data,
  metadata,
  "18S rRNA"
)

########## Show dbRDA results ##########
cat("\nCOI dbRDA model summary:\n")
print(coi_result$model_summary)

cat("\n18S rRNA dbRDA model summary:\n")
print(rrna_result$model_summary)

cat("\nCOI envfit results:\n")
print(coi_result$envfit)

cat("\n18S rRNA envfit results:\n")
print(rrna_result$envfit)

########## Figure 4C: hierarchical partitioning ##########
contribution_plot_data <- bind_rows(
  coi_result$contribution,
  rrna_result$contribution
)

contribution_plot_data$Marker <- factor(
  contribution_plot_data$Marker,
  levels = c(
    "COI",
    "18S rRNA"
  )
)

factor_order <- contribution_plot_data %>%
  group_by(Factor) %>%
  summarise(
    Mean_contribution = mean(
      Contribution_percent,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(Mean_contribution)
  ) %>%
  pull(Factor)

contribution_plot_data$Factor <- factor(
  contribution_plot_data$Factor,
  levels = factor_order
)

fig4c_y_upper <- max(
  contribution_plot_data$Contribution_percent,
  na.rm = TRUE
) / 0.8

fig4c_plot <- ggplot(
  contribution_plot_data,
  aes(
    x = Factor,
    y = Contribution_percent,
    fill = Marker
  )
) +
  geom_col(
    position = position_dodge(
      width = 0.78
    ),
    width = 0.70,
    colour = "black",
    linewidth = bar_outline_1pt
  ) +
  scale_fill_manual(
    values = marker_colours,
    breaks = c(
      "COI",
      "18S rRNA"
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      fig4c_y_upper
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  labs(
    x = NULL,
    y = "Independent contribution to R² (%)",
    fill = NULL
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    ),
    aspect.ratio = panel_aspect_ratio,
    legend.position = "none",
    legend.title = element_blank(),
    axis.title = element_text(
      colour = "black"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    )
  )

########## Display Figure 4A-C ##########
print(coi_result$plot)
print(rrna_result$plot)
print(fig4c_plot)
