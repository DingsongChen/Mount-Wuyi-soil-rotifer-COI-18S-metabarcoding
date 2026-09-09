rm(list = ls())

########## Load required packages ##########
library(vegan)
library(ggplot2)
library(cowplot)

########## Global settings ##########
set.seed(123)

mantel_permutations <- 9999
main_to_density_ratio <- 6

########## Helper functions ##########
format_p <- function(p) {
  if (is.na(p)) {
    return("P = NA")
  }

  if (p < 0.001) {
    return("P < 0.001")
  }

  paste0(
    "P = ",
    format(
      round(p, 3),
      nsmall = 3
    )
  )
}

safe_numeric_df <- function(dat) {
  dat <- as.data.frame(dat)
  dat[] <- lapply(
    dat,
    function(x) as.numeric(as.character(x))
  )
  dat[is.na(dat)] <- 0
  dat
}

clamp_value <- function(x, lower, upper) {
  pmin(
    pmax(x, lower),
    upper
  )
}

fill_missing_cols <- function(mat, all_cols) {
  missing_cols <- setdiff(
    all_cols,
    colnames(mat)
  )

  if (length(missing_cols) > 0) {
    add_mat <- matrix(
      0,
      nrow = nrow(mat),
      ncol = length(missing_cols),
      dimnames = list(
        rownames(mat),
        missing_cols
      )
    )

    mat <- cbind(
      mat,
      add_mat
    )
  }

  mat[
    ,
    all_cols,
    drop = FALSE
  ]
}

########## Original plotting function ##########
# Plotting parameters and patchwork layout are retained from
# the original Fig. S3 script without modification.
plot_concordance_density <- function(
    plot_df,
    x_lab,
    y_lab,
    label_text,
    point_color = "#5EAFCF",
    fixed_range = NULL,
    breaks = NULL,
    clamp_prediction = FALSE,
    main_to_density_ratio = 6
) {

  plot_df <- plot_df[
    complete.cases(plot_df),
    ,
    drop = FALSE
  ]

  if (is.null(fixed_range)) {

    range_all <- range(
      c(
        plot_df$x,
        plot_df$y
      ),
      na.rm = TRUE
    )

    range_span <- diff(range_all)

    if (range_span == 0) {
      range_span <- 1
    }

    plot_min <- range_all[1] - 0.05 * range_span
    plot_max <- range_all[2] + 0.05 * range_span

  } else {

    plot_min <- fixed_range[1]
    plot_max <- fixed_range[2]
  }

  if (is.null(breaks)) {
    breaks <- pretty(
      c(plot_min, plot_max),
      n = 5
    )
  }

  lm_fit <- lm(
    y ~ x,
    data = plot_df
  )

  fit_df <- data.frame(
    x = seq(
      plot_min,
      plot_max,
      length.out = 200
    )
  )

  fit_pred <- predict(
    lm_fit,
    newdata = fit_df,
    interval = "confidence",
    level = 0.95
  )

  fit_df$fit <- fit_pred[, "fit"]
  fit_df$lwr <- fit_pred[, "lwr"]
  fit_df$upr <- fit_pred[, "upr"]

  if (clamp_prediction) {
    fit_df$fit <- clamp_value(
      fit_df$fit,
      plot_min,
      plot_max
    )
    fit_df$lwr <- clamp_value(
      fit_df$lwr,
      plot_min,
      plot_max
    )
    fit_df$upr <- clamp_value(
      fit_df$upr,
      plot_min,
      plot_max
    )
  }

  p_main <- ggplot(
    plot_df,
    aes(
      x = x,
      y = y
    )
  ) +
    geom_point(
      size = 2.6,
      alpha = 0.55,
      color = point_color,
      fill = point_color,
      shape = 21,
      stroke = 0.4
    ) +
    geom_ribbon(
      data = fit_df,
      aes(
        x = x,
        ymin = lwr,
        ymax = upr
      ),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.6
    ) +
    geom_line(
      data = fit_df,
      aes(
        x = x,
        y = fit
      ),
      inherit.aes = FALSE,
      color = point_color,
      linewidth = 1.2
    ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50",
      linewidth = 0.8
    ) +
    annotate(
      "text",
      x = plot_min + 0.04 * (plot_max - plot_min),
      y = plot_max - 0.04 * (plot_max - plot_min),
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 5.5,
      family = "serif"
    ) +
    scale_x_continuous(
      limits = c(plot_min, plot_max),
      breaks = breaks,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(plot_min, plot_max),
      breaks = breaks,
      expand = c(0, 0)
    ) +
    coord_fixed(
      ratio = 1,
      xlim = c(plot_min, plot_max),
      ylim = c(plot_min, plot_max),
      expand = FALSE
    ) +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    theme_bw(
      base_size = 16
    ) +
    theme(
      text = element_text(
        family = "serif",
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold"
      ),
      axis.text = element_text(
        size = 13,
        color = "black"
      ),
      panel.grid = element_line(
        color = "grey90",
        linewidth = 0.4
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 1
      ),
      plot.margin = margin(
        0,
        0,
        0,
        0
      )
    )

  # Top marginal density: identical to the original template.
  p_top <- ggplot(
    plot_df,
    aes(x = x)
  ) +
    geom_density(
      color = point_color,
      fill = point_color,
      alpha = 0.35,
      linewidth = 1,
      adjust = 1
    ) +
    scale_x_continuous(
      limits = c(plot_min, plot_max),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    coord_cartesian(
      xlim = c(plot_min, plot_max)
    ) +
    theme_void() +
    theme(
      plot.margin = margin(
        0,
        0,
        0,
        0
      )
    )

  # Right marginal density: identical to the original template.
  # coord_flip() is retained so the density extends horizontally
  # from the y-axis range while sample values remain vertical.
  p_right <- ggplot(
    plot_df,
    aes(x = y)
  ) +
    geom_density(
      color = point_color,
      fill = point_color,
      alpha = 0.35,
      linewidth = 1,
      adjust = 1
    ) +
    scale_x_continuous(
      limits = c(plot_min, plot_max),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      expand = c(0, 0)
    ) +
    coord_flip(
      xlim = c(plot_min, plot_max)
    ) +
    theme_void() +
    theme(
      plot.margin = margin(
        0,
        0,
        0,
        0
      )
    )

  # Attach the original marginal plots to the actual plotting panel.
  # The marginal ranges therefore follow x = min:max and y = min:max
  # of the black plotting frame, rather than the outer axis text area.
  p_final <- cowplot::insert_xaxis_grob(
    p_main,
    ggplotGrob(p_top),
    height = grid::unit(
      1 / main_to_density_ratio,
      "null"
    ),
    position = "top"
  )

  p_final <- cowplot::insert_yaxis_grob(
    p_final,
    ggplotGrob(p_right),
    width = grid::unit(
      1 / main_to_density_ratio,
      "null"
    ),
    position = "right"
  )

  cowplot::ggdraw(p_final)
}

########## One-marker analysis ##########
run_rarefaction_robustness <- function(
    input_file,
    dataset_name,
    rarefaction_depth,
    point_color
) {

  # Reset the seed so each marker reproduces the original
  # single-marker workflow.
  set.seed(123)

  otu_raw <- read.csv(
    input_file,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )

  otu_raw <- safe_numeric_df(
    otu_raw
  )

  otu_raw <- otu_raw[
    rowSums(
      otu_raw,
      na.rm = TRUE
    ) > 0,
    ,
    drop = FALSE
  ]

  if (any(otu_raw < 0, na.rm = TRUE)) {
    stop(
      paste0(
        dataset_name,
        ": negative values were found in the input table."
      )
    )
  }

  otu_sample <- t(
    otu_raw
  )

  otu_sample <- otu_sample[
    rowSums(
      otu_sample,
      na.rm = TRUE
    ) > 0,
    ,
    drop = FALSE
  ]

  otu_sample <- otu_sample[
    ,
    colSums(
      otu_sample,
      na.rm = TRUE
    ) > 0,
    drop = FALSE
  ]

  sample_depths <- rowSums(
    otu_sample
  )

  if (any(sample_depths < rarefaction_depth)) {
    stop(
      paste0(
        dataset_name,
        ": at least one sample has fewer reads than the selected rarefaction depth."
      )
    )
  }

  cat(
    "\nDataset: ",
    dataset_name,
    "\nSamples: ",
    nrow(otu_sample),
    "\nASVs: ",
    ncol(otu_sample),
    "\nMinimum reads: ",
    min(sample_depths),
    "\nMaximum reads: ",
    max(sample_depths),
    "\nRarefaction depth: ",
    rarefaction_depth,
    "\n",
    sep = ""
  )

  ########## Rarefaction ##########
  otu_rarefied <- vegan::rrarefy(
    otu_sample,
    sample = rarefaction_depth
  )

  otu_rarefied <- otu_rarefied[
    rowSums(
      otu_rarefied,
      na.rm = TRUE
    ) > 0,
    colSums(
      otu_rarefied,
      na.rm = TRUE
    ) > 0,
    drop = FALSE
  ]

  common_samples <- intersect(
    rownames(otu_sample),
    rownames(otu_rarefied)
  )

  otu_original <- otu_sample[
    common_samples,
    ,
    drop = FALSE
  ]

  otu_rarefied <- otu_rarefied[
    common_samples,
    ,
    drop = FALSE
  ]

  all_asvs <- union(
    colnames(otu_original),
    colnames(otu_rarefied)
  )

  otu_original <- fill_missing_cols(
    otu_original,
    all_asvs
  )

  otu_rarefied <- fill_missing_cols(
    otu_rarefied,
    all_asvs
  )

  ########## Shannon robustness ##########
  original_shannon <- vegan::diversity(
    otu_original,
    index = "shannon"
  )

  rarefied_shannon <- vegan::diversity(
    otu_rarefied,
    index = "shannon"
  )

  shannon_pearson <- cor.test(
    original_shannon,
    rarefied_shannon,
    method = "pearson"
  )

  shannon_spearman <- cor.test(
    original_shannon,
    rarefied_shannon,
    method = "spearman"
  )

  shannon_result <- data.frame(
    Dataset = dataset_name,
    Rarefaction_depth = rarefaction_depth,
    Pearson_r = as.numeric(
      shannon_pearson$estimate
    ),
    Pearson_P = shannon_pearson$p.value,
    Spearman_rho = as.numeric(
      shannon_spearman$estimate
    ),
    Spearman_P = shannon_spearman$p.value
  )

  cat(
    "\nShannon robustness:\n"
  )

  print(
    shannon_result
  )

  shannon_label <- paste0(
    "Pearson r = ",
    round(
      as.numeric(
        shannon_pearson$estimate
      ),
      3
    ),
    "\n",
    format_p(
      shannon_pearson$p.value
    ),
    "\n",
    "Spearman \u03c1 = ",
    round(
      as.numeric(
        shannon_spearman$estimate
      ),
      3
    ),
    "\n",
    format_p(
      shannon_spearman$p.value
    )
  )

  p_shannon <- plot_concordance_density(
    plot_df = data.frame(
      x = original_shannon,
      y = rarefied_shannon
    ),
    x_lab = "Original Shannon diversity",
    y_lab = "Rarefied Shannon diversity",
    label_text = shannon_label,
    point_color = point_color,
    fixed_range = NULL,
    breaks = NULL,
    clamp_prediction = FALSE,
    main_to_density_ratio = main_to_density_ratio
  )

  ########## Bray-Curtis robustness ##########
  bc_original <- vegan::vegdist(
    otu_original,
    method = "bray"
  )

  bc_rarefied <- vegan::vegdist(
    otu_rarefied,
    method = "bray"
  )

  mantel_result <- vegan::mantel(
    bc_original,
    bc_rarefied,
    method = "spearman",
    permutations = mantel_permutations
  )

  mantel_result_table <- data.frame(
    Dataset = dataset_name,
    Rarefaction_depth = rarefaction_depth,
    Mantel_method = "Spearman",
    Mantel_r = as.numeric(
      mantel_result$statistic
    ),
    Mantel_P = mantel_result$signif,
    Permutations = mantel_permutations
  )

  cat(
    "\nBray-Curtis robustness:\n"
  )

  print(
    mantel_result_table
  )

  bc_label <- paste0(
    "Mantel r = ",
    round(
      as.numeric(
        mantel_result$statistic
      ),
      3
    ),
    "\n",
    format_p(
      mantel_result$signif
    )
  )

  p_bc <- plot_concordance_density(
    plot_df = data.frame(
      x = as.vector(
        bc_original
      ),
      y = as.vector(
        bc_rarefied
      )
    ),
    x_lab = "Original Rotifer Community Composition\n(Bray-Curtis Dissimilarity)",
    y_lab = "Rarefied Rotifer Community Composition\n(Bray-Curtis Dissimilarity)",
    label_text = bc_label,
    point_color = point_color,
    fixed_range = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.25
    ),
    clamp_prediction = TRUE,
    main_to_density_ratio = main_to_density_ratio
  )

  list(
    Shannon_result = shannon_result,
    BrayCurtis_result = mantel_result_table,
    Shannon_plot = p_shannon,
    BrayCurtis_plot = p_bc
  )
}

########## COI ##########
coi_result <- run_rarefaction_robustness(
  input_file = "COI_rotifer_unrarefied.csv",
  dataset_name = "COI",
  rarefaction_depth = 53,
  point_color = "#F87A52"
)

########## 18S rRNA ##########
rrna_result <- run_rarefaction_robustness(
  input_file = "18S_rotifer_unrarefied.csv",
  dataset_name = "18S rRNA",
  rarefaction_depth = 29,
  point_color = "#5EAFCF"
)

########## Display plots ##########
print(
  coi_result$Shannon_plot
)

print(
  coi_result$BrayCurtis_plot
)

print(
  rrna_result$Shannon_plot
)

print(
  rrna_result$BrayCurtis_plot
)
