rm(list = ls())

########## Load required packages ##########
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)

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
  "metadata_all.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

number_data <- read.csv(
  "Number.csv",
  row.names = 1,
  check.names = FALSE
)

########## Settings ##########
p_threshold <- 0.05

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

environment_columns <- c(
  "Elevation" = "Elevation",
  "pH" = "pH",
  "SWC" = "SWC",
  "NO3_N" = "NO3--N",
  "NH4_N" = "NH4+-N",
  "TN" = "TN",
  "DON" = "DON",
  "TC" = "TC",
  "DOC" = "DOC",
  "TP" = "TP",
  "AP" = "AP",
  "MBN" = "MBN",
  "MBC" = "MBC",
  "MBP" = "MBP",
  "AcP" = "AcP",
  "BG" = "BG",
  "NAG" = "NAG",
  "LAP" = "LAP"
)

display_names <- c(
  "Elevation" = "Elevation",
  "pH" = "pH",
  "SWC" = "SWC",
  "NO3_N" = "NO3-N",
  "NH4_N" = "NH4-N",
  "TN" = "TN",
  "DON" = "DON",
  "TC" = "TC",
  "DOC" = "DOC",
  "TP" = "TP",
  "AP" = "AP",
  "MBN" = "MBN",
  "MBC" = "MBC",
  "MBP" = "MBP",
  "AcP" = "AcP",
  "BG" = "BG",
  "NAG" = "NAG",
  "LAP" = "LAP"
)

########## Prepare metadata ##########
metadata_clean <- data.frame(
  Sample = as.character(
    metadata[["OTU ID"]]
  ),
  check.names = FALSE
)

for (internal_name in names(environment_columns)) {

  source_name <- environment_columns[[internal_name]]

  metadata_clean[[internal_name]] <- as.numeric(
    metadata[[source_name]]
  )
}

########## Calculate Shannon diversity ##########
calculate_shannon <- function(
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

  shannon <- vegan::diversity(
    t(abundance),
    index = "shannon"
  )

  tibble(
    Sample = sample_cols,
    Response = as.numeric(shannon),
    Dataset = marker_name
  ) %>%
    filter(
      is.finite(Response),
      Response > 0
    ) %>%
    left_join(
      metadata_clean,
      by = "Sample"
    )
}

coi_shannon <- calculate_shannon(
  coi_data,
  "COI Shannon"
)

rrna_shannon <- calculate_shannon(
  rrna_data,
  "18S rRNA Shannon"
)

########## Prepare morphological abundance ##########
number_df <- number_data %>%
  rownames_to_column(
    "Sample"
  )

if (!"Number" %in% colnames(number_df)) {
  stop(
    "Number.csv must contain a column named 'Number'."
  )
}

number_env <- number_df %>%
  transmute(
    Sample = as.character(Sample),
    Response = as.numeric(Number),
    Dataset = "Morphological abundance"
  ) %>%
  left_join(
    metadata_clean,
    by = "Sample"
  ) %>%
  filter(
    is.finite(Response)
  )

########## Regression helper ##########
fit_environment_models <- function(
    data,
    dataset_name,
    point_color,
    line_color,
    fill_color,
    response_label,
    force_zero = FALSE
) {

  linear_results <- list()
  quadratic_results <- list()
  linear_plots <- list()
  quadratic_plots <- list()

  response_max <- max(
    data$Response,
    na.rm = TRUE
  )

  if (
    grepl(
      "Shannon",
      dataset_name,
      fixed = TRUE
    )
  ) {
    y_max <- ceiling(
      response_max / 0.5
    ) * 0.5

    y_breaks <- seq(
      0,
      y_max,
      by = 0.5
    )
  }

  for (env in names(environment_columns)) {

    sub_data <- data %>%
      filter(
        complete.cases(
          Response,
          .data[[env]]
        )
      )

    if (nrow(sub_data) < 4) {
      next
    }

    x <- sub_data[[env]]
    y <- sub_data$Response

    linear_fit <- lm(
      y ~ x
    )

    linear_summary <- summary(
      linear_fit
    )

    linear_p <- linear_summary$coefficients[
      2,
      4
    ]

    linear_r2 <- linear_summary$r.squared

    quadratic_fit <- lm(
      y ~ poly(
        x,
        2
      )
    )

    quadratic_summary <- summary(
      quadratic_fit
    )

    quadratic_f <- quadratic_summary$fstatistic

    quadratic_p <- pf(
      quadratic_f[1],
      quadratic_f[2],
      quadratic_f[3],
      lower.tail = FALSE
    )

    quadratic_r2 <- quadratic_summary$r.squared

    env_label <- unname(
      display_names[[env]]
    )

    if (
      is.finite(linear_p) &&
      linear_p < p_threshold
    ) {

      linear_results[[env]] <- tibble(
        Dataset = dataset_name,
        Environment = env_label,
        Model = "Linear",
        n = nrow(sub_data),
        R2 = linear_r2,
        P = linear_p
      )

      p_linear <- ggplot(
        sub_data,
        aes(
          x = .data[[env]],
          y = Response
        )
      ) +
        geom_point(
          color = point_color,
          alpha = 0.7,
          size = 2
        ) +
        geom_smooth(
          method = "lm",
          formula = y ~ x,
          se = TRUE,
          color = line_color,
          fill = fill_color
        ) +
        labs(
          title = paste(
            dataset_name,
            "-",
            env_label
          ),
          subtitle = paste0(
            "R² = ",
            round(
              linear_r2,
              3
            ),
            ", p = ",
            format(
              round(
                linear_p,
                4
              ),
              nsmall = 4
            )
          ),
          x = env_label,
          y = response_label
        ) +
        theme_classic() +
        theme(
          panel.border = element_rect(
            color = "black",
            fill = NA
          ),
          panel.grid = element_blank(),
          axis.title = element_text(
            size = 20
          ),
          axis.text = element_text(
            size = 15
          )
        )

      if (
        grepl(
          "Shannon",
          dataset_name,
          fixed = TRUE
        )
      ) {
        p_linear <- p_linear +
          coord_cartesian(
            ylim = c(
              0,
              y_max
            )
          ) +
          scale_y_continuous(
            breaks = y_breaks
          )
      } else if (force_zero) {
        p_linear <- p_linear +
          expand_limits(
            y = 0
          )
      }

      linear_plots[[env]] <- p_linear
    }

    if (
      is.finite(quadratic_p) &&
      quadratic_p < p_threshold
    ) {

      quadratic_results[[env]] <- tibble(
        Dataset = dataset_name,
        Environment = env_label,
        Model = "Quadratic",
        n = nrow(sub_data),
        R2 = quadratic_r2,
        P = quadratic_p
      )

      p_quadratic <- ggplot(
        sub_data,
        aes(
          x = .data[[env]],
          y = Response
        )
      ) +
        geom_point(
          color = point_color,
          alpha = 0.7,
          size = 2
        ) +
        geom_smooth(
          method = "lm",
          formula = y ~ poly(
            x,
            2
          ),
          se = TRUE,
          color = line_color,
          fill = fill_color
        ) +
        labs(
          title = paste(
            dataset_name,
            "-",
            env_label
          ),
          subtitle = paste0(
            "R² = ",
            round(
              quadratic_r2,
              3
            ),
            ", p = ",
            format(
              round(
                quadratic_p,
                4
              ),
              nsmall = 4
            )
          ),
          x = env_label,
          y = response_label
        ) +
        theme_classic() +
        theme(
          panel.border = element_rect(
            color = "black",
            fill = NA
          ),
          panel.grid = element_blank(),
          axis.title = element_text(
            size = 20
          ),
          axis.text = element_text(
            size = 15
          )
        )

      if (
        grepl(
          "Shannon",
          dataset_name,
          fixed = TRUE
        )
      ) {
        p_quadratic <- p_quadratic +
          coord_cartesian(
            ylim = c(
              0,
              y_max
            )
          ) +
          scale_y_continuous(
            breaks = y_breaks
          )
      } else if (force_zero) {
        p_quadratic <- p_quadratic +
          expand_limits(
            y = 0
          )
      }

      quadratic_plots[[env]] <- p_quadratic
    }
  }

  list(
    Linear_results = bind_rows(
      linear_results
    ),
    Quadratic_results = bind_rows(
      quadratic_results
    ),
    Linear_plots = linear_plots,
    Quadratic_plots = quadratic_plots
  )
}

########## Run analyses ##########
coi_results <- fit_environment_models(
  data = coi_shannon,
  dataset_name = "COI Shannon",
  point_color = "#F87A52",
  line_color = "#F87A52",
  fill_color = "#F87A524D",
  response_label = "Shannon Index"
)

rrna_results <- fit_environment_models(
  data = rrna_shannon,
  dataset_name = "18S rRNA Shannon",
  point_color = "#5EAFCF",
  line_color = "#5EAFCF",
  fill_color = "#5EAFCF4D",
  response_label = "Shannon Index"
)

number_results <- fit_environment_models(
  data = number_env,
  dataset_name = "Morphological abundance",
  point_color = "#EB7FFF",
  line_color = "#EB7FFF",
  fill_color = "#F9D9FF",
  response_label = "Number",
  force_zero = TRUE
)

########## Print significant result tables ##########
print_result_group <- function(
    result,
    label
) {

  cat(
    "\n============================================================\n",
    label,
    "\n============================================================\n",
    sep = ""
  )

  cat(
    "\nLinear relationships (P < 0.05):\n"
  )

  if (
    nrow(
      result$Linear_results
    ) > 0
  ) {
    print(
      result$Linear_results
    )
  } else {
    cat(
      "None\n"
    )
  }

  cat(
    "\nQuadratic relationships (P < 0.05):\n"
  )

  if (
    nrow(
      result$Quadratic_results
    ) > 0
  ) {
    print(
      result$Quadratic_results
    )
  } else {
    cat(
      "None\n"
    )
  }
}

print_result_group(
  coi_results,
  "COI SHANNON"
)

print_result_group(
  rrna_results,
  "18S rRNA SHANNON"
)

print_result_group(
  number_results,
  "MORPHOLOGICAL ABUNDANCE"
)

########## Print significant plots by category ##########
print_plot_group <- function(
    result,
    label
) {

  cat(
    "\n============================================================\n",
    label,
    " - LINEAR PLOTS\n",
    "============================================================\n",
    sep = ""
  )

  if (
    length(
      result$Linear_plots
    ) == 0
  ) {
    cat(
      "No significant linear plots.\n"
    )
  } else {
    for (
      plot_name in names(
        result$Linear_plots
      )
    ) {
      cat(
        "\n",
        plot_name,
        "\n",
        sep = ""
      )

      print(
        result$Linear_plots[[plot_name]]
      )
    }
  }

  cat(
    "\n============================================================\n",
    label,
    " - QUADRATIC PLOTS\n",
    "============================================================\n",
    sep = ""
  )

  if (
    length(
      result$Quadratic_plots
    ) == 0
  ) {
    cat(
      "No significant quadratic plots.\n"
    )
  } else {
    for (
      plot_name in names(
        result$Quadratic_plots
      )
    ) {
      cat(
        "\n",
        plot_name,
        "\n",
        sep = ""
      )

      print(
        result$Quadratic_plots[[plot_name]]
      )
    }
  }
}

print_plot_group(
  coi_results,
  "COI SHANNON"
)

print_plot_group(
  rrna_results,
  "18S rRNA SHANNON"
)

print_plot_group(
  number_results,
  "MORPHOLOGICAL ABUNDANCE"
)
