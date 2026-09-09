rm(list = ls())

########## Load required packages ##########
library(psych)
library(dplyr)
library(pheatmap)

########## Read data ##########
coi_data <- read.csv(
  "COI_data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rrna_data <- read.csv(
  "18S_data.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  "metadata_all.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

########## Settings ##########
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

# All physicochemical variables used in Fig. S7.
# Elevation and soil layer are not included in this correlation heatmap.
environment_columns <- c(
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

environment_labels <- c(
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

coi_top8 <- c(
  "Adineta",
  "Macrotrachela",
  "Philodina",
  "Polyarthra",
  "Pleuretra",
  "Habrotrocha",
  "Lecane",
  "Anomopus"
)

rrna_top8 <- c(
  "Rotaria",
  "Habrotrocha",
  "Anomopus",
  "Philodina",
  "Lecane",
  "Floscularia",
  "Bryceella",
  "Mytilina"
)

heatmap_breaks <- seq(
  -0.45,
  0.45,
  length.out = 100
)

heatmap_colours <- colorRampPalette(
  c(
    "#4979b8",
    "white",
    "#e1302e"
  )
)(99)

########## Prepare environmental matrix ##########
metadata_env <- data.frame(
  Sample = as.character(
    metadata[["OTU ID"]]
  ),
  check.names = FALSE
)

for (internal_name in names(environment_columns)) {

  source_name <- environment_columns[[internal_name]]

  metadata_env[[internal_name]] <- as.numeric(
    metadata[[source_name]]
  )
}

########## Helper functions ##########
clean_genus <- function(x) {

  x <- trimws(
    as.character(x)
  )

  x <- sub(
    "^g__",
    "",
    x
  )

  x[
    is.na(x) |
      x == "" |
      x == "NA" |
      x == "Unassigned"
  ] <- NA_character_

  x
}

prepare_genus_matrix <- function(
    data,
    target_genera,
    marker_name
) {

  sample_cols <- setdiff(
    colnames(data),
    tax_cols
  )

  abundance <- data[
    ,
    sample_cols,
    drop = FALSE
  ]

  abundance[] <- lapply(
    abundance,
    as.numeric
  )

  genus_abundance <- data.frame(
    Genus = clean_genus(
      data$Genus
    ),
    abundance,
    check.names = FALSE
  ) %>%
    filter(
      !is.na(Genus)
    ) %>%
    group_by(
      Genus
    ) %>%
    summarise(
      across(
        all_of(sample_cols),
        ~ sum(
          .x,
          na.rm = TRUE
        )
      ),
      .groups = "drop"
    )

  missing_genera <- setdiff(
    target_genera,
    genus_abundance$Genus
  )

  if (length(missing_genera) > 0) {
    stop(
      paste0(
        marker_name,
        ": target genera not found: ",
        paste(
          missing_genera,
          collapse = ", "
        )
      )
    )
  }

  genus_abundance <- genus_abundance %>%
    filter(
      Genus %in% target_genera
    ) %>%
    mutate(
      Genus = factor(
        Genus,
        levels = target_genera
      )
    ) %>%
    arrange(
      Genus
    )

  genus_matrix <- as.matrix(
    genus_abundance[
      ,
      sample_cols,
      drop = FALSE
    ]
  )

  rownames(genus_matrix) <- as.character(
    genus_abundance$Genus
  )

  storage.mode(genus_matrix) <- "numeric"

  genus_matrix
}

run_genus_environment_correlation <- function(
    data,
    target_genera,
    marker_name
) {

  genus_matrix <- prepare_genus_matrix(
    data = data,
    target_genera = target_genera,
    marker_name = marker_name
  )

  common_samples <- intersect(
    colnames(genus_matrix),
    metadata_env$Sample
  )

  genus_sample <- t(
    genus_matrix[
      ,
      common_samples,
      drop = FALSE
    ]
  )

  environment_sample <- metadata_env[
    match(
      common_samples,
      metadata_env$Sample
    ),
    names(environment_columns),
    drop = FALSE
  ]

  rownames(environment_sample) <- common_samples

  complete_samples <- complete.cases(
    environment_sample
  )

  genus_sample <- genus_sample[
    complete_samples,
    ,
    drop = FALSE
  ]

  environment_sample <- environment_sample[
    complete_samples,
    ,
    drop = FALSE
  ]

  correlation <- psych::corr.test(
    genus_sample,
    environment_sample,
    method = "spearman",
    adjust = "none"
  )

  correlation_matrix <- correlation$r
  p_matrix <- correlation$p

  colnames(correlation_matrix) <- unname(
    environment_labels[
      colnames(correlation_matrix)
    ]
  )

  colnames(p_matrix) <- unname(
    environment_labels[
      colnames(p_matrix)
    ]
  )

  significance <- matrix(
    "",
    nrow = nrow(p_matrix),
    ncol = ncol(p_matrix),
    dimnames = dimnames(p_matrix)
  )

  significance[
    p_matrix < 0.01
  ] <- "**"

  significance[
    p_matrix >= 0.01 &
      p_matrix < 0.05
  ] <- "*"

  heatmap <- pheatmap::pheatmap(
    correlation_matrix,
    scale = "none",
    color = heatmap_colours,
    breaks = heatmap_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "white",
    display_numbers = significance,
    fontsize_number = 12,
    number_color = "white",
    cellwidth = 20,
    cellheight = 20,
    main = paste0(
      marker_name,
      ": Top 8 genera vs environmental factors"
    ),
    silent = TRUE
  )

  list(
    correlation = correlation_matrix,
    p_value = p_matrix,
    significance = significance,
    heatmap = heatmap
  )
}

########## COI ##########
coi_result <- run_genus_environment_correlation(
  data = coi_data,
  target_genera = coi_top8,
  marker_name = "COI"
)

########## 18S rRNA ##########
rrna_result <- run_genus_environment_correlation(
  data = rrna_data,
  target_genera = rrna_top8,
  marker_name = "18S rRNA"
)

########## Show results ##########
cat(
  "\nCOI Spearman correlations:\n"
)
print(
  round(
    coi_result$correlation,
    3
  )
)

cat(
  "\nCOI raw P values:\n"
)
print(
  round(
    coi_result$p_value,
    4
  )
)

cat(
  "\n18S rRNA Spearman correlations:\n"
)
print(
  round(
    rrna_result$correlation,
    3
  )
)

cat(
  "\n18S rRNA raw P values:\n"
)
print(
  round(
    rrna_result$p_value,
    4
  )
)

########## Display Fig. S7 heatmaps ##########
grid::grid.newpage()
grid::grid.draw(
  coi_result$heatmap$gtable
)

grid::grid.newpage()
grid::grid.draw(
  rrna_result$heatmap$gtable
)
