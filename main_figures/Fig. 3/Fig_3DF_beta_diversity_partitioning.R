rm(list = ls())

########## Load required packages ##########
library(betapart)
library(ggplot2)
library(ggtern)
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

########## Prepare sample columns ##########
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

common_samples <- intersect(
  coi_samples,
  rrna_samples
)

cat(
  "COI samples:",
  length(coi_samples),
  "\n"
)

cat(
  "18S rRNA samples:",
  length(rrna_samples),
  "\n"
)

cat(
  "Shared samples used for both markers:",
  length(common_samples),
  "\n"
)

########## Prepare ASV matrices ##########
prepare_asv_table <- function(
    data,
    sample_names
) {

  asv_table <- as.matrix(
    data[
      ,
      sample_names,
      drop = FALSE
    ]
  )

  storage.mode(asv_table) <- "numeric"

  asv_table[
    is.na(asv_table)
  ] <- 0

  asv_table <- asv_table[
    rowSums(
      asv_table,
      na.rm = TRUE
    ) > 0,
    ,
    drop = FALSE
  ]

  asv_table
}

coi_asv <- prepare_asv_table(
  coi_data,
  common_samples
)

rrna_asv <- prepare_asv_table(
  rrna_data,
  common_samples
)

########## Jaccard beta-diversity decomposition ##########
decompose_jaccard <- function(
    asv_table,
    marker_name
) {

  community <- t(
    asv_table
  )

  community <- community[
    rowSums(
      community,
      na.rm = TRUE
    ) > 0,
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

  community_pa <- ifelse(
    community > 0,
    1,
    0
  )

  beta_core <- betapart::betapart.core(
    community_pa
  )

  beta_pair <- betapart::beta.pair(
    beta_core,
    index.family = "jac"
  )

  pairwise_data <- data.frame(
    Similarity =
      1 - as.vector(
        beta_pair$beta.jac
      ),
    Turnover =
      as.vector(
        beta_pair$beta.jtu
      ),
    Nestedness =
      as.vector(
        beta_pair$beta.jne
      ),
    Marker = marker_name
  )

  pairwise_data
}

coi_beta <- decompose_jaccard(
  coi_asv,
  "COI"
)

rrna_beta <- decompose_jaccard(
  rrna_asv,
  "18S rRNA"
)

########## Show summary ##########
beta_summary <- bind_rows(
  coi_beta,
  rrna_beta
) %>%
  group_by(Marker) %>%
  summarise(
    Similarity = mean(
      Similarity,
      na.rm = TRUE
    ),
    Turnover = mean(
      Turnover,
      na.rm = TRUE
    ),
    Nestedness = mean(
      Nestedness,
      na.rm = TRUE
    ),
    Pairwise_comparisons = n(),
    .groups = "drop"
  )

print(beta_summary)

########## Plot settings ##########
marker_colors <- c(
  "COI" = "#F87A52",
  "18S rRNA" = "#5EAFCF"
)

plot_ternary <- function(
    data,
    marker_name
) {

  ggtern(
    data = data,
    aes(
      x = Similarity,
      y = Turnover,
      z = Nestedness
    )
  ) +
    geom_point(
      color = marker_colors[
        marker_name
      ],
      shape = 16,
      size = 2.0,
      alpha = 0.45
    ) +
    labs(
      x = expression(
        "Similarity (" * 1 - beta[jac] * ")"
      ),
      y = expression(
        "Turnover (" * beta[jtu] * ")"
      ),
      z = expression(
        "Nestedness (" * beta[jne] * ")"
      ),
      title = marker_name
    ) +
    theme_rgbw() +
    theme(
      text = element_text(
        family = "serif",
        color = "black"
      ),
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      legend.position = "none"
    )
}

########## Figure 3D: COI ##########
p_coi <- plot_ternary(
  coi_beta,
  "COI"
)

########## Figure 3F: 18S rRNA ##########
p_rrna <- plot_ternary(
  rrna_beta,
  "18S rRNA"
)

########## Display plots ##########
print(p_coi)
print(p_rrna)
