rm(list = ls())

########## Load required packages ##########
library(betapart)
library(vegan)
library(permute)
library(dplyr)
library(ggplot2)

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
analysis_seed <- 20260719

marker_colours <- c(
  "COI" = "#F87A52",
  "18S rRNA" = "#5EAFCF"
)

environment_vars <- c(
  "pH",
  "SWC",
  "NO3_N",
  "NH4_N",
  "TC",
  "TP",
  "AP",
  "MBN",
  "AcP",
  "LAP"
)

########## Prepare metadata ##########
metadata <- metadata %>%
  transmute(
    Sample = as.character(`OTU ID`),
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

########## Retain common samples ##########
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

metadata_model <- metadata %>%
  filter(
    Sample %in% common_samples
  ) %>%
  arrange(
    match(
      Sample,
      common_samples
    )
  ) %>%
  filter(
    complete.cases(
      across(
        all_of(environment_vars)
      )
    )
  )

final_samples <- metadata_model$Sample

cat(
  "Common samples used:",
  length(final_samples),
  "\n"
)

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

coi_comm <- prepare_community(
  coi_data,
  final_samples
)

rrna_comm <- prepare_community(
  rrna_data,
  final_samples
)

########## Standardize environmental variables ##########
environment <- metadata_model %>%
  select(
    all_of(environment_vars)
  ) %>%
  as.data.frame()

rownames(environment) <- final_samples

environment_scaled <- as.data.frame(
  scale(environment),
  check.names = FALSE
)

rownames(environment_scaled) <- final_samples

########## Incidence-based Jaccard decomposition ##########
decompose_jaccard <- function(community) {

  community_pa <- ifelse(
    community > 0,
    1,
    0
  )

  beta_pair <- betapart::beta.pair(
    community_pa,
    index.family = "jaccard"
  )

  list(
    Turnover = beta_pair$beta.jtu,
    Nestedness = beta_pair$beta.jne
  )
}

coi_beta <- decompose_jaccard(
  coi_comm
)

rrna_beta <- decompose_jaccard(
  rrna_comm
)

########## Generate one common permutation matrix ##########
set.seed(analysis_seed)

permutation_matrix <- shuffleSet(
  n = nrow(environment_scaled),
  nset = nperm,
  control = how(
    nperm = nperm
  )
)

########## Fit environmental dbRDA models ##########
fit_component_model <- function(
    response_dist,
    marker_name,
    component_name
) {

  model <- dbrda(
    response_dist ~ .,
    data = environment_scaled,
    sqrt.dist = FALSE,
    add = "lingoes"
  )

  overall_test <- anova(
    model,
    permutations = permutation_matrix,
    model = "reduced"
  )

  r2_result <- RsquareAdj(
    model
  )

  tibble(
    Marker = marker_name,
    Component = component_name,
    R2 = unname(
      r2_result$r.squared
    ),
    Adjusted_R2 = unname(
      r2_result$adj.r.squared
    ),
    Adjusted_R2_percent =
      unname(
        r2_result$adj.r.squared
      ) * 100,
    P_value = as.data.frame(
      overall_test,
      check.names = FALSE
    )[1, "Pr(>F)"]
  )
}

fig4f_data <- bind_rows(
  fit_component_model(
    coi_beta$Turnover,
    "COI",
    "Turnover"
  ),
  fit_component_model(
    coi_beta$Nestedness,
    "COI",
    "Nestedness-resultant"
  ),
  fit_component_model(
    rrna_beta$Turnover,
    "18S rRNA",
    "Turnover"
  ),
  fit_component_model(
    rrna_beta$Nestedness,
    "18S rRNA",
    "Nestedness-resultant"
  )
) %>%
  mutate(
    Marker = factor(
      Marker,
      levels = c(
        "COI",
        "18S rRNA"
      )
    ),
    Component = factor(
      Component,
      levels = c(
        "Turnover",
        "Nestedness-resultant"
      )
    ),
    P_label = case_when(
      is.na(P_value) ~ "P = NA",
      P_value < 0.001 ~ "P < 0.001",
      TRUE ~ paste0(
        "P = ",
        sprintf(
          "%.3f",
          P_value
        )
      )
    ),
    Label_y = ifelse(
      Adjusted_R2_percent >= 0,
      Adjusted_R2_percent + 0.25,
      Adjusted_R2_percent - 0.25
    )
  )

########## Show results ##########
cat("\nFig. 4F environmental dbRDA results:\n")
print(fig4f_data)

########## Figure 4F ##########
fig4f <- ggplot(
  fig4f_data,
  aes(
    x = Component,
    y = Adjusted_R2_percent,
    fill = Marker
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_col(
    position = position_dodge(
      width = 0.75
    ),
    width = 0.68
  ) +
  geom_text(
    aes(
      y = Label_y,
      label = P_label
    ),
    position = position_dodge(
      width = 0.75
    ),
    size = 3.3
  ) +
  scale_fill_manual(
    values = marker_colours,
    breaks = c(
      "COI",
      "18S rRNA"
    ),
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      7.5,
      by = 2.5
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.03
      )
    )
  ) +
  coord_cartesian(
    ylim = c(
      -2.0,
      8.1
    ),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = expression(
      paste(
        "Adjusted ",
        R^2,
        " (%)"
      )
    ),
    fill = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12
    ),
    legend.key.size = grid::unit(
      0.8,
      "cm"
    ),
    axis.title.y = element_text(
      size = 12
    ),
    axis.text.y = element_text(
      size = 11
    ),
    axis.text.x = element_text(
      size = 11,
      angle = 20,
      hjust = 1
    ),
    axis.line = element_line(
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      linewidth = 0.6
    ),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

print(fig4f)
