rm(list = ls())

########## Load required packages ##########
library(vegan)
library(dplyr)
library(tidyr)
library(afex)

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
  "metadata.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

########## Type III ANOVA settings ##########
options(
  contrasts = c(
    "contr.sum",
    "contr.poly"
  )
)

########## Calculate alpha diversity ##########
calculate_alpha <- function(
    data,
    marker_name,
    metadata
) {

  sample_cols <- intersect(
    as.character(
      metadata[["OTU ID"]]
    ),
    colnames(data)
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
    Marker = marker_name,
    Shannon = vegan::diversity(
      t(abundance),
      index = "shannon"
    ),
    ASV_richness = colSums(
      abundance > 0
    )
  )
}

alpha_data <- bind_rows(
  calculate_alpha(
    coi_data,
    "COI",
    metadata
  ),
  calculate_alpha(
    rrna_data,
    "18S rRNA",
    metadata
  )
) %>%
  left_join(
    metadata %>%
      transmute(
        Sample = as.character(
          `OTU ID`
        ),
        Layer = recode(
          as.character(Group),
          "T" = "Topsoil",
          "S" = "Subsoil"
        ),
        Elevation = Elevation
      ),
    by = "Sample"
  ) %>%
  mutate(
    PlotID = sub(
      "^[TS]",
      "",
      Sample
    )
  )

########## Retain complete plots ##########
complete_plots <- alpha_data %>%
  distinct(
    PlotID,
    Marker,
    Layer
  ) %>%
  group_by(
    PlotID
  ) %>%
  summarise(
    n_marker = n_distinct(
      Marker
    ),
    n_layer = n_distinct(
      Layer
    ),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(
    n_marker == 2,
    n_layer == 2,
    n_cells == 4
  ) %>%
  pull(
    PlotID
  )

alpha_complete <- alpha_data %>%
  filter(
    PlotID %in% complete_plots
  ) %>%
  mutate(
    PlotID = factor(
      PlotID
    ),
    Marker = factor(
      Marker,
      levels = c(
        "COI",
        "18S rRNA"
      )
    ),
    Layer = factor(
      Layer,
      levels = c(
        "Topsoil",
        "Subsoil"
      )
    ),
    Elevation = factor(
      Elevation,
      levels = sort(
        unique(
          Elevation
        )
      )
    )
  ) %>%
  arrange(
    Elevation,
    PlotID,
    Layer,
    Marker
  )

if (
  n_distinct(
    alpha_complete$PlotID
  ) != 20
) {
  stop(
    "The analysis must contain 20 plots with complete observations for both markers and both soil layers."
  )
}

cat(
  "\nComplete plots used in the repeated-measures ANOVA: ",
  n_distinct(
    alpha_complete$PlotID
  ),
  "\n",
  sep = ""
)

cat(
  "Observations per alpha-diversity index: ",
  nrow(
    alpha_complete
  ),
  "\n",
  sep = ""
)

########## Convert to long format ##########
alpha_long <- alpha_complete %>%
  pivot_longer(
    cols = c(
      Shannon,
      ASV_richness
    ),
    names_to = "Index",
    values_to = "Value"
  )

########## Type III repeated-measures ANOVA ##########
run_rm_anova <- function(
    data,
    index_name
) {

  dat <- data %>%
    filter(
      Index == index_name
    ) %>%
    droplevels()

  afex::aov_ez(
    id = "PlotID",
    dv = "Value",
    data = dat,
    within = c(
      "Marker",
      "Layer"
    ),
    between = "Elevation",
    type = 3,
    anova_table = list(
      correction = "none",
      es = "pes"
    )
  )
}

shannon_anova <- run_rm_anova(
  alpha_long,
  "Shannon"
)

richness_anova <- run_rm_anova(
  alpha_long,
  "ASV_richness"
)

########## Format ANOVA tables ##########
effect_labels <- c(
  "Marker" = "Marker",
  "Elevation" = "Elevation",
  "Layer" = "Soil layer",
  "Marker:Elevation" = "Marker × Elevation",
  "Marker:Layer" = "Marker × Soil layer",
  "Elevation:Layer" = "Elevation × Soil layer",
  "Marker:Elevation:Layer" = "Marker × Elevation × Soil layer"
)

format_anova_table <- function(
    model,
    index_label
) {

  tab <- as.data.frame(
    model$anova_table,
    check.names = FALSE
  )

  tab$Effect <- rownames(
    tab
  )

  rownames(
    tab
  ) <- NULL

  out <- tibble(
    Index = index_label,
    Effect = recode(
      tab$Effect,
      !!!effect_labels
    ),
    Num_Df = tab[["num Df"]],
    Den_Df = tab[["den Df"]],
    F_value = tab[["F"]],
    P_value = tab[["Pr(>F)"]],
    Partial_eta2 = tab[["pes"]]
  )

  out
}

shannon_results <- format_anova_table(
  shannon_anova,
  "Shannon diversity"
)

richness_results <- format_anova_table(
  richness_anova,
  "ASV richness"
)

anova_results <- bind_rows(
  shannon_results,
  richness_results
)

########## Display results ##########
cat(
  "\n============================================================\n",
  "SHANNON DIVERSITY\n",
  "============================================================\n",
  sep = ""
)

print(
  shannon_results,
  n = Inf
)

cat(
  "\n============================================================\n",
  "ASV RICHNESS\n",
  "============================================================\n",
  sep = ""
)

print(
  richness_results,
  n = Inf
)

cat(
  "\n============================================================\n",
  "COMBINED TABLE\n",
  "============================================================\n",
  sep = ""
)

print(
  anova_results,
  n = Inf
)

########## Key effect sizes reported in the manuscript ##########
key_effects <- anova_results %>%
  filter(
    Effect %in% c(
      "Marker",
      "Soil layer"
    )
  ) %>%
  select(
    Index,
    Effect,
    P_value,
    Partial_eta2
  )

cat(
  "\n============================================================\n",
  "KEY MAIN EFFECTS\n",
  "============================================================\n",
  sep = ""
)

print(
  key_effects,
  n = Inf
)
