rm(list = ls())

########## Load required packages ##########
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

########## Read data ##########
coi_data <- read.csv("COI_data.csv", check.names = FALSE)
rrna_data <- read.csv("18S_data.csv", check.names = FALSE)

########## Set taxonomic ranks ##########
rank_levels <- c("Class", "Order", "Family", "Genus", "Species")

########## Extract unique taxa ##########
get_taxa <- function(data, rank) {
  taxa <- trimws(as.character(data[[rank]]))
  taxa <- taxa[!is.na(taxa) & taxa != ""]
  sort(unique(taxa))
}

########## Compare shared and marker-specific taxa ##########
overlap_summary <- bind_rows(
  lapply(rank_levels, function(rank_name) {
    coi_set <- get_taxa(coi_data, rank_name)
    rrna_set <- get_taxa(rrna_data, rank_name)
    
    data.frame(
      Rank = rank_name,
      COI_only = length(setdiff(coi_set, rrna_set)),
      Shared = length(intersect(coi_set, rrna_set)),
      S18_only = length(setdiff(rrna_set, coi_set))
    )
  })
)

overlap_long <- overlap_summary %>%
  pivot_longer(
    cols = c(COI_only, Shared, S18_only),
    names_to = "Category",
    values_to = "N_taxa"
  ) %>%
  mutate(
    Category = recode(
      Category,
      "COI_only" = "COI-only",
      "S18_only" = "18S-only"
    ),
    Rank = factor(Rank, levels = rank_levels),
    Category = factor(
      Category,
      levels = c("COI-only", "Shared", "18S-only")
    )
  ) %>%
  group_by(Rank) %>%
  mutate(
    Proportion = N_taxa / sum(N_taxa),
    Label = ifelse(Proportion < 0.03, "", paste0(round(Proportion * 100), "%"))
  ) %>%
  ungroup()

########## Show results ##########
print(overlap_summary)
print(overlap_long)

########## Plot settings ##########
overlap_cols <- c(
  "COI-only" = "#F87A52",
  "Shared" = "grey75",
  "18S-only" = "#5EAFCF"
)

common_theme <- theme_bw(base_size = 15) +
  theme(
    text = element_text(family = "serif", color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 13, color = "black"),
    axis.text.x = element_text(
      size = 13,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    aspect.ratio = 1
  )

########## Plot Figure 1D ##########
p1d <- ggplot(
  overlap_long,
  aes(
    x = Rank,
    y = Proportion,
    fill = Category
  )
) +
  geom_col(
    width = 0.68,
    color = "black",
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = Label),
    position = position_stack(vjust = 0.5),
    size = 5,
    family = "serif",
    color = "black"
  ) +
  scale_fill_manual(
    values = overlap_cols,
    breaks = c("COI-only", "Shared", "18S-only")
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Taxonomic rank",
    y = "Proportion of taxa",
    fill = NULL
  ) +
  common_theme

print(p1d)