rm(list = ls())

########## Load required packages ##########
library(ggplot2)

########## Read data ##########
coi_data <- read.csv("COI_data.csv", check.names = FALSE)
rrna_data <- read.csv("18S_data.csv", check.names = FALSE)

########## Count detected taxa ##########
rank_levels <- c("Class", "Order", "Family", "Genus", "Species")

count_taxa <- function(data, rank) {
  taxa <- data[[rank]]
  taxa <- taxa[!is.na(taxa) & taxa != ""]
  length(unique(taxa))
}

taxa_counts <- rbind(
  data.frame(
    Rank = rank_levels,
    Marker = "COI",
    Detected_taxa = sapply(rank_levels, function(x) count_taxa(coi_data, x))
  ),
  data.frame(
    Rank = rank_levels,
    Marker = "18S",
    Detected_taxa = sapply(rank_levels, function(x) count_taxa(rrna_data, x))
  )
)

taxa_counts$Rank <- factor(
  taxa_counts$Rank,
  levels = rank_levels
)

taxa_counts$Marker <- factor(
  taxa_counts$Marker,
  levels = c("COI", "18S")
)

########## Show results ##########
print(taxa_counts)

########## Plot settings ##########
marker_cols <- c(
  "COI" = "#F87A52",
  "18S" = "#5EAFCF"
)

y_max <- max(taxa_counts$Detected_taxa, na.rm = TRUE)
y_upper <- ceiling(y_max * 1.12)

if (is.na(y_upper) || y_upper <= 0) {
  y_upper <- 1
}

if (y_upper == y_max) {
  y_upper <- y_max + 1
}

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
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    ),
    aspect.ratio = 1
  )

########## Plot Figure 1C ##########
p1c <- ggplot(
  taxa_counts,
  aes(
    x = Rank,
    y = Detected_taxa,
    group = Marker,
    color = Marker
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3.2) +
  scale_color_manual(values = marker_cols) +
  scale_y_continuous(
    limits = c(0, y_upper),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Taxonomic rank",
    y = "Number of detected taxa",
    color = NULL
  ) +
  common_theme

print(p1c)