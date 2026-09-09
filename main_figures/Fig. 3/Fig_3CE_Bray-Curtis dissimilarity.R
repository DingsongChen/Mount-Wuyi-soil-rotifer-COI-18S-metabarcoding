rm(list = ls())

########## Load required packages ##########
library(ggplot2)
library(vegan)

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
  metadata$`OTU ID`
)

########## Prepare community matrix ##########
prepare_community <- function(data, metadata) {
  
  sample_cols <- intersect(
    metadata$Sample,
    colnames(data)
  )
  
  community <- t(
    as.matrix(
      data[, sample_cols, drop = FALSE]
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
  
  metadata_marker <- metadata[
    match(
      rownames(community),
      metadata$Sample
    ),
    ,
    drop = FALSE
  ]
  
  if (
    !identical(
      rownames(community),
      metadata_marker$Sample
    )
  ) {
    stop(
      "Sample order does not match between the community matrix and metadata."
    )
  }
  
  list(
    community = community,
    metadata = metadata_marker
  )
}

coi <- prepare_community(
  coi_data,
  metadata
)

rrna <- prepare_community(
  rrna_data,
  metadata
)

########## Mantel analysis ##########
run_distance_decay <- function(
    community,
    metadata,
    marker_name
) {
  
  bray_dist <- vegdist(
    community,
    method = "bray"
  )
  
  elevation <- as.numeric(
    as.character(
      metadata$Elevation
    )
  )
  
  elevation_dist <- dist(
    elevation,
    method = "euclidean"
  )
  
  mantel_result <- mantel(
    bray_dist,
    elevation_dist,
    method = "spearman",
    permutations = 9999,
    na.rm = TRUE
  )
  
  plot_data <- data.frame(
    Bray = as.vector(
      bray_dist
    ),
    Elevational_distance = as.vector(
      elevation_dist
    )
  )
  
  cat(
    "\n",
    marker_name,
    " Mantel test\n",
    "R = ",
    round(
      mantel_result$statistic,
      4
    ),
    "\nP = ",
    mantel_result$signif,
    "\n",
    sep = ""
  )
  
  list(
    mantel = mantel_result,
    data = plot_data
  )
}

coi_result <- run_distance_decay(
  coi$community,
  coi$metadata,
  "COI"
)

rrna_result <- run_distance_decay(
  rrna$community,
  rrna$metadata,
  "18S rRNA"
)

########## Plotting function ##########
plot_distance_decay <- function(
    result,
    marker_name,
    low_color,
    high_color
) {
  
  plot_data <- result$data
  
  mantel_r <- round(
    result$mantel$statistic,
    4
  )
  
  mantel_p <- result$mantel$signif
  
  p_label <- ifelse(
    mantel_p < 0.001,
    "P < 0.001",
    paste0(
      "P = ",
      format(
        round(
          mantel_p,
          3
        ),
        nsmall = 3
      )
    )
  )
  
  ggplot(
    plot_data,
    aes(
      x = Elevational_distance,
      y = Bray
    )
  ) +
    geom_point(
      aes(
        fill = Elevational_distance / 1000
      ),
      size = 4,
      alpha = 0.75,
      colour = "black",
      shape = 21
    ) +
    geom_smooth(
      method = "lm",
      colour = "black",
      alpha = 0.2,
      level = 0.95
    ) +
    scale_fill_continuous(
      low = low_color,
      high = high_color,
      name = "Elevational gradient"
    ) +
    labs(
      x = "Elevational distance",
      y = "Bray-Curtis Dissimilarity",
      title = marker_name
    ) +
    annotate(
      "text",
      label = paste0(
        "R = ",
        mantel_r
      ),
      x = max(
        plot_data$Elevational_distance
      ) * 0.9,
      y = max(
        plot_data$Bray
      ) * 0.9,
      size = 5,
      colour = "black"
    ) +
    annotate(
      "text",
      label = p_label,
      x = max(
        plot_data$Elevational_distance
      ) * 0.9,
      y = max(
        plot_data$Bray
      ) * 0.85,
      size = 5,
      colour = "black"
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,
      axis.text = element_text(
        colour = "black",
        size = 12
      ),
      axis.title = element_text(
        face = "bold",
        size = 14,
        colour = "black"
      ),
      panel.grid = element_blank(),
      panel.border = element_rect(
        fill = NA,
        colour = "black"
      ),
      legend.position = "top",
      legend.text = element_text(
        size = 10,
        face = "bold"
      ),
      legend.title = element_text(
        size = 11,
        face = "bold"
      )
    )
}

########## Figure 3C: COI ##########
p_coi <- plot_distance_decay(
  result = coi_result,
  marker_name = "COI",
  low_color = "#F64007",
  high_color = "#F9CB42"
)

########## Figure 3E: 18S rRNA ##########
p_rrna <- plot_distance_decay(
  result = rrna_result,
  marker_name = "18S rRNA",
  low_color = "#163D7B",
  high_color = "#5EAFCF"
)

########## Display plots ##########
print(p_coi)
print(p_rrna)