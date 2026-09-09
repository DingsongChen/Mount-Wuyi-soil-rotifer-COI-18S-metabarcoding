rm(list = ls())
graphics.off()

########## Load required packages ##########
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(VennDiagram)

########## Read data ##########
data <- read.csv(
  "genus heatmap.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Keep compatibility with the earlier input file
data$Group[data$Group == "Bdelloidae"] <- "Bdelloidea"

########## Calculate relative abundance ##########
elevation_cols <- c(
  "300 m", "600 m", "900 m", "1200 m",
  "1600 m", "1800 m", "2100 m"
)

expr_rel_by_cluster <- data %>%
  group_by(Cluster) %>%
  mutate(
    across(
      all_of(elevation_cols),
      ~ (.x / sum(.x)) * 100
    )
  ) %>%
  ungroup()

########## Prepare heatmap matrix ##########
expr_matrix <- as.matrix(
  expr_rel_by_cluster[, elevation_cols]
)

rownames(expr_matrix) <- expr_rel_by_cluster$Label

# Row-wise Z-score standardization
expr_matrix <- t(
  scale(
    t(expr_matrix)
  )
)

expr_rel_by_cluster$Cluster <- factor(
  expr_rel_by_cluster$Cluster,
  levels = c(
    "Soil-18S rRNA",
    "Soil-COI"
  )
)

split_factor <- expr_rel_by_cluster$Cluster

########## Plot settings ##########
cluster_colors <- c(
  "Soil-18S rRNA" = "#C88688",
  "Soil-COI" = "#FCF8F9"
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c(
    "#3A84B5",
    "white",
    "#C54C43"
  )
)

group_shapes <- list(
  "Bdelloidea" = 21,
  "Monogononta" = 23
)

group_colors <- list(
  "Bdelloidea" = "#efca72",
  "Monogononta" = "#93cc82"
)

########## Open plotting device ##########
windows(
  width = 10,
  height = 10
)

########## Initialize circular heatmap ##########
circos.clear()

circos.par(
  start.degree = 60,
  gap.after = c(5, 30),
  track.margin = c(0, 0.01),
  cell.padding = c(0, 0, 0, 0)
)

########## Draw circular heatmap ##########
circos.heatmap(
  expr_matrix,
  split = split_factor,
  cluster = FALSE,
  bg.border = "grey20",
  bg.lwd = 1,
  cell.border = "white",
  cell.lwd = 0.5,
  col = col_fun,
  track.height = 0.3,
  rownames.side = "outside",
  rownames.cex = 0.9
)

########## Add elevation labels ##########
circos.track(
  track.index = get.current.track.index(),
  panel.fun = function(x, y) {
    
    if (CELL_META$sector.numeric.index == 2) {
      
      cn <- colnames(expr_matrix)
      n <- length(cn)
      
      circos.text(
        rep(
          CELL_META$cell.xlim[2],
          n
        ) + convert_x(1, "mm"),
        1:n - 0.5,
        cn,
        cex = 0.8,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  },
  bg.border = NA
)

########## Add rotifer-group track ##########
circos.track(
  ylim = c(0, 1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y) {
    
    current_sector <- CELL_META$sector.index
    
    sector_idx <- which(
      expr_rel_by_cluster$Cluster == current_sector
    )
    
    sector_groups <- expr_rel_by_cluster$Group[
      sector_idx
    ]
    
    n_points <- length(sector_groups)
    
    for (i in seq_len(n_points)) {
      
      x_pos <-
        CELL_META$xlim[1] +
        (
          CELL_META$xlim[2] -
            CELL_META$xlim[1]
        ) *
        (i - 0.5) / n_points
      
      current_group <- as.character(
        sector_groups[i]
      )
      
      circos.points(
        x = x_pos,
        y = 0.5,
        pch = group_shapes[[current_group]],
        cex = 1.3,
        col = "black",
        bg = group_colors[[current_group]]
      )
    }
    
    if (
      CELL_META$sector.numeric.index ==
      length(levels(split_factor))
    ) {
      
      circos.lines(
        c(
          CELL_META$cell.xlim[2],
          CELL_META$cell.xlim[2] +
            convert_x(1, "mm")
        ),
        c(0.5, 0.5),
        col = "black",
        lwd = 2
      )
      
      circos.text(
        CELL_META$cell.xlim[2] +
          convert_x(1.5, "mm"),
        0.5,
        "Order",
        cex = 0.9,
        adj = c(0, 0.5),
        facing = "inside"
      )
    }
  }
)

########## Add marker labels ##########
circos.track(
  ylim = c(0, 1),
  track.height = 0.05,
  bg.col = adjustcolor(
    cluster_colors[
      levels(split_factor)
    ],
    alpha.f = 0.8
  ),
  panel.fun = function(x, y) {
    
    current_sector <- CELL_META$sector.index
    
    if (
      current_sector == "Soil-18S rRNA"
    ) {
      
      text_y <- CELL_META$ylim[1] + 0.45
      text_facing <- "bending.outside"
      
    } else if (
      current_sector == "Soil-COI"
    ) {
      
      text_y <- CELL_META$ylim[2] - 0.5
      text_facing <- "bending.inside"
    }
    
    circos.text(
      CELL_META$xcenter,
      text_y,
      current_sector,
      facing = text_facing,
      cex = 0.9,
      adj = c(0.5, 0.5)
    )
  }
)

########## Prepare central Venn diagram ##########
cluster_sets <- split(
  data$Label,
  data$Cluster
)

cluster_sets <- cluster_sets[
  c(
    "Soil-18S rRNA",
    "Soil-COI"
  )
]

fill_colors <- cluster_colors[
  names(cluster_sets)
]

venn_plot <- venn.diagram(
  x = cluster_sets,
  filename = NULL,
  fill = fill_colors,
  alpha = 0.5,
  cex = 1,
  cat.cex = 0,
  cat.dist = 0.05,
  margin = 0.1,
  lty = "dashed",
  lwd = 1
)

########## Add central Venn diagram ##########
grid::pushViewport(
  grid::viewport(
    x = 0.5,
    y = 0.5,
    width = 0.3,
    height = 0.3,
    just = c(
      "center",
      "center"
    )
  )
)

grid::grid.draw(
  venn_plot
)

grid::popViewport()

grid::grid.text(
  "Genus annotation",
  x = 0.5,
  y = 0.42,
  gp = grid::gpar(
    fontsize = 13,
    fontface = "bold"
  )
)

########## Add legends ##########
heatmap_legend <- Legend(
  title = "Relative abundance",
  col_fun = col_fun,
  at = seq(
    -2,
    2,
    length.out = 5
  ),
  title_position = "leftcenter-rot",
  title_gp = grid::gpar(
    fontsize = 12
  ),
  labels_gp = grid::gpar(
    fontsize = 12
  )
)

group_legend <- Legend(
  title = "Order",
  title_position = "leftcenter-rot",
  title_gp = grid::gpar(
    fontsize = 12
  ),
  labels_gp = grid::gpar(
    fontsize = 12
  ),
  nrow = 2,
  labels = c(
    "Bdelloidea",
    "Monogononta"
  ),
  type = "points",
  pch = c(21, 23),
  size = grid::unit(
    4,
    "mm"
  ),
  legend_gp = grid::gpar(
    col = "black",
    fill = c(
      "#efca72",
      "#93cc82"
    )
  )
)

draw(
  heatmap_legend,
  x = grid::unit(
    0.9,
    "npc"
  ) -
    grid::unit(
      5,
      "mm"
    ),
  y = grid::unit(
    0.95,
    "npc"
  ) -
    grid::unit(
      5,
      "mm"
    ),
  just = c(
    "right",
    "top"
  )
)

draw(
  group_legend,
  x = grid::unit(
    0.95,
    "npc"
  ) -
    grid::unit(
      5,
      "mm"
    ),
  y = grid::unit(
    0.15,
    "npc"
  ) -
    grid::unit(
      5,
      "mm"
    ),
  just = c(
    "right",
    "top"
  )
)

########## Reset circular plotting parameters ##########
circos.clear()