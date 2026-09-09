rm(list = ls())

########## Load required packages ##########
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)

########## Enter morphology and marker detection data ##########
morph_detect <- tibble::tribble(
  ~Taxon,                           ~Morphology, ~COI, ~`18S rRNA`,
  "Adineta vaga",                   1,           1,    0,
  "Macrotrachela quadricornifera",  1,           1,    0,
  "Habrotrocha sp.",                1,           1,    1,
  "Lepadella sp.",                  1,           0,    1
)

########## Prepare plotting data ##########
plot_data <- morph_detect %>%
  pivot_longer(
    cols = c(Morphology, COI, `18S rRNA`),
    names_to = "Method",
    values_to = "Detected"
  ) %>%
  mutate(
    Method = factor(
      Method,
      levels = c("Morphology", "COI", "18S rRNA")
    ),
    Taxon = factor(
      Taxon,
      levels = c(
        "Lepadella sp.",
        "Habrotrocha sp.",
        "Adineta vaga",
        "Macrotrachela quadricornifera"
      )
    ),
    Fill_group = ifelse(
      Detected == 1,
      as.character(Method),
      "Not detected"
    ),
    Fill_group = factor(
      Fill_group,
      levels = c(
        "Morphology",
        "COI",
        "18S rRNA",
        "Not detected"
      )
    )
  )

########## Plot settings ##########
method_cols <- c(
  "Morphology" = "#E76BF3",
  "COI" = "#FC7854",
  "18S rRNA" = "#5AB1CD",
  "Not detected" = "grey92"
)

taxon_labels <- c(
  "Adineta vaga" = "<i>Adineta vaga</i>",
  "Macrotrachela quadricornifera" =
    "<i>Macrotrachela quadricornifera</i>",
  "Habrotrocha sp." = "<i>Habrotrocha</i> sp.",
  "Lepadella sp." = "<i>Lepadella</i> sp."
)

########## Plot Figure 1E ##########
p1e <- ggplot(
  plot_data,
  aes(x = Taxon, y = Method, fill = Fill_group)
) +
  geom_tile(
    color = "white",
    linewidth = 0.9,
    width = 0.92,
    height = 0.92
  ) +
  scale_fill_manual(
    values = method_cols,
    name = NULL
  ) +
  scale_x_discrete(labels = taxon_labels) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    axis.text.x = element_markdown(
      size = 9,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 10,
      color = "black"
    ),
    axis.ticks = element_blank(),
    plot.title = element_text(
      size = 11,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 9,
      hjust = 0
    ),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

print(p1e)