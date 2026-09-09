rm(list = ls())

########## Load required package ##########
library(vegan)

########## Settings ##########
coi_file <- "COI_rotifer_unrarefied.csv"
rrna_file <- "18S_rotifer_unrarefied.csv"

common_depth <- 29
n_iter <- 999
set.seed(123)

########## Read filtered, non-rarefied ASV matrices ##########
read_asv_matrix <- function(file) {

  x <- read.csv(
    file,
    row.names = 1,
    check.names = FALSE
  )

  x[] <- lapply(
    x,
    as.numeric
  )

  x <- as.matrix(x)
  storage.mode(x) <- "numeric"

  if (anyNA(x)) {
    stop(
      paste0(
        file,
        " contains non-numeric values."
      )
    )
  }

  if (any(x < 0)) {
    stop(
      paste0(
        file,
        " contains negative abundances."
      )
    )
  }

  x
}

coi_mat <- read_asv_matrix(
  coi_file
)

rrna_mat <- read_asv_matrix(
  rrna_file
)

########## Match samples ##########
common_samples <- intersect(
  colnames(coi_mat),
  colnames(rrna_mat)
)

if (length(common_samples) != 41) {
  stop(
    paste0(
      "Expected 41 matched samples, but found ",
      length(common_samples),
      "."
    )
  )
}

coi_mat <- coi_mat[
  ,
  common_samples,
  drop = FALSE
]

rrna_mat <- rrna_mat[
  ,
  common_samples,
  drop = FALSE
]

depth_check <- data.frame(
  Sample = common_samples,
  COI_reads = colSums(coi_mat),
  `18S_rRNA_reads` = colSums(rrna_mat),
  check.names = FALSE
)

if (
  any(
    depth_check$COI_reads < common_depth |
      depth_check$`18S_rRNA_reads` < common_depth
  )
) {
  stop(
    "At least one matched sample has fewer than 29 reads."
  )
}

cat(
  "\nMatched samples: ",
  length(common_samples),
  "\n",
  sep = ""
)

cat(
  "Common rarefaction depth: ",
  common_depth,
  " reads\n",
  sep = ""
)

cat(
  "Repeated rarefactions: ",
  n_iter,
  "\n",
  sep = ""
)

########## Rarefaction and alpha diversity ##########
calculate_alpha_once <- function(
    asv_mat,
    depth
) {

  rarefied <- vegan::rrarefy(
    t(asv_mat),
    sample = depth
  )

  data.frame(
    Sample = rownames(rarefied),
    Shannon = vegan::diversity(
      rarefied,
      index = "shannon"
    ),
    ASV_richness = rowSums(
      rarefied > 0
    ),
    row.names = NULL
  )
}

n_samples <- length(
  common_samples
)

coi_shannon <- matrix(
  NA_real_,
  nrow = n_samples,
  ncol = n_iter,
  dimnames = list(
    common_samples,
    NULL
  )
)

rrna_shannon <- coi_shannon
coi_richness <- coi_shannon
rrna_richness <- coi_shannon

for (i in seq_len(n_iter)) {

  coi_alpha <- calculate_alpha_once(
    coi_mat,
    common_depth
  )

  rrna_alpha <- calculate_alpha_once(
    rrna_mat,
    common_depth
  )

  coi_alpha <- coi_alpha[
    match(
      common_samples,
      coi_alpha$Sample
    ),
    ,
    drop = FALSE
  ]

  rrna_alpha <- rrna_alpha[
    match(
      common_samples,
      rrna_alpha$Sample
    ),
    ,
    drop = FALSE
  ]

  coi_shannon[
    ,
    i
  ] <- coi_alpha$Shannon

  rrna_shannon[
    ,
    i
  ] <- rrna_alpha$Shannon

  coi_richness[
    ,
    i
  ] <- coi_alpha$ASV_richness

  rrna_richness[
    ,
    i
  ] <- rrna_alpha$ASV_richness
}

########## Mean alpha diversity across 999 rarefactions ##########
sample_mean_results <- data.frame(
  Sample = common_samples,
  COI_Shannon = rowMeans(
    coi_shannon
  ),
  `18S_rRNA_Shannon` = rowMeans(
    rrna_shannon
  ),
  COI_ASV_richness = rowMeans(
    coi_richness
  ),
  `18S_rRNA_ASV_richness` = rowMeans(
    rrna_richness
  ),
  check.names = FALSE
)

########## Paired Wilcoxon tests ##########
shannon_test <- wilcox.test(
  sample_mean_results$`18S_rRNA_Shannon`,
  sample_mean_results$COI_Shannon,
  paired = TRUE,
  exact = FALSE
)

richness_test <- wilcox.test(
  sample_mean_results$`18S_rRNA_ASV_richness`,
  sample_mean_results$COI_ASV_richness,
  paired = TRUE,
  exact = FALSE
)

########## Iteration-level marker differences ##########
shannon_iteration_difference <- colMeans(
  rrna_shannon - coi_shannon
)

richness_iteration_difference <- colMeans(
  rrna_richness - coi_richness
)

########## Summary table ##########
make_summary <- function(
    metric,
    coi_values,
    rrna_values,
    iteration_difference,
    test
) {

  data.frame(
    Metric = metric,
    COI_mean = mean(
      coi_values
    ),
    `18S_rRNA_mean` = mean(
      rrna_values
    ),
    Mean_difference_18S_minus_COI = mean(
      rrna_values - coi_values
    ),
    Difference_2.5_percentile = unname(
      quantile(
        iteration_difference,
        0.025
      )
    ),
    Difference_97.5_percentile = unname(
      quantile(
        iteration_difference,
        0.975
      )
    ),
    Proportion_iterations_18S_greater_COI = mean(
      iteration_difference > 0
    ),
    Wilcoxon_V = unname(
      test$statistic
    ),
    Wilcoxon_P = test$p.value,
    check.names = FALSE
  )
}

table_s5 <- rbind(
  make_summary(
    metric = "Shannon diversity",
    coi_values = sample_mean_results$COI_Shannon,
    rrna_values = sample_mean_results$`18S_rRNA_Shannon`,
    iteration_difference = shannon_iteration_difference,
    test = shannon_test
  ),
  make_summary(
    metric = "ASV richness",
    coi_values = sample_mean_results$COI_ASV_richness,
    rrna_values = sample_mean_results$`18S_rRNA_ASV_richness`,
    iteration_difference = richness_iteration_difference,
    test = richness_test
  )
)

########## Display results ##########
cat(
  "\n============================================================\n",
  "TABLE S5: COMMON-DEPTH SENSITIVITY ANALYSIS\n",
  "============================================================\n",
  sep = ""
)

print(
  table_s5,
  row.names = FALSE
)

cat(
  "\nSample-level mean alpha diversity is stored in: sample_mean_results\n"
)
