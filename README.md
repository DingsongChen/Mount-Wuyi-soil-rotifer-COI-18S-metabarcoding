# COI and 18S rRNA gene metabarcoding of soil rotifers in Wuyishan National Park

This repository contains the R scripts and input data used to reproduce the main and supplementary analyses for the manuscript:

**Comparison of COI and 18S rRNA metabarcoding for characterizing the taxonomic profiles and diversity of soil rotifers**

The study compares COI and 18S rRNA gene metabarcoding for characterizing soil rotifer communities along elevational and soil-depth gradients in **Wuyishan National Park, southeastern China**. Samples were collected across seven elevations (300, 600, 900, 1200, 1600, 1800, and 2100 m) from topsoil (0–10 cm) and subsoil (10–20 cm).

## Repository structure

The repository is organized by figure and supplementary table. Each analysis folder contains the corresponding R script(s) and the input file(s) required to reproduce that analysis.

```text
.
├── README.md
├── Main_Figures/
│   ├── Fig_1/
│   ├── Fig_2/
│   ├── Fig_3/
│   └── Fig_4/
├── Supplementary_Figures/
│   ├── Fig_S3/
│   ├── Fig_S5/
│   ├── Fig_S6/
│   └── Fig_S7/
└── Supplementary_Tables/
    ├── Table_S4/
    └── Table_S5/
```

### Main figures

- **Fig. 1** — Taxonomic profiles detected by COI and 18S rRNA gene metabarcoding, including taxonomic richness, taxonomic overlap, and comparison with morphological observations.
- **Fig. 2** — Dominant genera, elevational distributions, and alpha diversity.
- **Fig. 3** — Community differentiation, PERMANOVA, distance-decay relationships, and beta-diversity partitioning.
- **Fig. 4** — Environmental associations, hierarchical partitioning, co-inertia analysis, and environmental explanation of beta-diversity components.

### Supplementary figures

- **Fig. S3** — Rarefaction robustness analysis comparing alpha diversity and Bray-Curtis dissimilarity before and after rarefaction.
- **Fig. S5** — ASV richness comparisons between markers, soil layers, and elevations.
- **Fig. S6** — Significant linear and quadratic relationships between alpha diversity or morphological abundance and environmental variables.
- **Fig. S7** — Spearman correlations between the dominant genera and physicochemical variables.

### Supplementary tables

- **Table S4** — Type III repeated-measures ANOVA for Shannon diversity and ASV richness.
- **Table S5** — Common-depth rarefaction sensitivity analysis using 41 matched samples, with both markers repeatedly rarefied to 29 reads.
- **Table S6** — PERMANOVA results are reproduced by the analysis code provided with **Fig. 3** and therefore do not require a separate script.

## Main input files

Several input files are reused across analyses.

- `COI_data.csv` — Final rarefied COI rotifer ASV table with taxonomic assignments.
- `18S_data.csv` — Final rarefied 18S rRNA gene rotifer ASV table with taxonomic assignments.
- `metadata.csv` — Sample metadata used for the main community and diversity analyses.
- `metadata_all.csv` — Expanded metadata table containing the full set of physicochemical variables; used only for Fig. S6 and Fig. S7.
- `COI_rotifer_unrarefied.csv` — Filtered but non-rarefied COI rotifer ASV abundance matrix.
- `18S_rotifer_unrarefied.csv` — Filtered but non-rarefied 18S rRNA gene rotifer ASV abundance matrix.
- `Number.csv` — Morphological abundance data used in Fig. S6.
- `genus heatmap.csv` — Genus-level input used for the elevational heatmap in Fig. 2.

Because each figure/table folder contains its own required input files, individual analyses can be run independently.

## Statistical analyses

The repository includes the analyses reported in the manuscript, including:

- Shannon diversity and ASV richness
- paired Wilcoxon tests
- Kruskal-Wallis tests and Dunn post hoc comparisons
- Type III repeated-measures ANOVA
- Bray-Curtis dissimilarity
- PERMANOVA
- Mantel tests
- PCoA
- SIMPER
- incidence-based Jaccard beta-diversity partitioning
- distance-decay relationships
- db-RDA
- hierarchical partitioning
- co-inertia analysis
- Spearman correlation analysis
- linear and quadratic regression
- repeated-rarefaction sensitivity analyses

For analyses involving multiple testing or permutations, the settings implemented in the corresponding R scripts reproduce those used in the manuscript.

## Software

Analyses were conducted in **R 4.5.0**.

The scripts use packages including:

```text
vegan
ggplot2
dplyr
tidyr
afex
rstatix
ggpubr
betapart
ade4
rdacca.hp
psych
pheatmap
circlize
cowplot
```

Additional packages required for individual figures are loaded within the corresponding scripts.

## Running the analyses

1. Download or clone this repository.
2. Open the folder for the figure or supplementary table you want to reproduce.
3. Set the R working directory to that folder, if needed.
4. Run the corresponding `.R` script.

The scripts use relative file names rather than computer-specific absolute paths. Input files should therefore remain in the same folder as the script unless the file paths are modified.

Most scripts display figures and statistical results directly in the R session rather than automatically saving output files.

## Sequencing data

Raw sequencing data are available from the NCBI Sequence Read Archive under BioProject:

**PRJNA1244230**

SRA accessions:

- `SRR39229812–SRR39229881`
- `SRR39230836–SRR39230905`

## Study design

Soil samples were collected in July 2024 along an elevational gradient in Wuyishan National Park. Seven elevations were sampled:

**300, 600, 900, 1200, 1600, 1800, and 2100 m**

At each elevation, five 20 × 20 m plots were established. Topsoil (0–10 cm) and subsoil (10–20 cm) were sampled separately.

Metabarcoding was conducted using:

- **COI:** mlCOIintF / jgHCO2198
- **18S rRNA gene:** TAReuk454FWD / REV3

Sequence processing was performed using QIIME 2 and DADA2, and taxonomic assignments were obtained against the NCBI-NT database.

## Notes on reproducibility

The repository contains the analysis-ready data used for the figures and supplementary tables. Some scripts use filtered, non-rarefied ASV matrices for rarefaction sensitivity analyses, whereas the primary diversity and community analyses use the final rarefied datasets.

For the common-depth sensitivity analysis in Table S5, the 41 samples shared by both markers were repeatedly rarefied to 29 reads for 999 iterations.

## Citation

If you use these data or scripts, please cite the associated manuscript once published.

A full citation will be added after publication.

## Contact

For questions regarding the data or analysis workflow, please open an issue in this repository.
