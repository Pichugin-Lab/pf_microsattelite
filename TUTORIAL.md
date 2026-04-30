# Repository Tour & Tutorial: *P. falciparum* Microsatellite Analysis

This tutorial walks you through every analysis and figure produced by
`microsattelite.qmd`. By the end you will be able to reproduce all
results from scratch on your own machine.

---

## Table of Contents

1. [Repository Overview](#1-repository-overview)
2. [Prerequisites](#2-prerequisites)
3. [Data](#3-data)
4. [Running the Analysis](#4-running-the-analysis)
5. [Step-by-Step Walkthrough](#5-step-by-step-walkthrough)
   - [5.1 Import & Clean Data](#51-import--clean-data)
   - [5.2 Convert to `genind` Object](#52-convert-to-genind-object)
   - [5.3 Population Diversity Statistics](#53-population-diversity-statistics)
   - [5.4 Locus-Specific Descriptors](#54-locus-specific-descriptors)
   - [5.5 Minimum Spanning Network](#55-minimum-spanning-network)
   - [5.6 Genotype Accumulation Curve](#56-genotype-accumulation-curve)
   - [5.7 Discriminant Analysis of Principal Components (DAPC)](#57-discriminant-analysis-of-principal-components-dapc)
   - [5.8 Final PCA Plot](#58-final-pca-plot)
   - [5.9 Population Statistics with Confidence Intervals](#59-population-statistics-with-confidence-intervals)
   - [5.10 Genetic Distance Matrices](#510-genetic-distance-matrices)
   - [5.11 Package Citations](#511-package-citations)
6. [Output Files & Figures](#6-output-files--figures)
7. [Troubleshooting](#7-troubleshooting)

---

## 1. Repository Overview

```
pf_microsattelite/
├── data/
│   └── 20240726_ms_data.csv   # Raw microsatellite genotype data
├── figures/                   # Generated figures (PNG/SVG)
│   ├── PCA.png
│   ├── PCA_narrow.png
│   ├── min-spanning-network.png
│   ├── population_barplot.png
│   ├── table_loci-unique-alleles.png
│   └── table_microsattelilte-diverity.png
├── microsattelite.qmd         # Main Quarto analysis document
├── TUTORIAL.md                # This file
└── README.md
```

The primary analysis document is `microsattelite.qmd`. It is a
[Quarto](https://quarto.org/) notebook written in R that performs
population-genetics analyses on *Plasmodium falciparum* microsatellite
data collected from Kenya, Ghana, Colombia, and Thailand/Cambodia. All
figures and summary tables are generated inline and several are saved to
`figures/`.

---

## 2. Prerequisites

### Software

| Software | Minimum version | Install |
|---|---|---|
| R | ≥ 4.3 | https://cran.r-project.org/ |
| Quarto | ≥ 1.4 | https://quarto.org/docs/get-started/ |
| RStudio (optional) | ≥ 2023.12 | https://posit.co/download/rstudio-desktop/ |

### R Packages

Install all required packages from an R session before rendering:

```r
install.packages(c(
  "tidyverse",   # data wrangling & plotting
  "ggplot2",     # plotting
  "ggrepel",     # non-overlapping labels
  "ggside",      # marginal histograms
  "ggforce",     # ellipse annotations
  "gt",          # formatted tables
  "dplyr",
  "tidyr",
  "scales",
  "writexl"      # Excel export
))

# Bioconductor / population-genetics packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages(c("ape", "pegas"))   # phylogenetics helpers

# adegenet and poppr are on CRAN but may need Bioconductor dependencies
install.packages(c("adegenet", "poppr"))

install.packages("igraph")            # network visualisation
```

> **Tip:** If you are on macOS you may need to install Xcode Command
> Line Tools (`xcode-select --install`) and the `gfortran` compiler for
> some packages to build from source.

---

## 3. Data

### File: `data/20240726_ms_data.csv`

The CSV contains one row per sample and the following columns:

| Column | Description |
|---|---|
| `Sample` | Unique sample identifier. Prefix encodes country of origin (K = Kenya, G = Ghana, C = Colombia, T = Thailand). Reference strains (3D7, 7G8, Dd2, FVO(NIH), NF54) are also present. |
| `Poly-A` | Polyclonality flag. Rows marked `POLYCLONAL` are excluded from downstream analysis. Rows marked `RERUN` are treated as missing. |
| `PFPK2`, `TAA81`, `ARA2`, `TAA87`, `TA40`, `TA42`, `2490`, `TA1`, `TAA60`, `TAA109`, `Pfg377` | Allele sizes (in base-pairs) for 11 microsatellite loci typed by fragment analysis. |
| `AMA1 Sequenced`, `MSP1 Sequenced`, `MSP2 Sequenced`, `CSP Fullgene Sequenced` | Flags indicating whether additional amplicon sequencing was performed (not used in this analysis). |
| `Additional Notes` | Free-text field (dropped during import). |

**Sample naming convention**

- A bare alphanumeric ID (e.g. `K1`, `G3`) indicates the **original
  vial** (uncloned parasite stock).
- An ID containing a hyphen (e.g. `K2-1`, `K2-2`) indicates a **cloned
  line** derived from the original stock.
- Samples with the word *"cloned"* in their name are also treated as
  clones.

The analysis uses 228 rows (including the header), yielding ~227
parasite isolates before quality filtering.

---

## 4. Running the Analysis

### Option A — Render the entire Quarto document (recommended)

```bash
# From the repository root
quarto render microsattelite.qmd
```

This produces `microsattelite.html` with all figures embedded inline,
and also writes PNG/SVG files to `figures/`.

> **Note:** Two code chunks are marked `#| eval: false` (the interactive
> `find.clusters()` calls). These require manual input in an interactive
> R session to choose the number of clusters (see §5.7). You will need
> to run those chunks by hand and then store the resulting `clust`
> object before the rest of the DAPC section can execute.

### Option B — Run interactively in RStudio

1. Open `microsattelite.qmd` in RStudio.
2. Click **Run All** (or use Ctrl+Alt+R / Cmd+Option+R) to execute all
   chunks sequentially.
3. For the `find.clusters()` chunks, follow the prompts in the R
   console to select the number of PCs and clusters (see §5.7 for
   recommended values).

---

## 5. Step-by-Step Walkthrough

### 5.1 Import & Clean Data

**Goal:** Load the raw CSV, recode missing/invalid values, and attach
metadata (country of origin, passage history).

```r
library(tidyverse)

# Read raw data
msdf <- read_csv('data/20240726_ms_data.csv')

# Treat "RERUN" entries as true NA
msdf <- msdf |>
  mutate(across(where(is.character), ~ na_if(., "RERUN")))

# Derive country and passage_history from the Sample column
msdf_metadata <- msdf |>
  mutate(
    country = case_when(
      str_detect(str_to_lower(Sample), "3d7|7g8|fvo|nf54|dd2") ~ "reference",
      str_starts(Sample, "K") ~ "Kenya",
      str_starts(Sample, "G") ~ "Ghana",
      str_starts(Sample, "C") ~ "Colombia",
      str_starts(Sample, "T") ~ "Thailand",
      TRUE ~ "unknown"
    ),
    passage_history = case_when(
      str_detect(str_to_lower(Sample), "cloned") ~ "clone",
      str_detect(str_to_lower(Sample), "original vial") ~ "original",
      str_detect(Sample, "-") ~ "clone",
      TRUE ~ "original"
    )
  )

# Quality filter: remove NA samples and polyclonal infections
msdf_filter <- msdf |>
  drop_na(Sample) |>
  filter(`Poly-A` != "POLYCLONAL") |>
  select(-`AMA1 Sequenced`, -`MSP1 Sequenced`,
         -`MSP2 Sequenced`, -`CSP Fullgene Sequenced`,
         -`Additional Notes`)

# Join country metadata back onto the filtered data
msdf_metadata_filtered <- inner_join(msdf_metadata, msdf_filter, by = "Sample")
msdf_filter <- msdf_filter |>
  left_join(select(msdf_metadata_filtered, Sample, country), by = "Sample")
```

**Key data frames produced:**

| Data frame | Description |
|---|---|
| `msdf_metadata` | Full dataset with country & passage_history columns added |
| `msdf_metadata_filtered` | Filtered dataset joined with metadata |
| `msdf_filter` | Filtered dataset with country column — used for all downstream analysis |

---

### 5.2 Convert to `genind` Object

**Goal:** Package the allele data into an `adegenet` `genind` object so
that population-genetics functions from `adegenet` and `poppr` can be
applied.

```r
library(adegenet)
library(poppr)

pfms <- df2genind(
  X        = msdf_filter[, c(2:13)],   # columns 2–13 are the 12 loci
  ind.names = msdf_filter$Sample,
  pop       = msdf_filter$country,
  NA.char   = "0",
  sep       = ",",
  type      = "codom",
  ploidy    = 1                         # haploid organism
)
```

- `ploidy = 1` is critical: *P. falciparum* is a haploid parasite.
- `NA.char = "0"` maps the integer zero to a missing genotype.
- `pop` assigns each sample to its geographic population for all
  group-level analyses.

---

### 5.3 Population Diversity Statistics

**Goal:** Compute a suite of diversity indices for each geographic
population and for the combined dataset. Produces **Table 1** in the
manuscript (`figures/table_microsattelilte-diverity.png`).

```r
library(gt)

df_poppr <- poppr(pfms, total = TRUE, missing = "ignore",
                  clonecorrect = FALSE)

# Formatted gt table (all columns except administrative ones)
gt_table <- df_poppr |>
  select(-c(File, MLG, eMLG, SE)) |>
  gt() |>
  tab_header(title = "Microsatellite Diversity by Population") |>
  cols_label(
    Pop    = "Population",
    N      = "Sample Size (N)",
    H      = "Shannon-Weiner Diversity Index (H)",
    G      = "Stoddard-Taylor's Index (G)",
    lambda = "Simpson's Index (λ)",
    E.5    = "Evenness (E.5)",
    Hexp   = "Expected Heterozygosity (Hexp)",
    Ia     = "Index of Association (Ia)",
    rbarD  = "Standardized Index of Association (rbarD)"
  ) |>
  fmt_number(columns = vars(N, H, G, lambda, E.5, Hexp, Ia, rbarD),
             decimals = 2)

gtsave(gt_table, "population_statistics.html")
```

**Diversity indices at a glance:**

| Index | Interpretation |
|---|---|
| H (Shannon-Weiner) | Overall genotypic diversity; higher = more diverse |
| G (Stoddard-Taylor) | Effective number of genotypes |
| λ (Simpson's) | Probability that two randomly chosen samples differ |
| E.5 (Evenness) | How evenly genotypes are distributed |
| Hexp (Expected heterozygosity) | Mean allelic diversity across loci |
| Ia / rbarD (Index of Association) | Linkage disequilibrium signal; near 0 = random mating |

---

### 5.4 Locus-Specific Descriptors

**Goal:** Count the number of unique alleles at each microsatellite
locus in each population. Produces the allele-count bar chart and
**Table 2** (`figures/table_loci-unique-alleles.png`).

```r
# Compute per-population locus tables
locus_table_kenya    <- locus_table(pfms, population = "Kenya",
                                    index = "simpson", lev = "allele")
locus_table_ghana    <- locus_table(pfms, population = "Ghana",
                                    index = "simpson", lev = "allele")
locus_table_thailand <- locus_table(pfms, population = "Thailand",
                                    index = "simpson", lev = "allele")
locus_table_colombia <- locus_table(pfms, population = "Colombia",
                                    index = "simpson", lev = "allele")

# Assemble into a wide data frame
alleles_combined <- data.frame(
  Locus    = rownames(locus_table_kenya),
  Kenya    = as.integer(locus_table_kenya[, "allele"]),
  Ghana    = as.integer(locus_table_ghana[, "allele"]),
  Thailand = as.integer(locus_table_thailand[, "allele"]),
  Colombia = as.integer(locus_table_colombia[, "allele"])
)

# Bar plot
alleles_long <- alleles_combined |>
  pivot_longer(cols = -Locus, names_to = "Population",
               values_to = "Allele_Count")

ggplot(alleles_long, aes(x = Locus, y = Allele_Count, fill = Population)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Kenya"    = "#FF5800",
                                "Ghana"    = "#0072B5FF",
                                "Thailand" = "#EE4C97FF",
                                "Colombia" = "#20854EFF")) +
  labs(title = "Allele Counts by Locus and Population",
       x = "Microsatellite Locus", y = "Unique Alleles") +
  theme_minimal()
```

---

### 5.5 Minimum Spanning Network

**Goal:** Visualise genetic relationships among samples using a minimum
spanning network (MSN) built on pairwise dissimilarity distances.
Produces `figures/min-spanning-network.png`.

```r
library(igraph)

# 1. Pairwise dissimilarity (proportion of differing alleles)
adist <- diss.dist(pfms, percent = TRUE)

# 2. Build the MSN
amsn <- poppr.msn(pfms, adist, showplot = FALSE)

# 3. Plot with reproducible layout
set.seed(500)
msn.p <- plot_poppr_msn(
  pfms, amsn,
  gadj     = 5,
  layfun   = layout_as_tree,
  mgl      = TRUE,
  vertex.color = c("#FF5800",   # Kenya
                   "#A1378B",   # Ghana
                   "#003A59",   # Thailand
                   "#416322FF") # Colombia
)
```

The node size is proportional to the number of samples with that
multilocus genotype (MLG); edge width reflects genetic distance.

---

### 5.6 Genotype Accumulation Curve

**Goal:** Determine the minimum number of loci needed to discriminate
all multilocus genotypes. This serves as a quality-control check that
12 loci provide sufficient resolution.

```r
genotype_curve(
  pfms,
  sample  = 100,   # 100 random permutations per locus count
  quiet   = FALSE,
  plot    = TRUE,
  drop    = FALSE, # keep monomorphic loci
  dropna  = TRUE
)
```

The plateau of the curve indicates that ~10 loci are sufficient to
resolve all genotypes in this dataset.

---

### 5.7 Discriminant Analysis of Principal Components (DAPC)

**Goal:** Detect genetic clusters without assuming Hardy-Weinberg
equilibrium and project samples onto discriminant axes. The DAPC is
implemented in `adegenet`.

#### Step 1 — Find the optimal number of clusters (interactive)

> This step **requires interactive input** and cannot run in
> non-interactive batch mode. Run it in an R console or RStudio.

```r
# Evaluate up to 10 clusters; follow the on-screen prompts to choose
# the number of retained PCs and clusters (recommended: 5 PCs, 5 clusters
# for this dataset)
clust <- find.clusters(pfms, max.n.clust = 10)
```

When prompted:
1. **"Choose the number of PCs to retain"** — type `5` and press Enter.
2. **"Choose the number of clusters"** — inspect the BIC plot; the
   elbow is typically at **5 clusters** for this dataset. Type `5` and
   press Enter.

#### Step 2 — Run DAPC

```r
# Again, follow the prompts to choose how many PCs and discriminant
# functions to retain (recommended: 5 PCs, 4 discriminant axes)
dapc1 <- dapc(pfms, clust$clust)
```

#### Step 3 — Quick scatter plot (base graphics)

```r
myCol <- c("red", "#827100", "forestgreen", "blue", "purple", "green")

scatter(dapc1,
        posi.da = "none", posi.pca = "topleft",
        bg = "white", cstar = 0, col = myCol,
        scree.pca = TRUE, cex = 2, solid = 0.7)
```

---

### 5.8 Final PCA Plot

**Goal:** Produce a publication-quality scatter plot of the DAPC
coordinates with confidence ellipses per population, marginal
histograms, and labelled reference strains. Saves
`figures/PCA.png`.

```r
library(ggplot2); library(ggrepel); library(ggside); library(ggforce)

# Extract coordinates
coords   <- dapc1$ind.coord
dapc_data <- data.frame(
  SampleID = rownames(coords),
  PC1      = coords[, 1],
  PC2      = coords[, 2],
  Region   = dapc1$grp
)

# Label only the five reference strains
reference_ids <- c("3D7", "7G8", "Dd2", "FVO(NIH)", "NF54")
dapc_data <- dapc_data |>
  mutate(Label = ifelse(SampleID %in% reference_ids, SampleID, ""))

dapc_data_no_ref <- dapc_data |> filter(!SampleID %in% reference_ids)

# Variance explained
pct_var   <- dapc1$pca.eig / sum(dapc1$pca.eig) * 100
pc1_var   <- round(pct_var[1], 2)
pc2_var   <- round(pct_var[2], 2)

region_colors <- c(
  "Thailand" = "#EE4C97FF", "Cambodia" = "#EE4C97FF",
  "Colombia" = "#20854EFF", "Ghana"    = "#0072B5FF",
  "Kenya"    = "#FF5800"
)

scatter_plot <- ggplot(dapc_data, aes(x = PC1, y = PC2, fill = Region)) +
  geom_mark_ellipse(data = dapc_data_no_ref,
                    aes(fill = Region), alpha = 0.1,
                    color = NA, expand = unit(2, "mm"),
                    show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.7, shape = 21, color = "black") +
  geom_label_repel(aes(label = Label), size = 3,
                   fontface = "bold", fill = "white",
                   box.padding = 0.35, point.padding = 0.5,
                   segment.color = "grey50", max.overlaps = 30) +
  scale_color_manual(values = region_colors) +
  scale_fill_manual(values = region_colors) +
  labs(x = paste("PC1 -", pc1_var, "% variance explained"),
       y = paste("PC2 -", pc2_var, "% variance explained"),
       color = "Region") +
  theme_bw() +
  theme(panel.grid  = element_blank(),
        panel.border = element_blank(),
        axis.line    = element_line(color = "black"))

# Add marginal histograms
scatter_plot_with_margins <- scatter_plot +
  geom_xsidehistogram(position = "stack") +
  geom_ysidehistogram(position = "stack") +
  theme(legend.position    = c(0.65, 0.08),
        legend.justification = c("left", "bottom"))

ggsave("figures/PCA.png", scatter_plot_with_margins,
       height = 5, width = 5, dpi = 500)
```

**What to look for:**
- Samples from the same geographic region should cluster together.
- Reference strains (3D7, Dd2, etc.) are labelled and typically
  scatter among their closest genetic relatives.
- The marginal histograms summarise the PC1/PC2 distribution per
  cluster.

---

### 5.9 Population Statistics with Confidence Intervals

**Goal:** Estimate bootstrap confidence intervals around the diversity
indices from §5.3 and visualise them as point-and-error-bar plots.
Saves `figures/poppr_ci.svg`.

```r
poppr_diversity_ci <- diversity_ci(
  pfms,
  n      = 100L,  # bootstrap iterations
  n.boot = 1L,
  ci     = 95,
  total  = FALSE,
  plot   = TRUE,
  center = TRUE
)

ggsave("figures/poppr_ci.svg", last_plot(), height = 4, width = 6)
```

The subsequent loop in the notebook plots each diversity metric (H, G,
λ, E.5) separately with error bars, using the merged observed and
estimated data frames.

---

### 5.10 Genetic Distance Matrices

**Goal:** Compute pairwise genetic distances between all samples using
Nei's and Rogers' distance metrics and optionally export the results.

```r
library(writexl)

# Nei's distance (standard genetic distance)
nei         <- nei.dist(pfms, warning = TRUE)
nei_df      <- as.data.frame(as.matrix(nei))

# Rogers' distance (chord distance; better for small samples)
rogers      <- rogers.dist(pfms)
rogers_df   <- as.data.frame(as.matrix(rogers))

# Save to CSV / Excel (uncomment as needed)
# write_csv(nei_df,    "results/microsattelite_nei.csv")
# write_xlsx(nei_df,   "results/microsattelite_nei.xlsx")
write_csv(rogers_df,  "results/microsattelite_rogers.csv")
write_xlsx(rogers_df, "results/microsattelite_rogers.xlsx")

# Additional distance metrics available in poppr/adegenet:
edwards.dist(pfms)
reynolds.dist(pfms)
provesti.dist(pfms)
```

> **Note:** A `results/` directory must exist before saving. Create it
> with `dir.create("results", showWarnings = FALSE)` if needed.

---

### 5.11 Package Citations

The final chunk generates formatted citations for all key packages:

```r
packages <- c("poppr", "adegenet", "ggplot2", "seqvisr")

for (pkg in packages) {
  print(citation(pkg))
}
```

Please cite these packages in any publication using this analysis
pipeline.

---

## 6. Output Files & Figures

| File | Section | Description |
|---|---|---|
| `microsattelite.html` | All | Rendered HTML report with all figures inline |
| `population_statistics.html` | §5.3 | `gt` table of population diversity indices |
| `figures/PCA.png` | §5.8 | DAPC scatter plot with marginal histograms |
| `figures/PCA_narrow.png` | §5.8 | Narrower version of the PCA plot |
| `figures/min-spanning-network.png` | §5.5 | Minimum spanning network |
| `figures/population_barplot.png` | §5.4 | Allele-count bar plot by population |
| `figures/table_loci-unique-alleles.png` | §5.4 | Table of unique alleles per locus |
| `figures/table_microsattelilte-diverity.png` | §5.3 | Diversity table screenshot |
| `figures/poppr_ci.svg` | §5.9 | Diversity CI plot |
| `results/microsattelite_rogers.csv` | §5.10 | Rogers pairwise distance matrix (CSV) |
| `results/microsattelite_rogers.xlsx` | §5.10 | Rogers pairwise distance matrix (Excel) |

---

## 7. Troubleshooting

**`find.clusters()` hangs or produces no output**
: This function opens an interactive graphical prompt. It requires an
  R session with graphics support (RStudio, X11, or macOS Quartz).
  It cannot be run in a headless terminal. Run the two `find.clusters`
  and `dapc` chunks manually in RStudio, then render the rest of the
  document.

**Font "helvetica" not found**
: The theme sets `base_family = "helvetica"`. On Linux this may not be
  installed. Replace with `"sans"` or install the font via your system
  package manager (`apt install fonts-tex-gyre` on Ubuntu).

**`gtsave()` produces an empty HTML file**
: Ensure the `webshot2` package is installed:
  ```r
  install.packages("webshot2")
  webshot2::install_chromium()
  ```

**`write_xlsx()` fails with "no such file or directory"**
: The `results/` directory does not exist. Create it first:
  ```r
  dir.create("results", showWarnings = FALSE)
  ```

**Package compilation errors on macOS**
: Install Xcode Command Line Tools and `gfortran`:
  ```bash
  xcode-select --install
  # Then download gfortran from https://github.com/fxcoudert/gfortran-for-macOS/releases
  ```
