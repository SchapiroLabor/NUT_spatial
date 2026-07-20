This repository consists of the preprocessing and downstream analysis scripts of the NUT Carcinoma MICS dataset, in collaboration between the Seitz and Schapiro groups at the Univeristy Hospital Tuebingen and University Hospital Heidelberg in Germany.


!["NUT Carcinoma tissue, inspired by Vincent van Gogh's The Starry Night, by Sophia Scheuermann"](starryNUT.jpeg)

## Overview

This is, to our knowledge, the first spatial immune profiling study of NUT carcinoma. Using a 72-marker MICS panel across a 13-patient FFPE cohort, we characterize the tumor immune microenvironment (TIME), with a focus on:

- Immune cell composition and the immunosupressive phenotype
- Spatial niche organization (including first-reported TLS-like structures in NC) 
- MDSC–CD8⁺ T cell spatial neighor preference 
- The association of the spatial and non-spatial features with clinical outcome

## Data

Relevant data is deposited to [TBD]: [TBD]

This includes:

- pre-processed images, as of the output of MACS iQ viewer
- the phenotyping labels, as of the output of the phenotyping process
- Processed single-cell table as .h5ad and .csv files

If you have any data-related questions, please contact the corresponding authors for access.


## Reproducing the analysis

In order to reproduce the analysis, please copy this repo to your local machine and follow the environment and configuration instructions. 

```bash
git clone https://github.com/SchapiroLabor/NUT_spatial
```

Make sure you have [git](https://git-scm.com/install/) installed on your machine

You can also find below the order of the scripts.
Paper figures are generated in **`figures/`**

## Environment

Analysis environments are managed with [pixi](https://pixi.sh), combining Python  and R. Each directory contains a separate pixi workspace member with its own environment for analysis or figure rendering

```bash
pixi install
```

## Configuration

All directory paths used across the pipeline (both Python and R scripts) are centralized in `config.yaml` at the repo root. Before running any part of the pipeline, **edit the paths in `config.yaml` to match your local machine or server** — no paths should be hardcoded elsewhere in the scripts.

```yaml
raw_dir: '/path/to/your/raw' #contain preprocessed macs iq data but not fully preprocessed
metadata_dir: '/path/to/your/metadata' #contains metadata
preprocessed_dir: '/path/to/your/preprocessed' # contains preprocessed data per ROI
results_dir: '/path/to/your/results' #results_dir contains anndata objects, generated tables, and figures folders
figure_1: "/path/to/your/results/figure_1" # contains figures 
```

Both the Python and R scripts read this same file, so you only need to update it in one place:

- **Python:** loaded via `yaml.safe_load(open("config.yaml"))`
- **R:** loaded via `yaml::read_yaml("config.yaml")`

Do not change folder names inside the scripts themselves — only the paths in `config.yaml` need to be adjusted to your setup.

## Pipeline

The analysis broadly follows this order:

1. **`preprocessing/`** — QC, marker normalization, construction of the unified `AnnData` object
2. **`proportions/`** — Cell type composition profiling (two-stage ROI → patient aggregation)
3. **`niches_clustering/`** — Spatial neighborhood clustering and biological neighborhood annotation 
4. **`spatial_distance/`** — Distance and tumor contour distance analyses
5. **`cozi/`** — Neighbor preference score using COZIpy (NEP)
6. **`survival_analysis/`** — Cox regression and Kaplan-Meier analyses linking spatial/compositional features to patient outcomes
7. **`figures/`** — Assembly of all manuscript and supplementary figures

## Manuscript

["TBD]