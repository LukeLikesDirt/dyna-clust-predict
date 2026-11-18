# dyna-clust-predict

Prediction of optimal sequence similarity cut-offs for classification
and clustering of metabarcoding data using **vsearch global alignment**
and **F-measure optimisation** as a confidence metric, following Vu *et
al.* (2022).

## Overview

`dyna-clust-predict` predicts optimal sequence similarity thresholds
across taxonomic ranks using a curated ITS reference database
(EUKARYOME).

The pipeline:

- Extracts ITS1 and ITS2 subregions using ITSx
- Computes global pairwise similarity using vsearch
- Optimises similarity cut-offs using F-measure
- Supports parallel processing for efficient large-scale analyses

The workflow is conceptually adapted from
[dnabarcoder](https://github.com/vuthuyduong/dnabarcoder), but:

- Uses **vsearch** for global alignment
- Supports scalable parallel computation
- Integrates similarity prediction within an R-based workflow

## Pipeline steps

  `scripts/01_reformat_ITS.sh`         Download EUKARYOME, reformat headers,
                                       extract taxonomy, canonicalise the lineage
                                       (one parent per taxon)

  `scripts/02_derep_and_clean.sh`      Length-filter, standardise infraspecific
                                       annotations, and LCA-dereplicate

 `scripts/03_extract_subregions.sh`   Extract ITS1 and ITS2 using ITSx

 `scripts/04_remove_complexes.sh`     Restrict to complete-span sequences and

                                       remove species-level complexes

  `scripts/05_prepare_subsets.sh`      Generate balanced prediction subsets

  `scripts/06_compute_sim.sh`          *(Optional)* Pre-compute similarity
                                       matrices

  `scripts/07a_predict_cutoffs.sh`     Predict optimal similarity cut-offs
                                       (all regions sequentially in one job)

  `scripts/07b_launch_parallel.sh`     Submit one SLURM job per region
                                       (full_ITS, ITS1, ITS2) in parallel

  `scripts/07c_predict_cutoffs_parallel.sh`  Worker script for a single region
                                             (called by `07b_launch_parallel.sh`)

  `scripts/08_consolidate_cutoffs.sh`  Fill gaps and repair monotonicity in
                                       each region's nested cutoff table

All scripts must be run from the **project root directory**.

### R modules

  `R/utils.R`               Shared utility functions (e.g., `is_identified()`)

  `R/reformat.R`            FASTA header parsing and taxonomy extraction

  `R/canonicalise_lineage.R` Gives every named taxon exactly one parent
                            (`canonicalise_lineage()` in `R/utils.R`)

  `R/check_annotations.R`   Infraspecific annotation standardisation

  `R/append_completeness.R` Annotates each ID's `its_complete` status from
                            ITSx's raw (pre-dereplication) ITS1/ITS2 lengths

  `R/remove_complexes.R`    Hub-guarded species-complex detection and removal

  `R/subset.R`              Balanced taxon subset generation

  `R/compute_sim.R`         Pairwise similarity computation using vsearch

  `R/predict.R`             Cut-off prediction (parallel or sequential)

  `R/consolidate_cutoffs.R` Fill gaps and repair monotonicity via a
                            confidence-ranked fallback chain (self ->
                            ancestor taxa -> eukaryome-wide global)

## Directory structure

    data/
    ├── full_ITS/    # Full ITS sequences, taxonomy, ID files, predictions
    ├── ITS1/        # ITS1 sequences, taxonomy, ID files, predictions
    └── ITS2/        # ITS2 sequences, taxonomy, ID files, predictions

    output/
    ├── full_ITS/    # eukaryome_cutoffs.txt -- final consolidated cutoffs
    ├── ITS1/
    └── ITS2/

## Environment setup

A conda environment specification is provided in `environment.yml`.

[mamba](https://mamba.readthedocs.io/) is strongly recommended over
`conda` for faster dependency resolution.

### Local installation

```bash
mamba env create -f environment.yml
conda activate dyna_clust_predict
```

### HPC (SLURM)

Run from the project root:

```bash
sbatch scripts/create_dyna_clust_env.sh
conda activate dyna_clust_predict
```

## Quick start

From the project root:

```bash
sbatch scripts/01_reformat_ITS.sh
sbatch scripts/02_derep_and_clean.sh
sbatch scripts/03_extract_subregions.sh
sbatch scripts/04_remove_complexes.sh
sbatch scripts/05_prepare_subsets.sh
sbatch scripts/06_compute_sim.sh   # Optional

# Step 07 — choose one:
sbatch scripts/07a_predict_cutoffs.sh              # All regions in one job
bash   scripts/07b_launch_parallel.sh              # One job per region (parallel)

sbatch scripts/08_consolidate_cutoffs.sh           # Fill gaps, repair monotonicity
```

> **Note:** Step 06 is optional. Similarity can be computed on-the-fly
> in Step 07, which is preferred for large datasets and parallel
> execution.
>
> `07a` runs all three regions sequentially in a single SLURM job.
> `07b` submits three independent jobs (one per region, each running
> `07c_predict_cutoffs_parallel.sh`) so they run in parallel — faster
> overall but uses more nodes.

# Key parameters

Production values, set at the top of each script. The R scripts take the same
settings as flags (e.g. `MIN_OVERLAP` → `--min_overlap`); run
`Rscript R/<script>.R --help` for every option.

## 01 Reformat

No tunable parameters. Genus names shared across kingdoms get a `(Kingdom)`
suffix (`data/homonym_manifest.txt`). Taxa with more than one parent get their
most frequent parent (`data/lineage_canonicalisation.txt`).

## 02 Clean and dereplicate

| Parameter | Value | What it does |
|---|---|---|
| `MIN_LENGTH` | 250 | Minimum sequence length (bp) |
| `MAX_LENGTH` | 1500 | Maximum sequence length (bp) |
| `EXCLUDE_AMBIGUOUS` | TRUE | Drop sequences with ambiguous bases |

## 03 Extract subregions

| Parameter | Value | What it does |
|---|---|---|
| `MIN_SUBREGION_LENGTH` | 50 | Minimum ITS1 and ITS2 length (bp) for a sequence to count as complete |

## 04 Remove complexes

| Parameter | Value | What it does |
|---|---|---|
| `THRESHOLD` | 1.0 | Similarity at which two species form a complex |
| `IDDEF` | 2 | vsearch identity definition: terminal gaps excluded |
| `MIN_OVERLAP` | 250 / 50 / 50 | Minimum overlap (bp) for a pair to count: full ITS / ITS1 / ITS2 |
| `MAX_HUB_SPECIES` | 3 | Ignore sequences that match this many species or more |
| `MIN_SPECIES_PER_PARENT` | 2 | Skip genera with fewer species |
| `REQUIRE_COMPLETE` | yes | Keep only complete sequences |
| `MAX_SEQ_NO` | 20,000 | Downsample larger genera during detection |

Each complex is reduced to its best-sampled species
(`data/<region>/complex_manifest.txt`).

## 05 Prepare subsets

| Parameter | Value | What it does |
|---|---|---|
| `MIN_SUBGROUPS` | 10 | Minimum child taxa per dataset |
| `MIN_SEQUENCES` | 30 | Minimum sequences per dataset |
| `MAX_SEQUENCES` | 20,000 | Downsample larger datasets by whole child taxa, never by single sequences |
| `MAX_PROPORTION` | 0.7 | Maximum share of a dataset from one child taxon |
| `MAX_KINGDOM_PROPORTION` | 0.5 | Maximum share of a global dataset from one kingdom |
| `MAX_SEQS_PER_GROUP` | 10 | Maximum sequences per species (species datasets only) |
| `MAX_SINGLETON_PROPORTION` | 0.7 | Maximum share of singleton species (species datasets only) |

A singleton scores a free Dice = 1.0 at threshold 1.0, so singleton-heavy
datasets have their cutoff pinned near 1.0. Whole-taxon sampling and the
singleton cap limit this.

## 06 Compute similarity (optional)

| Parameter | Value | What it does |
|---|---|---|
| `MIN_SIM` | 0 | Lowest similarity written to the `.sim` file |
| `IDDEF`, `MIN_OVERLAP` | as 04 | |

## 07 Predict cutoffs

| Parameter | Value | What it does |
|---|---|---|
| `START_THRESH` | 0.8 species, 0.7 genus, 0.6 family, 0.5 higher | First threshold of the sweep |
| `END_THRESH` | 1 | Last threshold of the sweep |
| `STEP` | 0.001 | Sweep step |
| `MIN_MULTISEQ_GROUPS` | 10 | Minimum taxa with 2+ sequences to report a cutoff |
| `CUTOFF_ROUND_TO` | 0.005 | Round the reported cutoff to this step |
| `IDDEF`, `MIN_OVERLAP` | as 04 | |
| `REDO` | yes | Recompute instead of reusing cached results |

Writes to `data/<region>/`: `eukaryome.predicted` (F-measure at every
threshold) and `eukaryome_raw_cutoffs.txt` (best cutoff per dataset, with its
F-measure as `confidence`).

## 08 Consolidate cutoffs

No tunable parameters. Takes each cutoff from the dataset or one of its
ancestors, whichever has the highest confidence (the global value only as a
last resort), then raises cutoffs so they never decrease from coarse to fine
ranks. Writes
`output/<region>/eukaryome_cutoffs.txt`; `source` shows where each value came
from and `clamped` marks raised values.

## Similarity definition

Steps 04, 06 and 07 all use `IDDEF=2` and the same `MIN_OVERLAP`. Cutoffs are
only valid when sequences are compared with the same settings.

## Tests

`tests/fresh_pipeline_test/fresh_test.R` runs steps 05 and 07 on real data with
the values set in the scripts. It writes only to `--work`.

## Cross-repo harmonisation

`dyna-clust-predict-am` and `dyna-clust-predict-ecm` share this repo's core R
files. Files listed in `tools/sync_manifest.tsv` are sent to them as pull
requests on every push to `main` (needs the `SIBLING_REPOS_PAT` secret). Shell
scripts never sync. Check for drift with `tools/check_drift.sh`.

## Citation

Vu, D., Nilsson, R. H., & Verkley, G. J. (2022). Dnabarcoder: An open‐source software package for analysing and predicting DNA sequence similarity cutoffs for fungal sequence identification. Molecular Ecology Resources, 22(7), 2793-2809 https://doi.org/10.1111/1755-0998.13651
