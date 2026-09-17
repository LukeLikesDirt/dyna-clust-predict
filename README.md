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

-   Extracts ITS1 and ITS2 subregions using ITSx
-   Computes global pairwise similarity using vsearch
-   Optimises similarity cut-offs using F-measure
-   Supports parallel processing for efficient large-scale analyses

The workflow is conceptually adapted from
[dnabarcoder](https://github.com/vuthuyduong/dnabarcoder), but:

-   Uses **vsearch** for global alignment
-   Supports scalable parallel computation
-   Integrates similarity prediction within an R-based workflow

## Pipeline steps

  `scripts/01_reformat_ITS.sh`         Download EUKARYOME, reformat headers,
                                       extract taxonomy

  `scripts/02_derep_and_clean.sh`      Length-filter, standardise infraspecific
                                       annotations, and LCA-dereplicate

  `scripts/03_extract_subregions.sh`   Extract ITS1 and ITS2 using ITSx;
                                       annotate completeness (`its_complete`)

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

``` bash
mamba env create -f environment.yml
conda activate dyna_clust_predict
```

### HPC (SLURM)

Run from the project root:

``` bash
sbatch scripts/create_dyna_clust_env.sh
conda activate dyna_clust_predict
```

## Quick start

From the project root:

``` bash
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

## Taxonomy reformatting

Used in: `reformat.R` (via `01_reformat_ITS.sh`)

Numbered Tedersoo et al. (2024, MycoKeys 107) placeholder codes (e.g.
`Densosporales.fam02`, `Glomeraceae.gen05`) are preserved as usable
taxonomic groups rather than discarded, with `.` converted to `_` so the
labels survive downstream taxonomy-string parsing.

Genus names that collide across more than one identified kingdom (e.g. a
plant genus and an unrelated animal genus sharing the same name) are
disambiguated by appending `(Kingdom)`, matching the convention EUKARYOME's
own curators already use for some known homonyms (e.g. `Achlya(Metazoa)` vs
`Achlya(Straminipila)`) — so pipeline-generated and upstream-supplied
disambiguation are indistinguishable downstream. Disambiguation happens
before species-name construction, so the species binomial inherits the
suffix automatically. Every resolved collision is written to a manifest so
a future EUKARYOME release that adds or drops a homonym shows up as a
reviewable diff rather than a silent rename.

    --manifest_out  FILE  Output path for the homonym manifest [default: data/homonym_manifest.txt]

Manifest columns: `original_genus`, `kingdom`, `n_sequences`,
`disambiguated_genus`.

## Completeness annotation

Used in: `append_completeness.R` (via `03_extract_subregions.sh`)

Flags each ID `its_complete = TRUE` iff its raw (pre-dereplication) ITS1 and
ITS2 extractions both exist and are each at least `--min_length` bp.
"Present in both ITS1/ITS2 outputs" alone is not sufficient — it only
confirms ITSx assigned *some* label, not that the label represents a
genuine, non-truncated region.

    --min_length  INT  Minimum length (bp) for each subregion [default: 50]

## Complex removal

Used in: `remove_complexes.R` (via `04_remove_complexes.sh`)

Detects species that cannot be discriminated by the barcode marker
(single-linkage clustering at `--threshold` identity, within each
`--parent_rank` taxon) and removes all but the best-sampled species per
complex. `--require_complete` (default `yes`) restricts to `its_complete`
rows first and drops incomplete rows from the output entirely — truncated
records inflate similarity broadly, not just for complex detection.

Not a parity port of dnabarcoder's `-removecomplexes`: that algorithm's
unguarded single-linkage merge is fragile to database contamination — one
mislabelled record can transitively chain unrelated, well-separated species
into a false mega-complex. `--max_hub_species` excludes any sequence whose
threshold-identity matches span at least that many distinct species from the
complex graph before clustering, guarding against exactly this.

    --rank                   STR    Rank at which complexes are detected [default: species]
    --parent_rank             STR    Rank within which complexes are detected [default: genus]
    --threshold               FLOAT  Identity threshold defining a complex edge [default: 1.0]
    --max_hub_species         INT    Exclude a sequence from complex detection if its matches
                                     span >= this many distinct species [default: 3]
    --min_species_per_parent  INT    Skip parent groups with fewer distinct species [default: 2]
    --require_complete        yes/no Restrict to its_complete rows; drop the rest [default: yes]
    --max_seq_no              INT    Max sequences per parent group; excess is randomly
                                     downsampled [default: 20000]
    --iddef                   0-4    vsearch pairwise identity definition [default: 2]

Example:

``` bash
Rscript R/remove_complexes.R \
  --fasta_in data/full_ITS/eukaryome_ITS.fasta \
  --classification_in data/full_ITS/eukaryome_ITS.classification \
  --fasta_out data/full_ITS/eukaryome_ITS_nocomplex.fasta \
  --classification_out data/full_ITS/eukaryome_ITS_nocomplex.classification \
  --manifest_out data/full_ITS/complex_manifest.txt \
  --threshold 1.0 \
  --max_hub_species 3 \
  --n_cpus 32
```

## Sequence selection

Used in: `subset.R` (via `05_prepare_subsets.sh`)

These parameters control taxonomic balance and sampling constraints
prior to similarity prediction.

    --min_subgroups   INT    Minimum unique child taxa per parent taxon (default: 10)
    --min_sequences   INT    Minimum sequences per parent taxon after proportion cap (default: 30)
    --max_sequences   INT    Maximum sequences per parent taxon; excess are balanced round-robin downsampled across child taxa (default: 25000)
    --max_proportion  FLOAT  Maximum fraction a child taxon may represent (default: 0.5)
    --max_kingdom_proportion  FLOAT  Maximum fraction of the STEP 2 global pool that the
                                     dominant kingdom may represent (default: 0.5). Independent
                                     of --max_proportion, which caps the target rank's own
                                     dominant clade -- inert for global pools since no single
                                     target-rank clade dominates, unlike kingdom composition.
    --max_seqs_per_group  INT   Max sequences any single species may contribute (default:
                                 0/off; production runs use 10). Applied only when
                                 target_rank == species -- a no-op at every coarser rank, where
                                 groups are naturally few and legitimately need deep sequence
                                 support. Species over the cap are thinned, never dropped. A
                                 species with many sequences and real internal diversity can't
                                 all be pairwise-identical near threshold 1.0, so it fragments
                                 there, creating useful downward pressure on the optimal
                                 cutoff -- but only if enough of its sequences survive to
                                 reveal that diversity, and a small number of species are
                                 oversampled in EUKARYOME for non-biological reasons (3.3% of
                                 Fungi species carry 26% of all sequences).
    --max_singleton_proportion  FLOAT  Max fraction of species (by count, not sequence count)
                                       that may be singletons (default: 1/off; production runs
                                       use 0.7). Applied only when target_rank == species, for
                                       the same reason as --max_seqs_per_group above. A
                                       singleton scores a free Dice = 1.0 at threshold 1.0, so
                                       datasets dominated by singleton species have their
                                       optimum pinned near 1.0 regardless of biology; excess
                                       singleton species are randomly dropped (never
                                       multi-sequence ones) to bring the fraction under the cap.

Both caps replace the previous row-level random downsample (Filter 4, and the
equivalent step in `global_prediction_filter`) with a group-aware downsample
that keeps a species' sequences together or drops it entirely -- a row-level
downsample silently fragments large species back down to one surviving
sequence once the sampling fraction is steep, manufacturing apparent
singletons that undo both caps above (confirmed on kingdom Fungi: raw pool is
50.3% singleton by species, but naive row-level downsampling from 118,171 to
5,000 sequences came out 89.0% singleton, with 76% of those "singleton"
species actually having >= 2 sequences in the real data).

Example:

``` bash
Rscript R/subset.R \
  --fasta_in input.fasta \
  --classification_in input.classification \
  --min_subgroups 10 \
  --min_sequences 30 \
  --max_sequences 25000 \
  --max_proportion 0.5 \
  --max_seqs_per_group 10 \
  --max_singleton_proportion 0.7 \
  --output_dir output
```

## Similarity prediction

Used in: `predict.R` (via `07a_predict_cutoffs.sh` / `07c_predict_cutoffs_parallel.sh`)

These parameters select the rank combination to predict:

    --rank          STR    Target rank(s), comma-separated (e.g. species,genus) [required]
    --higher_rank   STR    Parent rank(s) for local prediction, comma-separated (e.g. genus,family).
                           Omit to run a single global prediction across all sequences. [default: ""]

These parameters control the similarity sweep:

    --start_threshold   FLOAT   Starting similarity threshold (default: 0.0)
    --end_threshold     FLOAT   Ending similarity threshold (default: 1.0)
    --step              FLOAT   Threshold step size (default: 0.001)
    --min_cutoff        FLOAT   Min cutoff value to report in output (default: 0.0)

The same filtering thresholds used during sequence selection (`subset.R`) also apply here, but as
strict dataset filters rather than balancing rules — datasets that fail a threshold are
skipped entirely rather than downsampled:

    --min_group_no    INT    Min unique child taxa required to report a cutoff (default: 10)
    --min_seq_no      INT    Min sequences required to report a cutoff (default: 30)
    --max_seq_no      INT    Max sequences per dataset; excess is randomly sampled (default: 25000)
    --max_proportion  FLOAT  Skip datasets where the dominant group exceeds this fraction (default: 1.0)
    --min_multiseq_groups INT  Min groups with >= 2 sequences required to report a cutoff
                               (default: 0/off; production runs use 10). A group of size 1
                               scores a free Dice = 1.0 at threshold 1.0, so datasets
                               dominated by singleton groups have their optimum pinned at
                               1.0 regardless of biology -- min_group_no alone does not
                               catch this, since it counts groups of any size.

Threshold-selection controls:

    --tie_tolerance  FLOAT  Widens the tied-optimum-threshold selection from exact F-measure
                            equality to fmeasures >= best_f - tie_tolerance, still picking
                            the middle of the tied range (default: 0/off). Not used in
                            production -- exact-equality tie-breaking plus
                            --cutoff_round_to's coarser reporting grid was judged simpler,
                            without a meaningful accuracy tradeoff.
    --iddef          0-4    vsearch pairwise identity definition (default: 2, vsearch's own
                            default). Investigation found --iddef 1 costs 0.05-0.15
                            F-measure even on completeness-filtered data, so the default is
                            unchanged in production; exposed for future experimentation.
    --cutoff_round_to  FLOAT  Round the reported cut-off to the nearest multiple of this
                              value (default: 0/off; production runs use 0.005). The
                              confidence value is re-synced to the F-measure actually
                              observed at the rounded threshold, not the unrounded optimum,
                              so it always describes the reported cut-off. Combined with
                              subset.R's --max_seqs_per_group/--max_singleton_proportion,
                              chosen so kingdom Fungi's species-level cutoff lands at a
                              conventional, defensible ~0.97 rather than reporting spurious
                              precision like 0.9683.

Execution controls:

    --run_parallel   yes/no   Parallel dataset processing via furrr/future (default: yes)
    --n_cpus         INT      Workers (parallel) or vsearch threads (sequential) (default: all−1)
    --tmp_dir        DIR      Directory for temporary vsearch output (default: ./tmp)

Example:

``` bash
Rscript R/predict.R \
  --input data/full_ITS/eukaryome_ITS_nocomplex.fasta \
  --classification data/full_ITS/eukaryome_ITS_nocomplex.classification \
  --rank species \
  --higher_rank genus \
  --start_threshold 0.9 \
  --end_threshold 1.0 \
  --step 0.001 \
  --min_multiseq_groups 10 \
  --cutoff_round_to 0.005 \
  --run_parallel yes \
  --n_cpus 80 \
  --out data/full_ITS \
  --prefix eukaryome_ITS
```

Writes three files to `--out` (e.g. for `--prefix eukaryome`):

    eukaryome.predicted            Full F-measure traces per threshold (JSON)
    eukaryome_raw_cutoffs.json     Best cutoff per dataset, this run only (JSON)
    eukaryome_raw_cutoffs.txt      Same, tab-delimited

"Raw" because each row is one dataset's direct computation only, with gaps
wherever `subset.R`'s filters excluded a parent taxon. The pipeline's actual
deliverable is `consolidate_cutoffs.R`'s gap-filled, monotonicity-repaired
`eukaryome_cutoffs.txt` (next section) -- a different file, not this one.

## Cutoff consolidation

Used in: `consolidate_cutoffs.R` (via `08_consolidate_cutoffs.sh`)

`predict.R`'s nested cutoff table has one row per `(higher_rank, dataset,
rank)`, but a cell only exists where that dataset's subset passed
`subset.R`'s filters -- so some parent taxa are missing values at some
target ranks, and independently-computed ranks can occasionally violate the
constraint that similarity must increase from phylum to species (each rank
nests inside the one above it).

`consolidate_cutoffs.R` resolves every `(higher_rank, dataset, target rank)`
cell that has at least one direct computation somewhere in its lineage by
comparing all available candidates -- the dataset's own value and each
ancestor taxon's value at the same target rank (walking the real taxonomic
lineage derived from the classification file) -- and keeping whichever has
the highest confidence (F-measure), **excluding non-self candidates that
bring less multi-sequence evidence than self** (see `--min_multiseq_groups`:
a candidate with few groups of >= 2 sequences can still report an
artificially high, singleton-driven confidence, so it is not allowed to
override a self value backed by more real evidence merely on confidence).
The eukaryome-wide **global** value is a true last resort, used only when
self and every ancestor produced nothing at all -- it never competes on
confidence, because it pools every kingdom together and isn't a like-for-like
estimate for any one taxon (a marginally higher confidence there just
reflects a much larger, unrelated pool, not a better answer). It then clamps
each dataset's own resolved row to be non-decreasing from its coarsest to its
finest target rank; a row raised by that clamp has its `confidence` nulled
out, since nothing was actually measured at the raised threshold.

    --cutoffs_in         FILE   <prefix>_raw_cutoffs.txt for one region -- predict.R's raw,
                                pre-consolidation output [required]
    --classification_in  FILE   Region classification file, for lineage lookup only [required]
    --output             FILE   Output path for the consolidated table [required]

Requires the region's global (no `--higher_rank`) predictions to have
already been run via `07a`/`07c`, since the global cutoffs are the
gap-filler of last resort.

Output columns extend `predict.R`'s own (`rank`, `higher_rank`, `dataset`,
`cut-off`, `confidence`, `sequence number`, `group number`, `multiseq group
number`, `max proportion`) with:

    source                What supplied the winning value: self / <rank>:<name> / global
    clamped                TRUE if the monotonicity step raised this value (confidence is
                           then NA -- nothing was measured at the raised cutoff)
    original_cutoff        self's own direct value, if one existed (regardless of whether
                           self won -- NA if this taxon had no direct computation at all)
    original_confidence    Its confidence, if one existed

Example:

``` bash
Rscript R/consolidate_cutoffs.R \
  --cutoffs_in data/full_ITS/eukaryome_raw_cutoffs.txt \
  --classification_in data/full_ITS/eukaryome_ITS_nocomplex.classification \
  --output output/full_ITS/eukaryome_cutoffs.txt
```

## Cross-repo harmonisation

This is one of three sibling repos sharing a common R pipeline core:
`dyna-clust-predict` (this repo, the general eukaryote database),
`dyna-clust-predict-am` (AM/Glomeromycota-focused), and
`dyna-clust-predict-ecm` (ectomycorrhizal). They share `R/utils.R`,
`predict.R`, `reformat.R`, `compute_sim.R`, `dereplicate_lca.R`,
`consolidate_cutoffs.R`, and `remove_complexes.R` verbatim; AM has
intentionally diverged on `subset.R` and `check_annotations.R` for its own
sampling and correction logic. Note that only files listed in
`tools/sync_manifest.tsv` sync — `scripts/*.sh` (including the numbered
pipeline steps and where each new flag/parameter is wired in with its
production default) never syncs, so a sibling repo's own pipeline scripts
must be updated by hand to pick up any change made only at the shell-script
level here.

`tools/sync_manifest.tsv` declares which shared files sync to which
sibling. On every push to `main` touching a manifest-listed path,
[.github/workflows/sync-to-siblings.yml](.github/workflows/sync-to-siblings.yml)
opens a pull request in each applicable sibling with the updated file(s) —
never a direct push, since a change correct for the general pipeline can
still be wrong for a dataset-specific subset. Requires a `SIBLING_REPOS_PAT`
repository secret (a PAT with write access to both sibling repos); nothing
propagates until that secret exists.

To check for drift locally (e.g. from the HPC, without waiting on CI):

```bash
tools/check_drift.sh          # both siblings
tools/check_drift.sh am       # just dyna-clust-predict-am
```

Sync direction is one-way (main → siblings). A fix that originates in a
sibling is a manual PR into main, not something the workflow handles.

`dyna-clust-predict-ecm` has an unrelated git root commit (it was created by
copying files rather than forking), so `git merge`/`pull` from main will
never work there directly — the sync workflow applies file content
directly rather than merging, so this is a non-issue in practice.

## Citation

Vu, D., Nilsson, R. H., & Verkley, G. J. (2022). Dnabarcoder: An open‐source software package for analysing and predicting DNA sequence similarity cutoffs for fungal sequence identification. Molecular Ecology Resources, 22(7), 2793-2809 https://doi.org/10.1111/1755-0998.13651
