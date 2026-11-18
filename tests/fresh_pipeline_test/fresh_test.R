#!/usr/bin/env Rscript
# fresh_test.R — Fresh test of the pipeline's R stages (05 sampling, 07 cut-off
# prediction) on one region's real data, using the parameters that the
# production scripts themselves carry.
#
# Why a driver and not the scripts: 05/07a/07c need bash >= 4 (associative
# arrays), conda and SLURM, and macOS ships bash 3.2. So this reads the
# constants and the Rscript flag lists straight out of scripts/05_prepare_
# subsets.sh and scripts/07a_predict_cutoffs.sh and runs the same R code with
# them. Only paths, the worker count and the temp dir are overridden; any other
# flag that cannot be resolved from the script is an error, so the test cannot
# drift from the scripts. Everything is written under --work; data/ is only read.
#
# Stages (--stage, comma-separated):
#   05      R/subset.R with 05's flags, fixed seed, then pool-vs-sample chunk
#           statistics for the species-target files (where the singleton cap acts)
#   07      R/predict.R with 07a's flags on the species-target (and a few other)
#           datasets of at most --max_chunk sequences. Larger ones are left out
#           because predict.R holds the whole all-pairs matrix (about 25-35 GB at
#           20,000 sequences, more than a laptop has)
#   xcheck  the fast sweep below against predict.R's own F-measure trace on the
#           largest datasets from stage 07 (all thresholds must agree)
#   heavy   the largest datasets (--heavy), too big for predict.R here, with the
#           same vsearch settings and an equivalent sweep: vsearch --id (start
#           threshold - 0.001) keeps only the pairs the sweep can use (identical
#           F at every threshold, ~1000x smaller)
#
# Usage (project root):
#   Rscript tests/fresh_pipeline_test/fresh_test.R --stage 05,07,xcheck,heavy \
#       --work tmp/fresh_pipeline_test
# Needs vsearch and seqkit on PATH, and R packages optparse, data.table,
# jsonlite, igraph.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
  library(igraph)
})
source("R/utils.R")
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

opt <- parse_args(OptionParser(option_list = list(
  make_option("--stage", default = "05,07,xcheck,heavy"),
  make_option("--region", default = "full_ITS"),
  make_option("--work", default = "tmp/fresh_pipeline_test"),
  make_option("--results", default = "tests/fresh_pipeline_test/results"),
  make_option("--seed", type = "integer", default = 42L),
  make_option("--combos", default = "species:gen,species:fam,species:ord,species:cls,species:phy,species:kng,genus:fam,family:ord",
              help = "target:parent-abbreviation pairs for stage 07"),
  make_option("--max_chunk", type = "integer", default = 2000L),
  make_option("--xcheck_n", type = "integer", default = 4L),
  make_option("--heavy", default = "Fungi@kng,Viridiplantae@kng,Metazoa@kng,global",
              help = "species-target datasets for stage heavy: <taxon>@<parent abbr> or global"),
  make_option("--n_cpus", type = "integer", default = 4L, help = "predict.R workers (stage 07)"),
  make_option("--threads", type = "integer", default = 8L, help = "vsearch threads (stages xcheck, heavy)")
)))
stages  <- strsplit(opt$stage, ",")[[1]]
region  <- opt$region
work    <- opt$work
dir.create(work, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$results, recursive = TRUE, showWarnings = FALSE)
id_dir  <- file.path(work, region)
dir.create(id_dir, recursive = TRUE, showWarnings = FALSE)

# ── Read constants and flag lists out of the production scripts ───────────────

parse_bash <- function(path) {
  lines <- sub("\\s+#.*$", "", readLines(path))
  consts <- list(); arrays <- list()
  for (ln in lines) {
    m <- regmatches(ln, regexec('^readonly ([A-Z_0-9]+)=(?:"([^"]*)"|([^"[:space:]]*))\\s*$', ln, perl = TRUE))[[1]]
    if (length(m)) consts[[m[2]]] <- if (nzchar(m[3])) m[3] else m[4]
    m <- regmatches(ln, regexec('^([A-Z_0-9]+)\\["([^"]+)"\\]="?([^"]*)"?\\s*$', ln))[[1]]
    if (length(m)) arrays[[m[2]]][[m[3]]] <- m[4]
  }
  list(consts = consts, arrays = arrays)
}

# The flag list of the Rscript call whose header matches `header`, and (if given)
# whose text contains `must_contain` (07a has a local and a global call).
parse_block <- function(path, header, must_contain = NULL) {
  lines <- readLines(path)
  for (i in grep(header, lines)) {
    j <- i; blk <- character()
    repeat {
      blk <- c(blk, lines[j])
      if (!grepl("\\\\\\s*$", lines[j])) break
      j <- j + 1L
    }
    txt <- paste(sub("\\\\\\s*$", "", blk), collapse = " ")
    if (is.null(must_contain) || grepl(must_contain, txt, fixed = TRUE)) {
      m <- regmatches(txt, gregexpr('--[a-z_]+\\s+("[^"]*"|\\S+)', txt))[[1]]
      return(lapply(m, function(x) {
        flag <- sub("^--([a-z_]+)\\s+.*$", "\\1", x)
        list(flag = flag, value = sub('^--[a-z_]+\\s+', "", x))
      }))
    }
  }
  stop("no Rscript call matching ", header, " in ", path)
}

resolve_value <- function(v, p, ctx) {
  v <- gsub('^"|"$', "", v)
  if (!grepl("\\$", v)) return(v)
  m <- regmatches(v, regexec("^\\$\\{?([A-Z_0-9]+)\\}?$", v))[[1]]
  if (length(m) && !is.null(p$consts[[m[2]]])) return(p$consts[[m[2]]])
  m <- regmatches(v, regexec("^\\$\\{([A-Z_0-9]+)\\[\\$([a-z_]+)\\]\\}$", v))[[1]]
  if (length(m) && !is.null(p$arrays[[m[2]]][[ctx[[m[3]]]]])) return(p$arrays[[m[2]]][[ctx[[m[3]]]]])
  if (grepl("^\\$\\{?st\\}?$", v)) return(p$arrays$START_THRESH[[ctx$target]])
  m <- regmatches(v, regexec("^\\$\\{?([a-z_]+)\\}?$", v))[[1]]
  if (length(m) && !is.null(ctx[[m[2]]])) return(as.character(ctx[[m[2]]]))
  stop("cannot resolve '", v, "' from the script")
}

# argv for an Rscript call: every flag resolved from the script, except those in
# `override`. An overridden flag that the script does not have is an error too.
build_args <- function(block, p, ctx, override) {
  flags <- vapply(block, function(b) b$flag, character(1))
  stopifnot(all(names(override) %in% flags))
  unlist(lapply(block, function(b) {
    v <- if (b$flag %in% names(override)) override[[b$flag]] else resolve_value(b$value, p, ctx)
    c(paste0("--", b$flag), v)
  }))
}

run_r <- function(script, args, log, env = character()) {
  cat("  $ Rscript", script, paste(args, collapse = " "), "\n")
  st <- system2("Rscript", c(script, args), stdout = log, stderr = log, env = env)
  if (st != 0) stop("failed (", st, "): see ", log)
  invisible(st)
}

p05 <- parse_bash("scripts/05_prepare_subsets.sh")
p07 <- parse_bash("scripts/07a_predict_cutoffs.sh")
fasta_in <- p07$arrays$REGION_FASTA[[region]]
class_in <- p07$arrays$REGION_CLASS[[region]]
abbr_to_rank <- setNames(names(p07$arrays$RANK_ABBR), unlist(p07$arrays$RANK_ABBR))
cat("Region ", region, "\n  FASTA ", fasta_in, "\n  class ", class_in, "\n", sep = "")
cat("  05: ", paste(sprintf("%s=%s", names(p05$consts), unlist(p05$consts)), collapse = " "), "\n", sep = "")
cat("  07: IDDEF=", p07$consts$IDDEF, " MIN_OVERLAP=", p07$arrays$MIN_OVERLAP[[region]],
    " MIN_MULTISEQ_GROUPS=", p07$consts$MIN_MULTISEQ_GROUPS,
    " CUTOFF_ROUND_TO=", p07$consts$CUTOFF_ROUND_TO, "\n", sep = "")

cls_all <- NULL
get_cls <- function() {
  if (is.null(cls_all)) cls_all <<- fread(class_in, sep = "\t", colClasses = "character")
  cls_all
}

# ── Stage 05 ──────────────────────────────────────────────────────────────────

chunk_stats <- function() {
  cls <- get_cls()[is_identified(species)]
  out <- list()
  for (f in list.files(id_dir, pattern = "^species_pred_id_.*\\.txt$", full.names = TRUE)) {
    key <- sub("^species_pred_id_(.*)\\.txt$", "\\1", basename(f))
    sel <- cls[id %in% readLines(f)]
    if (key == "global") {
      sel[, parent := "all"]; pool <- copy(cls)[, parent := "all"]
    } else {
      pr <- abbr_to_rank[[key]]
      sel[, parent := get(pr)]; pool <- copy(cls)[, parent := get(pr)][is_identified(get(pr))]
    }
    st <- function(d) d[, {
      s <- .SD[, .N, by = species]
      .(n_seq = .N, n_spp = nrow(s), n_multiseq = sum(s$N >= 2), single = mean(s$N == 1), dom = max(s$N) / .N)
    }, by = parent]
    m <- merge(st(pool)[parent %in% sel$parent], st(sel), by = "parent", suffixes = c("_pool", "_sel"))
    m[, parent_rank := key]
    out[[key]] <- m
  }
  rbindlist(out)
}

if ("05" %in% stages) {
  cat("\n== Stage 05: subset.R ==\n")
  seed_file <- file.path(work, "seed.R")
  writeLines(sprintf("set.seed(%d)", opt$seed), seed_file)
  block <- parse_block("scripts/05_prepare_subsets.sh", 'Rscript "\\$SUBSET"')
  args  <- build_args(block, p05, list(), list(classification_in = class_in, output_dir = id_dir))
  run_r("R/subset.R", args, file.path(work, sprintf("subset_%s.log", region)),
        env = paste0("R_PROFILE_USER=", seed_file))
  cs <- chunk_stats()
  fwrite(cs, file.path(opt$results, sprintf("stage05_species_chunks_%s.tsv", region)), sep = "\t")
  cap <- as.numeric(p05$consts$MAX_SINGLETON_PROPORTION)
  cs[, over_cap := single_pool > cap + 1e-9]
  summ <- cs[, .(chunks = .N, chunks_pool_over_cap = sum(over_cap),
                 chunks_no_multiseq_species = sum(n_multiseq_sel == 0),
                 chunks_lt10_multiseq = sum(n_multiseq_sel < 10),
                 seq_pool = sum(n_seq_pool), seq_sample = sum(n_seq_sel)), by = parent_rank]
  summ <- summ[order(match(parent_rank, c("gen", "fam", "ord", "cls", "phy", "kng", "global")))]
  fwrite(summ, file.path(opt$results, sprintf("stage05_summary_%s.tsv", region)), sep = "\t")
  cat("\nSpecies-target chunks (singleton cap ", cap, "):\n", sep = ""); print(summ)
}

# ── Stage 07 ──────────────────────────────────────────────────────────────────

combo_dir <- function(combo) file.path(work, "pred", region, gsub(":", "_", combo))

if ("07" %in% stages) {
  cat("\n== Stage 07: predict.R (datasets of at most ", opt$max_chunk, " sequences) ==\n", sep = "")
  block <- parse_block("scripts/07a_predict_cutoffs.sh", 'Rscript "\\$PREDICT"', must_contain = "--higher_rank")
  cls <- get_cls()
  for (combo in strsplit(opt$combos, ",")[[1]]) {
    target <- sub(":.*", "", combo); pkey <- sub(".*:", "", combo)
    parent <- abbr_to_rank[[pkey]]
    id_file <- file.path(id_dir, sprintf("%s_pred_id_%s.txt", target, pkey))
    if (!file.exists(id_file)) { cat("  skip ", combo, ": no ID file\n", sep = ""); next }
    sel   <- cls[id %in% readLines(id_file)]
    sizes <- sel[, .N, by = c(parent)]
    keep  <- sizes[N <= opt$max_chunk][[parent]]
    ids   <- sel[get(parent) %in% keep]$id
    cat(sprintf("\n-- %s: %d of %d datasets (%d of %d sequences)\n", combo, length(keep),
                nrow(sizes), length(ids), nrow(sel)))
    if (length(ids) == 0) next
    d <- combo_dir(combo); unlink(d, recursive = TRUE); dir.create(file.path(d, "tmp"), recursive = TRUE)
    writeLines(ids, file.path(d, "ids.txt"))
    system2("seqkit", c("grep", "-f", file.path(d, "ids.txt"), fasta_in), stdout = file.path(d, "in.fasta"), stderr = FALSE)
    fwrite(cls[id %in% ids], file.path(d, "in.classification"), sep = "\t")
    args <- build_args(block, p07, list(region = region, target = target, parent = parent),
                       list(input = file.path(d, "in.fasta"), classification = file.path(d, "in.classification"),
                            n_cpus = as.character(opt$n_cpus), tmp_dir = file.path(d, "tmp"), out = d))
    run_r("R/predict.R", args, file.path(d, "predict.log"))
  }
  # collate
  res <- list(); allds <- list()
  for (combo in strsplit(opt$combos, ",")[[1]]) {
    d <- combo_dir(combo); f <- file.path(d, "eukaryome_raw_cutoffs.txt"); pj <- file.path(d, "eukaryome.predicted")
    if (!file.exists(pj)) next
    target <- sub(":.*", "", combo)
    pred <- fromJSON(pj, simplifyVector = FALSE)[[target]]
    allds[[combo]] <- data.table(combo = combo, dataset = names(pred),
      n_seq = vapply(pred, function(e) as.numeric(e[["sequence number"]] %||% NA), numeric(1)),
      n_groups = vapply(pred, function(e) as.numeric(e[["group number"]] %||% NA), numeric(1)),
      n_multiseq = vapply(pred, function(e) as.numeric(e[["multiseq group number"]] %||% NA), numeric(1)),
      cutoff_all = vapply(pred, function(e) as.numeric(e[["cut-off"]] %||% NA), numeric(1)))
    if (file.exists(f)) res[[combo]] <- cbind(combo = combo, fread(f, sep = "\t"))
  }
  ds  <- rbindlist(allds, fill = TRUE); reported <- rbindlist(res, fill = TRUE)
  fwrite(reported, file.path(opt$results, sprintf("stage07_reported_%s.tsv", region)), sep = "\t")
  gate <- as.integer(p07$consts$MIN_MULTISEQ_GROUPS)
  ds[, reported := paste(combo, dataset) %in% paste(reported$combo, reported$dataset)]
  summ <- ds[, .(datasets = .N, reported = sum(reported),
                 dropped_by_multiseq_gate = sum(!reported & n_multiseq < gate),
                 median_cutoff = median(cutoff_all[reported]),
                 q10 = quantile(cutoff_all[reported], 0.1), q90 = quantile(cutoff_all[reported], 0.9),
                 pct_at_or_above_0.995 = round(100 * mean(cutoff_all[reported] >= 0.995), 1)), by = combo]
  fwrite(summ, file.path(opt$results, sprintf("stage07_summary_%s.tsv", region)), sep = "\t")
  cat("\nStage 07 summary (reported = passed predict.R's own filters, incl. --min_multiseq_groups ", gate, "):\n", sep = "")
  print(summ)
}

# ── Fast sweep (stages xcheck and heavy) ──────────────────────────────────────
# F-measure at every threshold exactly as predict_dataset() computes it: single-
# linkage clusters at score >= t, per class the best Dice over clusters, weighted
# by class size, rounded to 4 d.p. before the tie-break (middle of the tied range,
# tie_tolerance 0), then the optimum rounded to --cutoff_round_to and re-synced.

fast_sweep <- function(edges, class_of, csz, n, start_t, end_t, step_t, round_to) {
  ts <- seq(round(start_t, 4), round(end_t, 4), by = step_t)
  setorder(edges, -score)
  tot <- sum(csz); fs <- numeric(length(ts))
  for (k in seq_along(ts)) {
    e <- edges[score >= ts[k] - 1e-9]
    g <- make_empty_graph(n, directed = FALSE)
    if (nrow(e)) g <- add_edges(g, as.vector(rbind(e$i, e$j)))
    memb <- components(g)$membership
    dt <- data.table(root = memb, cls = class_of)[, .N, by = .(root, cls)]
    rsz <- tabulate(memb, nbins = max(memb))
    dt[, d := 2 * N / (csz[cls] + rsz[root])]
    best <- dt[, .(b = max(d)), by = cls]
    fs[k] <- round(sum(csz[best$cls] * best$b) / tot, 4)
  }
  b <- max(fs); tied <- which(fs >= b)
  pos <- tied[ceiling(length(tied) / 2)]
  raw_t <- ts[pos]
  if (round_to > 0) {
    rt <- round(round(raw_t / round_to) * round_to, 4)
    pos <- which.min(abs(ts - rt))
  }
  list(ts = ts, fs = fs, cutoff_raw = raw_t, cutoff = ts[pos], f = fs[pos], f_at_end = fs[length(fs)])
}

# vsearch exactly as predict.R's compute_sim() makes it, plus --id <start>.
fast_dataset <- function(ids, sp_of, tag, target = "species") {
  d <- file.path(work, "fast"); dir.create(d, showWarnings = FALSE)
  start_t <- as.numeric(p07$arrays$START_THRESH[[target]]); end_t <- as.numeric(p07$consts$END_THRESH)
  step_t  <- as.numeric(p07$consts$STEP); round_to <- as.numeric(p07$consts$CUTOFF_ROUND_TO)
  iddef   <- as.integer(p07$consts$IDDEF); mo <- as.integer(p07$arrays$MIN_OVERLAP[[region]])
  writeLines(ids, f_ids <- file.path(d, paste0(tag, ".ids")))
  fa <- file.path(d, paste0(tag, ".fasta")); sim <- file.path(d, paste0(tag, ".sim"))
  system2("seqkit", c("grep", "-f", f_ids, fasta_in), stdout = fa, stderr = FALSE)
  # vsearch prints identity to one decimal but tests --id on the unrounded value,
  # so a pair printing as 80.0 can have identity 0.7996 and be dropped by --id 0.8.
  # Start 0.001 lower: any pair printing at or above the first threshold is kept,
  # and the F trace then equals predict.R's at every threshold (checked on Insecta).
  id_floor <- format(round(start_t - 0.001, 4))
  accept <- if (mo > 0) sprintf("--id %s --mincols %d", id_floor, mo) else sprintf("--id %s", id_floor)
  cmd <- sprintf("vsearch --allpairs_global '%s' %s --userout '%s' --userfields query+target+id --iddef %d --threads %d",
                 fa, accept, sim, iddef, opt$threads)
  cat("  $", cmd, "\n"); t0 <- Sys.time()
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  raw <- fread(sim, header = FALSE, col.names = c("q", "t", "pid"), colClasses = c("character", "character", "numeric"))
  idx <- setNames(seq_along(ids), ids); sp <- sp_of[ids]
  classes <- unique(sp); class_of <- match(sp, classes); csz <- tabulate(class_of, nbins = length(classes))
  edges <- data.table(i = idx[raw$q], j = idx[raw$t], score = round(raw$pid / 100, 4))
  r <- fast_sweep(edges, class_of, csz, length(ids), start_t, end_t, step_t, round_to)
  r$n <- length(ids); r$n_spp <- length(classes); r$n_multiseq <- sum(csz >= 2)
  r$pct_singleton <- round(100 * mean(csz == 1), 1); r$max_prop <- round(max(csz) / length(ids), 4)
  r$n_edges <- nrow(raw); r$elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")))
  unlink(c(fa, sim, f_ids)); r
}

# ── Stage xcheck ──────────────────────────────────────────────────────────────

if ("xcheck" %in% stages) {
  cat("\n== Stage xcheck: fast sweep vs predict.R's own trace ==\n")
  cls <- get_cls(); rows <- list()
  cand <- list()
  for (combo in strsplit(opt$combos, ",")[[1]]) {
    if (sub(":.*", "", combo) != "species") next
    pj <- file.path(combo_dir(combo), "eukaryome.predicted")
    if (!file.exists(pj)) next
    pred <- fromJSON(pj, simplifyVector = FALSE)$species
    for (nm in names(pred)) cand[[length(cand) + 1]] <- data.table(combo = combo, dataset = nm,
      n = as.numeric(pred[[nm]][["sequence number"]]), reported = !is.null(pred[[nm]][["fmeasures"]]))
  }
  cand <- rbindlist(cand)[order(-n)][!duplicated(dataset)][seq_len(min(opt$xcheck_n, .N))]
  for (i in seq_len(nrow(cand))) {
    combo <- cand$combo[i]; nm <- cand$dataset[i]
    parent <- abbr_to_rank[[sub(".*:", "", combo)]]
    ids <- readLines(file.path(combo_dir(combo), "ids.txt"))
    ids <- cls[id %in% ids & get(parent) == nm]$id
    sp_of <- setNames(cls$species, cls$id)
    f <- fast_dataset(ids, sp_of, paste0("x", i))
    pred <- fromJSON(file.path(combo_dir(combo), "eukaryome.predicted"), simplifyVector = FALSE)$species[[nm]]
    real <- vapply(sprintf("%.4f", f$ts), function(k) as.numeric(pred$fmeasures[[k]]), numeric(1))
    rows[[i]] <- data.table(combo = combo, dataset = nm, n = f$n, thresholds = length(f$ts),
      max_abs_diff = max(abs(real - f$fs)), cutoff_predict = as.numeric(pred[["cut-off"]]), cutoff_fast = f$cutoff,
      f_predict = as.numeric(pred[["confidence"]]), f_fast = f$f)
    cat(sprintf("  %-12s %-22s n=%5d  max |dF| over %d thresholds = %g  cutoff %.3f vs %.3f\n", combo, nm, f$n,
                length(f$ts), rows[[i]]$max_abs_diff, rows[[i]]$cutoff_predict, rows[[i]]$cutoff_fast))
  }
  xc <- rbindlist(rows)
  fwrite(xc, file.path(opt$results, sprintf("xcheck_%s.tsv", region)), sep = "\t")
}

# ── Stage heavy ───────────────────────────────────────────────────────────────

if ("heavy" %in% stages) {
  cat("\n== Stage heavy: datasets too large for predict.R on this machine ==\n")
  cls <- get_cls()[is_identified(species)]; sp_of <- setNames(cls$species, cls$id)
  out_f <- file.path(opt$results, sprintf("heavy_%s.tsv", region)); rows <- list()
  for (spec in strsplit(opt$heavy, ",")[[1]]) {
    if (spec == "global") { f <- file.path(id_dir, "species_pred_id_global.txt"); ids <- readLines(f) }
    else {
      nm <- sub("@.*", "", spec); pkey <- sub(".*@", "", spec); parent <- abbr_to_rank[[pkey]]
      ids <- readLines(file.path(id_dir, sprintf("species_pred_id_%s.txt", pkey)))
      ids <- cls[id %in% ids & get(parent) == nm]$id
    }
    cat("\n-- ", spec, ": ", length(ids), " sequences\n", sep = "")
    f <- fast_dataset(ids, sp_of, gsub("[^A-Za-z0-9]", "_", spec))
    rows[[spec]] <- data.table(dataset = spec, n = f$n, n_spp = f$n_spp, n_multiseq = f$n_multiseq,
      pct_singleton = f$pct_singleton, max_prop = f$max_prop, edges_ge_start = f$n_edges,
      cutoff_raw = f$cutoff_raw, cutoff = f$cutoff, f_at_cutoff = f$f, f_at_1 = f$f_at_end, elapsed_sec = f$elapsed)
    print(rows[[spec]])
    fwrite(rbindlist(rows), out_f, sep = "\t")
  }
}
cat("\nDone.\n")
