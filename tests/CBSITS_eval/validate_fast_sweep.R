# Validates fast_sweep() against sweep_f() from test_fungi_grid.R on one small pool.
suppressPackageStartupMessages({ library(data.table); library(igraph) })
source("R/utils.R")
src <- readLines("tests/CBSITS_eval/test_fungi_grid.R")
s <- grep("^sweep_f <- function", src); e <- grep("^run_vsearch_f <- function", src) - 1
start_t <- 0.90; end_t <- 1.0; step_t <- 0.001
eval(parse(text = src[s:e]))
args_save <- commandArgs; 
eval(parse(text = paste(readLines("tests/CBSITS_eval/test_random_vs_species_selection.R")[
  grep("^fast_sweep <- function", readLines("tests/CBSITS_eval/test_random_vs_species_selection.R")):
  (grep("^run_one <- function", readLines("tests/CBSITS_eval/test_random_vs_species_selection.R")) - 1)], collapse = "\n")))
fasta_all <- "data/full_ITS/eukaryome_ITS_nocomplex.fasta"
cls <- fread("data/full_ITS/eukaryome_ITS_nocomplex.classification", sep = "\t", colClasses = "character")
cls <- cls[is_identified(species) & kingdom == "Fungi"]
set.seed(1)
for (rep in 1:2) {
  ids <- cls[sample(.N, 1500)]$id; sp <- setNames(cls$species, cls$id)[ids]
  d <- "tmp/validate_sweep"; dir.create(d, showWarnings = FALSE); writeLines(ids, f <- file.path(d, "ids")); fa <- file.path(d, "x.fa")
  system(sprintf("seqkit grep -f '%s' '%s' > '%s' 2>/dev/null", f, fasta_all, fa))
  system(sprintf("vsearch --allpairs_global '%s' --acceptall --userout %s/all.sim --userfields query+target+id --iddef 2 --threads 10", fa, d), ignore.stdout = TRUE, ignore.stderr = TRUE)
  system(sprintf("vsearch --allpairs_global '%s' --id 0.90 --userout %s/id.sim --userfields query+target+id --iddef 2 --threads 10", fa, d), ignore.stdout = TRUE, ignore.stderr = TRUE)
  idx <- setNames(seq_along(ids), ids); classes <- unique(sp); class_of <- match(sp, classes); cs <- tabulate(class_of, length(classes))
  rd <- function(p) { r <- fread(p, header = FALSE, col.names = c("q","t","pid")); data.table(i = idx[r$q], j = idx[r$t], score = r$pid/100) }
  a <- sweep_f(rd(file.path(d, "all.sim")), class_of, cs, length(ids))
  b <- fast_sweep(rd(file.path(d, "id.sim")), class_of, cs, length(ids))
  cat(sprintf("rep %d: old cut=%.3f f=%.4f | new cut=%.3f f=%.4f | match=%s\n", rep, a$cut, a$f, b$cut, b$f, a$cut == b$cut && a$f == b$f))
}
