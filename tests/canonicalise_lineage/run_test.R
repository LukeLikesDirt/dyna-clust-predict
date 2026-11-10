#!/usr/bin/env Rscript
# run_test.R — Correctness test for canonicalise_lineage() (R/utils.R) and the
# R/canonicalise_lineage.R step that 01_reformat_ITS.sh runs after reformat.R.
#
# Fixture design -- one scenario per group of records:
#   r1-r4   genus Gx under family FamA x3 and FamB x1 -> all FamA (modal parent)
#   r5-r6   FamB's other genus Gy keeps FamB, and FamB keeps its own order
#   r7-r9   genus Gi under an incertae family x2 and FamC x1 -> FamC (an
#           identified parent beats an incertae one, even when outnumbered)
#   r10-r11 genus Gt under FamD x1 and FamE x1 -> FamD (tie: alphabetical)
#   r12-r13 genus Go with its family missing in one record -> filled from the
#           other record (and that record's order and class follow)
#   r14     genus Gn whose family is missing everywhere -> left alone (orphan)
#   r15-r16 phylum Glomeromycota under kingdom _pseudogene in r15: the marker
#           record neither votes nor changes; r16 keeps Fungi
#   r17     family FamA's order also differs in this one record (OrdZ vs OrdA)
#           -> OrdA, so a change at family carries up through order
#
# Usage: Rscript tests/canonicalise_lineage/run_test.R
# Note:  Must be run from the project root directory.

suppressPackageStartupMessages(library(data.table))

source("R/utils.R")

n_pass <- 0L
n_fail <- 0L
check <- function(label, ok) {
  if (isTRUE(ok)) { cat("PASS  ", label, "\n"); n_pass <<- n_pass + 1L }
  else            { cat("FAIL  ", label, "\n"); n_fail <<- n_fail + 1L }
}

# identical() on data.tables also compares data.table's internal attributes
same <- function(a, b) identical(as.data.frame(a), as.data.frame(b))

rec <- function(id, k, p, c, o, f, g, s, extra = "x") {
  data.table(id = id, kingdom = k, phylum = p, class = c, order = o, family = f,
             genus = g, species = s, its_complete = extra)
}

fx <- rbindlist(list(
  rec("r1",  "Fungi", "Asco", "Sord", "OrdA", "FamA", "Gx", "Gx a"),
  rec("r2",  "Fungi", "Asco", "Sord", "OrdA", "FamA", "Gx", "Gx b"),
  rec("r3",  "Fungi", "Asco", "Sord", "OrdA", "FamA", "Gx", "Gx c"),
  rec("r4",  "Fungi", "Asco", "Leca", "OrdB", "FamB", "Gx", "Gx d"),
  rec("r5",  "Fungi", "Asco", "Leca", "OrdB", "FamB", "Gy", "Gy a"),
  rec("r6",  "Fungi", "Asco", "Leca", "OrdB", "FamB", "Gy", "Gy b"),
  rec("r7",  "Fungi", "Asco", "Sord", "OrdA", "OrdA_fam_Incertae_sedis", "Gi", "Gi a"),
  rec("r8",  "Fungi", "Asco", "Sord", "OrdA", "OrdA_fam_Incertae_sedis", "Gi", "Gi b"),
  rec("r9",  "Fungi", "Asco", "Sord", "OrdC", "FamC", "Gi", "Gi c"),
  rec("r10", "Fungi", "Asco", "Sord", "OrdE", "FamE", "Gt", "Gt a"),
  rec("r11", "Fungi", "Asco", "Sord", "OrdD", "FamD", "Gt", "Gt b"),
  rec("r12", "Fungi", "Asco", "Sord", "OrdA", "FamA", "Go", "Go a"),
  rec("r13", "Fungi", "Asco", "unidentified", "unidentified", "unidentified", "Go", "Go b"),
  rec("r14", "Fungi", "Asco", "Sord", "OrdA", "unidentified", "Gn", "Gn a"),
  rec("r15", "_pseudogene", "Glomeromycota", "Glomeromycetes", "Entrophosporales",
      "Entrophosporaceae", "Entrophospora", "Entrophospora sp01"),
  rec("r16", "Fungi", "Glomeromycota", "Glomeromycetes", "Entrophosporales",
      "Entrophosporaceae", "Entrophospora", "Entrophospora infrequens"),
  rec("r17", "Fungi", "Asco", "Sord", "OrdZ", "FamA", "Gz", "Gz a")
))

res <- canonicalise_lineage(fx)
out <- res$classification
row <- function(i) out[id == i]

cat("-- modal parent, and the change carries up the lineage --\n")
check("r4: Gx re-parented to FamA (3 of 4 records)", row("r4")$family == "FamA")
check("r4: ... and inherits FamA's order and class", row("r4")$order == "OrdA" && row("r4")$class == "Sord")
check("r5/r6: Gy keeps FamB and FamB keeps OrdB",
      all(out[id %in% c("r5", "r6"), family] == "FamB") && all(out[id %in% c("r5", "r6"), order] == "OrdB"))
check("r17: FamA's minority order OrdZ -> OrdA", row("r17")$order == "OrdA")

cat("\n-- identified beats incertae; ties alphabetical --\n")
check("r7/r8: Gi's incertae family replaced by FamC (identified, 1 of 3)",
      all(out[id %in% c("r7", "r8"), family] == "FamC"))
check("r7/r8: ... and take FamC's order", all(out[id %in% c("r7", "r8"), order] == "OrdC"))
check("r10/r11: tie FamD/FamE goes to FamD", all(out[id %in% c("r10", "r11"), family] == "FamD"))
check("r10: ... and FamD's order", row("r10")$order == "OrdD")

cat("\n-- missing parents --\n")
check("r13: Go's missing family filled from r12 (FamA)", row("r13")$family == "FamA")
check("r13: ... and its missing order and class follow", row("r13")$order == "OrdA" && row("r13")$class == "Sord")
check("r14: Gn has no identified family anywhere: left unidentified", row("r14")$family == "unidentified")
check("n_orphans counts r14 at genus", res$n_orphans[["genus"]] == 1L)

cat("\n-- marker kingdom --\n")
check("r15: _pseudogene record untouched", identical(unlist(row("r15")), unlist(fx[id == "r15"])))
check("r16: Glomeromycota stays in Fungi", row("r16")$kingdom == "Fungi")

cat("\n-- shape and bookkeeping --\n")
check("same ids in the same order", identical(out$id, fx$id))
check("other columns kept", identical(out$its_complete, fx$its_complete))
check("species never change", identical(out$species, fx$species))
check("result is a tree (marker rows aside)", nrow(find_nesting_violations(out[!startsWith(kingdom, "_")])) == 0L)
check("input was not a tree", nrow(find_nesting_violations(fx[!startsWith(kingdom, "_")])) > 0L)
changed <- fx$id[vapply(seq_len(nrow(fx)), function(i) !identical(unlist(fx[i]), unlist(out[i])), logical(1))]
check("n_rows_changed counts the re-parented records",
      res$n_rows_changed == length(changed) && setequal(changed, c("r4", "r7", "r8", "r10", "r13", "r17")))
m <- res$manifest
check("manifest lists Gx with FamA chosen and 1 record changed",
      nrow(m[rank == "genus" & child == "Gx" & chosen_parent == "FamA" & n_records_changed == 1L]) == 1L)
check("manifest lists Gt with both parents", m[rank == "genus" & child == "Gt", parents] == "FamD(1); FamE(1)")
again <- canonicalise_lineage(out)
check("a tree is left unchanged", same(again$classification, out) && again$n_rows_changed == 0L && nrow(again$manifest) == 0L)

cat("\n-- refuses to move a record to another kingdom --\n")
xk <- rbindlist(list(rec("k1", "Metazoa", "PhyM", "C1", "O1", "F1", "G1", "G1 a"),
                     rec("k2", "Metazoa", "PhyM", "C1", "O1", "F1", "G1", "G1 b"),
                     rec("k3", "Viridiplantae", "PhyV", "C2", "O2", "F1", "G2", "G2 a")))
check("a family name shared by two kingdoms stops the run",
      inherits(tryCatch(canonicalise_lineage(xk), error = function(e) e), "error"))

cat("\n-- R/canonicalise_lineage.R end to end --\n")
dir.create("tests/canonicalise_lineage/tmp", showWarnings = FALSE, recursive = TRUE)
f_in  <- "tests/canonicalise_lineage/tmp/in.classification"
f_out <- "tests/canonicalise_lineage/tmp/out.classification"
f_man <- "tests/canonicalise_lineage/tmp/manifest.txt"
fwrite(fx, f_in, sep = "\t")
log <- system2("Rscript", c("R/canonicalise_lineage.R", "--classification_in", f_in,
                            "--classification_out", f_out, "--manifest_out", f_man),
               stdout = TRUE, stderr = TRUE)
check("script exits 0", is.null(attr(log, "status")))
got <- fread(f_out, sep = "\t", quote = "", colClasses = "character", na.strings = NULL)
check("script output equals canonicalise_lineage()", same(got, out))
check("script writes the manifest", nrow(fread(f_man, sep = "\t")) == nrow(m))
unlink("tests/canonicalise_lineage/tmp", recursive = TRUE)

cat(sprintf("\nRESULT: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0) quit(status = 1)
