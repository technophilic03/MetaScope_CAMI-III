#!/usr/bin/env Rscript

# Build a cellular-only CAMI III gold standard.
#
# The 2025 target_reference holds no viral sequences, and filter_reference
# removes plasmid reads, so MetaScope can never report the Viruses (10239) or
# "other entries" (2787854) subtrees. In the human-gut-toy gold standard those
# two carry 36% to 51% of the profile. Leaving them in would make every
# cellular taxon read about 2x too abundant once metascope_to_cami.R rescales
# its own profile to 100.
#
# Every kept row is divided by that sample's own "cellular root" percentage.
# One factor per sample keeps the structure below it intact. Rescaling each
# rank to 100 instead would erase a real signal: order sums to 60, not 100,
# because some reads stop above that rank.
#
# Usage:
#   Rscript code/gold_standard_cellular.R --out FILE profile_0.txt [profile_1.txt ...]

suppressPackageStartupMessages({
  library(stringr)
})

cellular_root <- "131567"

args <- commandArgs(trailingOnly = TRUE)
out_at <- match("--out", args)
if (is.na(out_at) || out_at == length(args) || length(args) < 3) {
  stop(
    "Usage: gold_standard_cellular.R --out FILE profile_0.txt [profile_1.txt ...]"
  )
}
out_path <- args[out_at + 1]
inputs <- args[-c(out_at, out_at + 1)]

# The rows carry trailing empty columns that strsplit() drops, so every row is
# padded back to the width the "@@TAXID" line declares.
pad_to <- function(x, n) c(x, rep("", max(0, n - length(x))))

one_sample <- function(path) {
  lines <- readLines(path)
  lines <- lines[lines != ""]
  header <- lines[str_starts(lines, "@")]
  body <- lines[!str_starts(lines, "@")]

  width <- length(str_split(header[str_starts(header, "@@")][1], "\t")[[1]])
  fields <- lapply(str_split(body, "\t"), pad_to, n = width)

  taxpath <- vapply(fields, `[`, character(1), 3)
  # fixed(): "131567|" read as a regex is an alternation that matches everything.
  keep <- taxpath == cellular_root |
    str_starts(taxpath, fixed(str_c(cellular_root, "|")))

  root <- which(vapply(fields, `[`, character(1), 1) == cellular_root)
  if (length(root) != 1) {
    stop(path, ": expected exactly one row for taxid ", cellular_root)
  }
  root_pct <- as.numeric(fields[[root]][5])
  if (!is.finite(root_pct) || root_pct <= 0) {
    stop(path, ": cellular root percentage is ", fields[[root]][5])
  }
  scale <- 100 / root_pct

  kept <- fields[keep]
  kept <- lapply(kept, function(f) {
    f[5] <- sprintf("%.4f", as.numeric(f[5]) * scale)
    f
  })

  sample_id <- str_remove(header[str_starts(header, "@SampleID")][1], "^@SampleID:")
  species_sum <- sum(as.numeric(vapply(
    Filter(function(f) f[2] == "species", kept), `[`, character(1), 5
  )))
  message(sprintf(
    "  sample %-3s  rows %5d -> %5d   cellular root %7.4f -> 100   species sum %.4f",
    sample_id, length(fields), length(kept), root_pct, species_sum
  ))

  c(header, vapply(kept, str_c, character(1), collapse = "\t"))
}

message("Filtering ", length(inputs), " gold-standard profiles to the ",
        cellular_root, " subtree")
out <- unlist(lapply(inputs, one_sample), use.names = FALSE)
writeLines(out, out_path)
message("Wrote ", out_path)
