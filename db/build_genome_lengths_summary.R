# Build the taxid -> genome length table that the metascope_to_cami converters
# need, from NCBI's own per-assembly genome_size rather than by summing the
# reference index.
#
# WHY THIS EXISTS, AND WHAT WAS WRONG BEFORE
#
# build_genome_lengths.R walks the bowtie2 indices with bowtie2-inspect and does
#   group_by(taxid) |> summarise(length = sum(length))
# on the premise, stated in its own comment, that "contigs and plasmids of one
# genome share a taxid, so the sum over a taxid is that genome". That holds for
# one assembly per taxid. It does not hold for this reference, which carries
# 150,209 assemblies: E. coli has 2,472 of them and S. aureus 4,138, so the sum
# is thousands of genomes, not one.
#
# The resulting table gives E. coli a 1,457 Mb "genome" and K. pneumoniae
# 1,323 Mb. CAMI defines PERCENTAGE as cell abundance, and the converters get
# there by dividing read fraction by genome length, so a length inflated 317x
# makes that organism look 317x rarer. Measured on CAMI II strain-madness
# sample_0, the five truth species agreed with the gold standard at r = 0.81 as
# raw read fractions, 0.52 after dividing by the summed lengths, and 0.99 after
# dividing by the lengths this script writes. The old script's own comment warned
# about it: "Read the numbers before trusting the profile."
#
# WHY MEDIAN, AND WHY BOTH TAXID COLUMNS
#
# assembly_summary.txt carries one row per assembly with its own genome_size, so
# a per-taxid median is a representative genome rather than a pile of them. The
# median, not the mean, because a single mis-assembled or contaminated entry
# should not drag the estimate. Each assembly contributes to both its own taxid
# (column 6) and its species_taxid (column 7): a strain's genome size is that
# assembly, and a species' is the median over its assemblies. MetaScope reports
# whichever taxid accessionTaxa.sql resolved, so both need covering.
#
# WHY A FALLBACK
#
# 402 taxids in genome_lengths_new.tsv have no assembly_summary row -- mostly
# plasmid-only entries, median 0.04 Mb. For those the old summed length IS a
# single genome, and 401 of the 402 are under 20 Mb, so it is kept. The
# converters stop on a taxid with no length, so a gap here is a failed job
# rather than a silent error, and the fallback is what keeps coverage whole.
#
#   Rscript --vanilla db/build_genome_lengths_summary.R \
#     db/CAMI-III_reference/assembly_summary.txt \
#     db/genome_lengths_new.tsv \
#     db/genome_lengths_summary.tsv
#
# Reads two text files, so unlike build_genome_lengths.R it needs no
# bowtie2-inspect pass and no SQLite query, and runs in seconds on a login node.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: build_genome_lengths_summary.R <assembly_summary.txt> <fallback.tsv> <out.tsv>")
}
summary_path  <- args[1]
fallback_path <- args[2]
out_path      <- args[3]
for (f in c(summary_path, fallback_path)) {
  if (!file.exists(f)) stop("not found: ", f)
}

# The CAMI copy carries no header and no comment lines, but guard anyway.
raw <- read_tsv(
  summary_path,
  col_names = FALSE,
  comment = "#",
  col_types = cols(.default = col_character()),
  progress = FALSE
)
message(sprintf("Read %d assemblies from %s", nrow(raw), basename(summary_path)))

asm <- raw |>
  transmute(
    taxid         = .data$X6,
    species_taxid = .data$X7,
    genome_size   = suppressWarnings(as.numeric(.data$X26))
  ) |>
  filter(!is.na(.data$genome_size), .data$genome_size > 0)
message(sprintf("  %d with a usable genome_size", nrow(asm)))

# One row per (taxid, assembly size), counting each assembly at both its own
# taxid and its species.
per_taxid <- bind_rows(
  select(asm, taxid, genome_size),
  asm |> filter(.data$species_taxid != .data$taxid) |>
    transmute(taxid = .data$species_taxid, genome_size = .data$genome_size)
) |>
  group_by(.data$taxid) |>
  summarise(length = median(.data$genome_size), n_assemblies = n(), .groups = "drop")
message(sprintf("  %d taxids with a median genome size", nrow(per_taxid)))

fallback <- read_tsv(
  fallback_path,
  col_names = c("taxid", "length"),
  col_types = cols(taxid = col_character(), length = col_double()),
  progress = FALSE
) |>
  filter(!.data$taxid %in% per_taxid$taxid)
message(sprintf(
  "  %d taxids fall back to the summed length (median %.3f Mb, %d above 20 Mb)",
  nrow(fallback),
  median(fallback$length) / 1e6,
  sum(fallback$length > 20e6)
))
if (any(fallback$length > 20e6)) {
  big <- fallback |> filter(.data$length > 20e6) |> arrange(desc(.data$length))
  message("  fallbacks above 20 Mb are multi-assembly and still overstated:")
  message(paste(sprintf("    %-12s %8.1f Mb", big$taxid, big$length / 1e6), collapse = "\n"))
}

out <- bind_rows(
  select(per_taxid, taxid, length),
  select(fallback, taxid, length)
) |>
  arrange(.data$taxid)

if (any(out$length <= 0)) stop("a non-positive length reached the output")
if (anyDuplicated(out$taxid) > 0) stop("a taxid appears twice in the output")

# A bacterium over 20 Mb is a red flag, so report the tail the way the old
# script did. It should now be small and dominated by real eukaryotes.
message(sprintf(
  "%d taxids; length quantiles (Mb): %s",
  nrow(out),
  paste(round(quantile(out$length, c(0, .25, .5, .75, 1)) / 1e6, 3), collapse = " ")
))
message(sprintf(
  "  %d above 20 Mb (%.2f%% of taxids)",
  sum(out$length > 20e6), 100 * mean(out$length > 20e6)
))

write_tsv(out, out_path, col_names = FALSE)
message("Wrote ", out_path)
