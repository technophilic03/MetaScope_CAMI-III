# Build the taxid -> genome length table that code/metascope_to_cami.R needs.
#
# CAMI scores a profile as genome (cell) abundance, but MetaScope reports a read
# fraction. Dividing by genome length converts one into the other, so the length
# must describe the genomes the reads could actually align to.
#
# The bowtie2 indices are that exact set of sequences, and accessionTaxa.sql is
# the exact accession-to-taxid map metascope_id() used. Taking both from the same
# files is what makes the taxids line up with the MetaScope output.
#
# WHY EVERY LIBRARY, NOT ONE
# align_target_bowtie() is given every basename in libraries.txt, so a read can
# land in any of them. A length table built from one index would divide those
# read counts by a fraction of the reference and inflate the taxa it covers.
#
# This matters most for plasmids. CAMI plasmid accessions resolve to their host
# species taxid, not to a catch-all, so plasmid reads are counted against the
# host. If plasmid_01 is in libraries.txt for alignment but absent here, every
# plasmid-rich species is inflated. Borreliella burgdorferi carries more than
# ten. Keep the two lists identical.
#
#   Rscript --vanilla db/build_genome_lengths.R \
#     /projectsp/f_wj183_1/work/Yaoan/CAMI-III/db/bowtie_indices \
#     /projectsp/f_wj183_1/work/Yaoan/CAMI-III/db/manifests/libraries.txt \
#     db/accessionTaxa.sql \
#     db/genome_lengths.tsv
#
# This queries a 28 GB SQLite database over roughly ten million accessions, so
# run it on a compute node, never on the login node:
#   srun --partition=p_wj183_1 --cpus-per-task=4 --mem=128G --time=12:00:00 \
#        --pty Rscript --vanilla db/build_genome_lengths.R ...

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(taxonomizr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    "Usage: build_genome_lengths.R <lib_dir> <libraries.txt> <accessionTaxa.sql> <out.tsv>"
  )
}
lib_dir <- args[1]
libraries_path <- args[2]
taxa_db <- args[3]
out_path <- args[4]

for (f in c(lib_dir, libraries_path, taxa_db)) {
  if (!file.exists(f)) stop("not found: ", f)
}

libs <- readLines(libraries_path)
libs <- libs[nzchar(libs)]
if (length(libs) == 0) stop("no library names in ", libraries_path)

# bowtie2-inspect -s prints one "Sequence-N <name> <length>" line per reference
# sequence, tab separated. The name is the whole FASTA header, description
# included ("NC_009906.1 Plasmodium vivax chromosome 1, ..."), so only its
# leading token is the accession that accessionTaxa.sql is keyed on.
read_library <- function(lib) {
  index <- file.path(lib_dir, lib)
  inspect <- system2("bowtie2-inspect", c("-s", shQuote(index)), stdout = TRUE)
  status <- attr(inspect, "status")
  if (!is.null(status) && status != 0) {
    stop("bowtie2-inspect exited ", status, " for ", index)
  }
  parts <- str_split_fixed(inspect[str_starts(inspect, "Sequence-")], "\t", 3)
  if (nrow(parts) == 0) stop("bowtie2-inspect returned no sequences for ", index)

  # An NA accession survives accessionToTaxa()'s own order check, so a name the
  # regex fails to parse comes back looking exactly like a taxid that is absent
  # from the database. Separate the two here, while the offending name is still
  # in hand. Checking per library also keeps the full headers, which are long and
  # number in the millions, out of the combined table.
  accession <- str_extract(parts[, 2], "^\\S+")
  if (anyNA(accession)) {
    bad <- parts[which(is.na(accession))[1], 2]
    stop(sprintf(
      "%s: %d of %d sequence names have no leading non-space token; first is %s",
      lib,
      sum(is.na(accession)),
      length(accession),
      encodeString(bad, quote = "\"")
    ))
  }

  out <- tibble(
    library = lib,
    accession = accession,
    length = as.numeric(parts[, 3])
  )
  message(sprintf(
    "  %-22s %9d sequences  %8.2f Gbp",
    lib, nrow(out), sum(out$length) / 1e9
  ))
  out
}

message("Reading sequence lengths from ", length(libs), " bowtie2 indices")
seqs <- bind_rows(lapply(libs, read_library))
message(sprintf(
  "  %-22s %9d sequences  %8.2f Gbp",
  "TOTAL", nrow(seqs), sum(seqs$length) / 1e9
))

# One accession in two libraries would have its length added twice, so its taxon
# would look twice as large and half as abundant. The manifests are disjoint by
# file, so this should never fire. Stop rather than guess which copy to drop.
dups <- seqs |>
  count(accession, name = "n") |>
  filter(n > 1)
if (nrow(dups) > 0) {
  offender <- filter(seqs, accession == dups$accession[1])
  stop(sprintf(
    "%d accessions appear in more than one library; %s is in %s",
    nrow(dups),
    dups$accession[1],
    str_c(offender$library, collapse = ", ")
  ))
}

message("Resolving accessions to taxids")
probe <- head(seqs$accession, 20)
if (all(is.na(accessionToTaxa(probe, taxa_db)))) {
  stop(
    "none of the first 20 accessions resolve against ", taxa_db, "\n",
    "  first three: ", str_c(head(probe, 3), collapse = ", "), "\n",
    "  clean-looking accessions here mean the database is incomplete ",
    "(nucl_wgs was never loaded -- see CLAUDE.md)"
  )
}

seqs$taxid <- as.character(accessionToTaxa(seqs$accession, taxa_db))
unresolved <- sum(is.na(seqs$taxid))
missing_bases <- sum(seqs$length[is.na(seqs$taxid)]) / sum(seqs$length)
message(sprintf(
  "  %d of %d accessions have no taxid (%.3f%% of reference bases)",
  unresolved,
  nrow(seqs),
  100 * missing_bases
))

# Contigs and plasmids of one genome share a taxid, so the sum over a taxid is
# that genome.
lengths_tbl <- seqs |>
  filter(!is.na(taxid)) |>
  group_by(taxid) |>
  summarise(length = sum(length), .groups = "drop") |>
  arrange(taxid)

# A taxid carrying several assemblies sums to several genomes and makes that
# taxon look too big, which then reads as too few cells. Bacterial genomes stop
# well under 20 Mb, so this tail measures how far the full RefSeq is from one
# assembly per taxid. Read the numbers before trusting the profile.
oversize <- filter(lengths_tbl, length > 20e6)
message(sprintf(
  "  %d taxids; length quantiles (Mb): %s",
  nrow(lengths_tbl),
  str_c(
    round(quantile(lengths_tbl$length, c(0, .25, .5, .75, 1)) / 1e6, 3),
    collapse = " "
  )
))
message(sprintf(
  "  %d taxids above 20 Mb (%.2f%% of taxids, %.2f%% of assigned bases)",
  nrow(oversize),
  100 * nrow(oversize) / nrow(lengths_tbl),
  100 * sum(oversize$length) / sum(lengths_tbl$length)
))
if (nrow(oversize) > 0) {
  worst <- head(arrange(oversize, desc(length)), 10)
  message("  largest:")
  message(str_c(
    sprintf("    %-12s %8.1f Mb", worst$taxid, worst$length / 1e6),
    collapse = "\n"
  ))
}

write_tsv(lengths_tbl, out_path, col_names = FALSE)
message("Wrote ", out_path)
