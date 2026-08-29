# Build the taxid -> genome length table that code/metascope_to_cami.R needs.
#
# CAMI scores a profile as genome (cell) abundance, but MetaScope reports a read
# fraction. Dividing by genome length converts one into the other, so the length
# must describe the genomes the reads could actually align to.
#
# The bowtie2 index is that exact set of sequences, and accessionTaxa.sql is the
# exact accession-to-taxid map metascope_id() used. Taking both from the same
# two files is what makes the taxids line up with the MetaScope output.
#
#   Rscript --vanilla db/build_genome_lengths.R \
#     /home/yl2800/wejlab/reflib/2025_reference/bowtie_indices/target_reference \
#     db/accessionTaxa.sql \
#     db/genome_lengths.tsv
#
# This queries a 28 GB SQLite database over roughly a million accessions, so run
# it on a compute node, never on the login node:
#   srun --partition=p_wj183_1 --cpus-per-task=4 --mem=32G --time=04:00:00 \
#        --pty Rscript --vanilla db/build_genome_lengths.R ...

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(taxonomizr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: build_genome_lengths.R <bowtie2 index basename> <accessionTaxa.sql> <out.tsv>")
}
index <- args[1]
taxa_db <- args[2]
out_path <- args[3]

if (!file.exists(taxa_db)) {
  stop("accessionTaxa.sql not found: ", taxa_db)
}

# bowtie2-inspect -s prints one "Sequence-N <name> <length>" line per reference
# sequence. bowtie2 cuts a FASTA header at the first space, so the name is the
# bare accession, which is what accessionToTaxa() expects.
message("Reading sequence lengths from the bowtie2 index")
inspect <- system2("bowtie2-inspect", c("-s", shQuote(index)), stdout = TRUE)
parts <- str_split_fixed(inspect[str_starts(inspect, "Sequence-")], "\t", 3)
seqs <- tibble(accession = parts[, 2], length = as.numeric(parts[, 3]))
if (nrow(seqs) == 0) {
  stop("bowtie2-inspect returned no sequences for ", index)
}
message(sprintf("  %d sequences, %.2f Gbp total", nrow(seqs), sum(seqs$length) / 1e9))

message("Resolving accessions to taxids")
seqs$taxid <- as.character(accessionToTaxa(seqs$accession, taxa_db))
unresolved <- sum(is.na(seqs$taxid))
message(sprintf(
  "  %d of %d accessions have no taxid (%.3f%% of reference bases)",
  unresolved,
  nrow(seqs),
  100 * sum(seqs$length[is.na(seqs$taxid)]) / sum(seqs$length)
))

# Contigs of one genome share a taxid, so the sum over a taxid is that genome.
lengths_tbl <- seqs |>
  filter(!is.na(taxid)) |>
  group_by(taxid) |>
  summarise(length = sum(length), .groups = "drop") |>
  arrange(taxid)

# A taxid carrying several assemblies would sum to several genomes and make that
# taxon look too big. Bacterial genomes stop well under 20 Mb, so a long tail
# here means the reference is not one genome per taxid.
oversize <- filter(lengths_tbl, length > 20e6)
message(sprintf(
  "  %d taxids; length quantiles (Mb): %s",
  nrow(lengths_tbl),
  str_c(round(quantile(lengths_tbl$length, c(0, .25, .5, .75, 1)) / 1e6, 3),
        collapse = " ")
))
message(sprintf("  %d taxids above 20 Mb", nrow(oversize)))
if (nrow(oversize) > 0) {
  message("  largest: ", str_c(head(oversize$taxid[order(-oversize$length)], 10),
                               collapse = ", "))
}

write_tsv(lengths_tbl, out_path, col_names = FALSE)
message("Wrote ", out_path)
