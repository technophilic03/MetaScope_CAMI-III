# Map every accession in one CAMI III library to a taxid, for makeblastdb.
#
# metascope_blast() asks blastn for the staxid column, and blastn can only
# answer it from taxids stored inside the BLAST database. makeblastdb reads them
# from -taxid_map, which wants one "accession<TAB>taxid" line per sequence.
# Without that file every hit comes back with no staxid, and every join that
# metascope_blast() and blast_reassignment() make on staxid drops out.
#
# WHY accessionTaxa.sql AND NOT assembly_summary.txt
# assembly_summary.txt column 6 also carries a taxid, one per assembly, and it
# covers every genome with no database lookup at all. It is still the wrong
# source here. metascope_id() resolves the aligned accessions through
# accessionTaxa.sql, so a BLAST database built from the other file can label an
# accession with a strain taxid where the alignment step labelled it with the
# species. blast_reassignment() cannot tell that disagreement apart from a wrong
# assignment. One taxonomy for both steps is what makes the comparison mean
# anything.
#
#   Rscript --vanilla db/blast_taxid_map.R <accessions.txt> <accessionTaxa.sql> <out.tsv>
#
# accessions.txt holds one bare accession per line, each one distinct.
# db/build_blast_db.sbatch writes it from the FASTA headers of one manifest and
# rejects the library if an accession repeats, because -parse_seqids would then
# refuse to build at all.

suppressPackageStartupMessages({
  library(taxonomizr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: blast_taxid_map.R <accessions.txt> <accessionTaxa.sql> <out.tsv>")
}
acc_path <- args[1]
taxa_db <- args[2]
out_path <- args[3]

for (f in c(acc_path, taxa_db)) {
  if (!file.exists(f)) stop("not found: ", f)
}

accession <- readLines(acc_path)
accession <- accession[nzchar(accession)]
if (length(accession) == 0) stop("no accessions in ", acc_path)

# A ">" left on from a FASTA header makes accessionToTaxa() miss, and that miss
# looks exactly like an accession the database does not hold. Catch it here,
# while the offending string is still in hand.
if (any(startsWith(accession, ">"))) {
  stop(sum(startsWith(accession, ">")), " accessions still carry a leading '>'")
}
if (any(grepl("[[:space:]]", accession))) {
  stop(sum(grepl("[[:space:]]", accession)), " accessions contain whitespace")
}

message(sprintf("%d accessions read", length(accession)))

# An empty result and a broken database look the same downstream, so fail here
# instead. A complete accessionTaxa.sql resolves nearly every RefSeq accession.
probe <- head(accession, 20)
if (all(is.na(accessionToTaxa(probe, taxa_db)))) {
  stop(
    "none of the first 20 accessions resolve against ", taxa_db, "\n",
    "  first three: ", paste(head(probe, 3), collapse = ", "), "\n",
    "  clean-looking accessions here mean the database is incomplete ",
    "(nucl_wgs was never loaded -- see CLAUDE.md)"
  )
}

taxid <- accessionToTaxa(accession, taxa_db)
unresolved <- sum(is.na(taxid))
message(sprintf(
  "%d of %d accessions have no taxid (%.3f%%)",
  unresolved, length(taxid), 100 * unresolved / length(taxid)
))

# An unresolved accession still belongs in the BLAST database; it simply carries
# no taxid, and blastn then reports staxid 0 for it. Leaving it out of the map is
# how makeblastdb is told that.
keep <- !is.na(taxid)
out <- data.frame(accession = accession[keep], taxid = taxid[keep])
write.table(out, out_path, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
message(sprintf("Wrote %d mappings to %s", nrow(out), out_path))
