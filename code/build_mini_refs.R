suppressPackageStartupMessages({
    library(Biostrings)
    library(Rbowtie2)
    library(taxonomizr)
})

stem <- "/projectsp/f_wj183_1/work/Yaoan/CAMI-III"
mini <- file.path(stem, "mini")

full_sql <- file.path(stem, "db", "accessionTaxa.sql")
names_dmp <- file.path(stem, "db", "names.dmp")
nodes_dmp <- file.path(stem, "db", "nodes.dmp")

sample_name <- "mini_01"
read_len <- 100L
insert <- 300L
threads <- Sys.getenv("SLURM_CPUS_PER_TASK", "4")

# phiX sits in the target and the filter library, so the filter step removes a
# known number of pairs instead of zero.
refs <- data.frame(
    accession = c("NC_000913.3", "NC_007795.1", "NC_000908.2", "NC_001422.1"),
    label = c("ecoli", "saureus", "mgenitalium", "phix"),
    pairs = c(12000L, 6000L, 2000L, 1000L),
    in_filter = c(FALSE, FALSE, FALSE, TRUE),
    stringsAsFactors = FALSE
)

fa_dir <- file.path(mini, "fasta")
idx_dir <- file.path(mini, "bowtie_indices")
blast_dir <- file.path(mini, "blast")
reads_dir <- file.path(mini, "reads")
for (d in c(fa_dir, idx_dir, blast_dir, reads_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## 1. Genomes -----------------------------------------------------------------

efetch <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                 "?db=nuccore&id=%s&rettype=fasta&retmode=text")
refs$fasta <- file.path(fa_dir, paste0(refs$label, ".fa"))
for (i in seq_len(nrow(refs))) {
    if (file.exists(refs$fasta[i]) && file.size(refs$fasta[i]) > 1000) next
    message("Downloading ", refs$accession[i])
    utils::download.file(sprintf(efetch, refs$accession[i]), refs$fasta[i],
                         quiet = TRUE)
}
if (any(file.size(refs$fasta) < 1000)) stop("A genome download came back empty")

# The BLAST taxid map keys on the FASTA id, so it must equal the accession.
fa_ids <- vapply(refs$fasta, function(f) sub("^>([^ ]+).*", "\\1", readLines(f, n = 1)),
                 character(1), USE.NAMES = FALSE)
if (!identical(fa_ids, refs$accession)) {
    stop("FASTA ids ", paste(fa_ids, collapse = ", "), " do not match the accessions")
}
message("Genome sizes (bp): ",
        paste(sprintf("%s=%d", refs$label,
                      vapply(refs$fasta, function(f) sum(width(readDNAStringSet(f))),
                             numeric(1), USE.NAMES = FALSE)), collapse = " "))

## 2. Taxonomy ----------------------------------------------------------------

refs$taxid <- taxonomizr::accessionToTaxa(refs$accession, full_sql)
if (anyNA(refs$taxid)) {
    stop("No taxid in ", full_sql, " for ",
         paste(refs$accession[is.na(refs$taxid)], collapse = ", "))
}
message("Taxids: ", paste(sprintf("%s=%d", refs$accession, refs$taxid), collapse = " "))

message("Subsetting ", nodes_dmp)
nodes_raw <- readLines(nodes_dmp)
node_id <- as.integer(sub("^([^\t]*)\t.*$", "\\1", nodes_raw, perl = TRUE, useBytes = TRUE))
node_parent <- as.integer(sub("^[^\t]*\t\\|\t([^\t]*)\t.*$", "\\1", nodes_raw,
                              perl = TRUE, useBytes = TRUE))
parent_of <- node_parent
names(parent_of) <- as.character(node_id)

walk_to_root <- function(taxid) {
    out <- integer(0)
    cur <- taxid
    while (!is.na(cur) && cur != 1L && !(cur %in% out)) {
        out <- c(out, cur)
        cur <- parent_of[[as.character(cur)]]
    }
    c(out, 1L)
}
lineage <- sort(unique(unlist(lapply(refs$taxid, walk_to_root))))
message("Lineage holds ", length(lineage), " taxids")

nodes_mini <- file.path(mini, "nodes.dmp")
writeLines(nodes_raw[node_id %in% lineage], nodes_mini)
rm(nodes_raw, node_id, node_parent, parent_of)

message("Subsetting ", names_dmp)
names_raw <- readLines(names_dmp)
name_id <- as.integer(sub("^([^\t]*)\t.*$", "\\1", names_raw, perl = TRUE, useBytes = TRUE))
names_mini <- file.path(mini, "names.dmp")
writeLines(names_raw[name_id %in% lineage], names_mini)
rm(names_raw, name_id)

# taxonomizr drops line one and demands exactly four columns per line.
a2t <- file.path(mini, "mini.accession2taxid")
writeLines(c("accession\taccession.version\ttaxid\tgi",
             sprintf("%s\t%s\t%d\t0", sub("\\..*$", "", refs$accession),
                     refs$accession, refs$taxid)), a2t)

mini_sql <- file.path(mini, "accessionTaxa_mini.sql")
if (file.exists(mini_sql)) file.remove(mini_sql)
taxonomizr::read.names.sql(names_mini, mini_sql)
taxonomizr::read.nodes.sql(nodes_mini, mini_sql)
taxonomizr::read.accession2taxid(a2t, mini_sql, vocal = FALSE, overwrite = TRUE)

round_trip <- taxonomizr::accessionToTaxa(refs$accession, mini_sql)
if (!identical(as.integer(round_trip), as.integer(refs$taxid))) {
    stop("Mini taxonomy DB returned ", paste(round_trip, collapse = ", "))
}
refs$species <- taxonomizr::getTaxonomy(refs$taxid, mini_sql, desiredTaxa = "species")[, 1]
# phiX only has to be filtered out, so it never needs a species name.
missing_species <- is.na(refs$species) & !refs$in_filter
if (any(missing_species)) {
    stop("No species name for ", paste(refs$taxid[missing_species], collapse = ", "))
}
message("Species: ", paste(refs$species, collapse = " | "))

## 3. Bowtie2 indices ---------------------------------------------------------

Rbowtie2::bowtie2_build(references = refs$fasta,
                        bt2Index = file.path(idx_dir, "mini_target"),
                        paste("--threads", threads), overwrite = TRUE)
Rbowtie2::bowtie2_build(references = refs$fasta[refs$in_filter],
                        bt2Index = file.path(idx_dir, "mini_filter"),
                        paste("--threads", threads), overwrite = TRUE)

## 4. BLAST database ----------------------------------------------------------

blast_fna <- file.path(blast_dir, "mini_nt.fna")
if (file.exists(blast_fna)) file.remove(blast_fna)
file.create(blast_fna)
file.append(blast_fna, refs$fasta)

taxid_map <- file.path(blast_dir, "mini_nt.taxid_map")
writeLines(sprintf("%s %d", refs$accession, refs$taxid), taxid_map)

blast_db <- file.path(blast_dir, "mini_nt")
status <- system2("makeblastdb",
                  c("-in", blast_fna, "-dbtype", "nucl", "-parse_seqids",
                    "-taxid_map", taxid_map, "-blastdb_version", "5",
                    "-title", "mini_nt", "-out", blast_db))
if (status != 0) stop("makeblastdb exited with status ", status)

## 5. Reads -------------------------------------------------------------------

set.seed(1)
draw_pairs <- function(fasta, n) {
    genome <- readDNAStringSet(fasta)[[1]]
    starts <- sample.int(length(genome) - insert, n, replace = TRUE)
    list(
        r1 = DNAStringSet(Views(genome, start = starts, width = read_len)),
        r2 = reverseComplement(
            DNAStringSet(Views(genome, start = starts + insert - read_len,
                               width = read_len)))
    )
}
# c() on a *named* list of DNAStringSets silently returns a list, so drop the
# names Map() attaches before concatenating.
drawn <- unname(Map(draw_pairs, refs$fasta, refs$pairs))
r1 <- do.call(c, lapply(drawn, `[[`, "r1"))
r2 <- do.call(c, lapply(drawn, `[[`, "r2"))
n_pairs <- sum(refs$pairs)
stopifnot(length(r1) == n_pairs, length(r2) == n_pairs,
          all(width(r1) == read_len), all(width(r2) == read_len))

ids <- sprintf("r%06d", seq_along(r1))
write_fastq <- function(seqs, ids, path) {
    con <- gzfile(path, "wb")
    on.exit(close(con))
    writeLines(as.vector(rbind(paste0("@", ids), as.character(seqs), "+",
                               strrep("I", read_len))), con)
}
fq1 <- file.path(reads_dir, paste0(sample_name, "_R1.fq.gz"))
fq2 <- file.path(reads_dir, paste0(sample_name, "_R2.fq.gz"))
write_fastq(r1, paste0(ids, "/1"), fq1)
write_fastq(r2, paste0(ids, "/2"), fq2)

# bowtie2 aborts on a quality line shorter than its sequence line, so check here.
check_fastq <- function(path, n) {
    lines <- readLines(path)
    if (length(lines) != 4L * n) {
        stop(path, " holds ", length(lines), " lines, expected ", 4L * n)
    }
    if (!all(lines[seq(3, length(lines), 4)] == "+")) stop(path, " has a bad third line")
    seqs <- lines[seq(2, length(lines), 4)]
    quals <- lines[seq(4, length(lines), 4)]
    if (!all(nchar(seqs) == read_len) || !all(nchar(quals) == read_len)) {
        stop(path, " has a record whose sequence and quality lengths differ")
    }
}
check_fastq(fq1, n_pairs)
check_fastq(fq2, n_pairs)

## 6. Expected answer ---------------------------------------------------------

target_pairs <- sum(refs$pairs[!refs$in_filter])
expected <- data.frame(species = refs$species[!refs$in_filter],
                       pairs = refs$pairs[!refs$in_filter])
expected$proportion <- round(expected$pairs / target_pairs, 4)
write.csv(expected, file.path(mini, "expected_profile.csv"), row.names = FALSE)

message("\nBuilt ", mini)
message("Read pairs written: ", length(r1), " (", sum(refs$pairs[refs$in_filter]),
        " phiX pairs the filter step must remove)")
message("Expected profile after filtering, over ", target_pairs, " pairs:")
print(expected, row.names = FALSE)
