#!/usr/bin/env Rscript

# Convert MetaPhlAn 4 relative-abundance profiles into a CAMI II taxonomic
# profile.
#
# No genome length correction, as for mOTUs: MetaPhlAn estimates abundance from
# marker gene coverage, so it already reports cells rather than reads, which is
# what CAMI defines PERCENTAGE as.
#
# Usage:
#   Rscript code/metaphlan_to_cami.R --taxdump DIR --out FILE in1.metaphlan.tsv ...
#
#   --taxdump   directory holding nodes.dmp and merged.dmp
#   --out       output profile; one @SampleID block per input file
#
# A MetaPhlAn row is most of a CAMI row already. Its NCBI_tax_id column is a
# root-to-leaf taxid path with the same "|" separator and one slot per rank,
# left empty where that rank has no NCBI taxon: that is TAXPATH. Its clade_name
# with the k__/p__/... prefixes removed is TAXPATHSN. The rank comes from the
# prefix of the last clade, never from nodes.dmp, so the 2024 NCBI rename of
# superkingdom to domain cannot reach the output.
#
# Three kinds of MetaPhlAn row do not become CAMI rows. Each is checked, not
# assumed, and counted in the per-sample report:
#
#   1. t__ rows. A t__ clade is an SGB, a species-level genome bin with no NCBI
#      taxid, and its taxid path ends in an empty slot. Its abundance is already
#      inside its parent s__ row, so the SGBs of each species must sum to that
#      species, and are then dropped. Adding them to the species would count
#      them twice.
#   2. Rows whose own taxid slot is empty: unnamed SGB species such as
#      s__GGB28965_SGB41675, and the GGB genera above them. Their abundance is
#      already inside the nearest ancestor that has a taxid.
#   3. Rows whose taxid the dump does not hold even after following merged.dmp.
#      MetaPhlAn's taxids are frozen in its database, so a taxid the taxonomy
#      has merged since is renumbered rather than dropped.
#
# Nothing is renormalised. Abundance with no taxid at a rank is left out of that
# rank, as the other converters and the gold standard do, so a rank may sum to
# less than 100.
#
# The sample ID is the #SampleID line, which run_metaphlan.sh sets with
# --sample_id. OPAL pairs a prediction to a gold standard by sample ID.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
})

# The CAMI II eight-rank vocabulary, which the strain-madness gold standard
# uses. MetaPhlAn's k__ level is the domain level, which CAMI II calls
# superkingdom.
cami_ranks <- c(
  "superkingdom", "phylum", "class", "order",
  "family", "genus", "species", "strain"
)
ranks_header <- str_c(cami_ranks, collapse = "|")
profile_version <- "0.9.1"
rank_of_prefix <- c(
  k = "superkingdom",
  p = "phylum",
  c = "class",
  o = "order",
  f = "family",
  g = "genus",
  s = "species"
)
# What nodes.dmp calls each rank, for the consistency count only. The domain
# level is "domain" in a dump from 2024 on and "superkingdom" before.
ncbi_ranks_for <- list(
  superkingdom = c("domain", "superkingdom"),
  phylum = "phylum",
  class = "class",
  order = "order",
  family = "family",
  genus = "genus",
  species = "species"
)
# run_metaphlan.sh asserts this header too. Checking it here as well means this
# script cannot be fed a renamed or reordered column either.
expected_header <- "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species"
# MetaPhlAn writes five decimals, and each printed value is rounded, so a
# parent and the sum of its children can disagree by a few units in the fifth
# decimal. The checks allow for that and no more.
input_tol <- 1e-4
# No flooring, unlike the other converters: their percentages come out of
# divisions with many decimals, but these arrive with five, which six decimals
# print exactly.
percent_digits <- 6L

parse_args <- function(argv) {
  flags <- c("--taxdump" = "taxdump", "--out" = "out")
  opts <- list(taxdump = NA_character_, out = NA_character_, inputs = character(0))
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (key %in% names(flags)) {
      if (i == length(argv)) {
        stop("Missing value after ", key, call. = FALSE)
      }
      opts[[flags[[key]]]] <- argv[[i + 1L]]
      i <- i + 2L
    } else {
      opts$inputs <- c(opts$inputs, key)
      i <- i + 1L
    }
  }
  for (need in c("taxdump", "out")) {
    if (is.na(opts[[need]])) {
      stop("Missing --", need, call. = FALSE)
    }
  }
  if (length(opts$inputs) == 0) {
    stop("No input .metaphlan.tsv file given", call. = FALSE)
  }
  absent <- opts$inputs[!file.exists(opts$inputs)]
  if (length(absent) > 0) {
    stop("Input file not found: ", str_c(absent, collapse = ", "), call. = FALSE)
  }
  opts
}

# NCBI .dmp files separate fields with "\t|\t" and end lines with "\t|", so a
# plain tab split puts the wanted values at odd column positions.
read_dmp <- function(path, keep, new_names) {
  if (!file.exists(path)) {
    stop("Taxonomy file not found: ", path, call. = FALSE)
  }
  read_delim(
    path,
    delim = "\t",
    col_names = FALSE,
    quote = "",
    col_select = all_of(keep),
    col_types = cols(.default = col_character()),
    progress = FALSE
  ) |>
    set_names(new_names) |>
    mutate(across(everything(), str_trim))
}

# A taxid can be merged more than once across releases, so follow each chain to
# its end. Empty slots match nothing and pass through unchanged.
follow_merged <- function(taxids, merged) {
  out <- taxids
  for (step in seq_len(20L)) {
    nxt <- merged$new[match(out, merged$old)]
    moved <- !is.na(nxt)
    if (!any(moved)) {
      return(out)
    }
    out[moved] <- nxt[moved]
  }
  stop("merged.dmp has a chain longer than 20 steps", call. = FALSE)
}

read_metaphlan <- function(path) {
  lines <- read_lines(path, progress = FALSE)
  comments <- lines[str_starts(lines, "#")]
  id_line <- str_subset(comments, "^#SampleID\t")
  if (length(id_line) != 1) {
    stop(path, ": expected one #SampleID line, found ", length(id_line), call. = FALSE)
  }
  if (!expected_header %in% comments) {
    stop(path, ": column header is not ", expected_header, call. = FALSE)
  }
  body <- lines[!str_starts(lines, "#") & lines != ""]
  fields <- str_split_fixed(body, fixed("\t"), 4)
  rows <- tibble(
    clade = fields[, 1],
    taxpath = fields[, 2],
    abundance = suppressWarnings(as.numeric(fields[, 3]))
  )
  if (anyNA(rows$abundance)) {
    stop(path, ": relative_abundance is not a number in every row", call. = FALSE)
  }
  list(
    sample_id = str_remove(id_line, "^#SampleID\t"),
    database = comments[[1]],
    rows = rows
  )
}

profile_for_sample <- function(path, nodes, merged) {
  input <- read_metaphlan(path)
  sample_id <- input$sample_id
  rows <- input$rows

  clades <- str_split(rows$clade, fixed("|"))
  slots <- str_split(rows$taxpath, fixed("|"))
  uneven <- which(lengths(clades) != lengths(slots))
  if (length(uneven) > 0) {
    stop(
      sample_id,
      ": clade path and taxid path differ in length at ",
      str_c(rows$clade[uneven], collapse = ", "),
      call. = FALSE
    )
  }
  last_clade <- map_chr(clades, last)
  if (!all(str_detect(last_clade, "^[kpcofgst]__"))) {
    stop(
      sample_id,
      ": row without a k__..t__ clade, such as UNCLASSIFIED, which only a run",
      " without --skip_unclassified_estimation writes: ",
      str_c(head(rows$clade[!str_detect(last_clade, "^[kpcofgst]__")], 3), collapse = ", "),
      call. = FALSE
    )
  }
  prefix <- str_sub(last_clade, 1L, 1L)
  own_taxid <- map_chr(slots, last)

  # 1. SGBs must add up to their species before they are dropped.
  is_sgb <- prefix == "t"
  sgb_sums <- tibble(
    clade = str_remove(rows$clade[is_sgb], "\\|t__[^|]*$"),
    sgb = rows$abundance[is_sgb]
  ) |>
    group_by(clade) |>
    summarise(sgb = sum(sgb), .groups = "drop") |>
    left_join(select(rows, clade, species = abundance), by = "clade")
  orphan <- filter(sgb_sums, is.na(species))
  if (nrow(orphan) > 0) {
    stop(sample_id, ": SGB with no species row above it: ", orphan$clade[[1]], call. = FALSE)
  }
  off <- filter(sgb_sums, abs(sgb - species) > input_tol)
  if (nrow(off) > 0) {
    stop(
      sample_id,
      ": SGBs do not sum to their species at ",
      str_c(off$clade, " (", off$sgb, " vs ", off$species, ")", collapse = "; "),
      call. = FALSE
    )
  }

  # 2. Rows with no taxid of their own.
  no_taxid <- !is_sgb & own_taxid == ""

  # 3. Follow merged.dmp on every slot, then require the dump to hold the
  # result: a dead taxid inside a path would put a node into OPAL's tree that
  # no taxonomy has.
  keep <- !is_sgb & !no_taxid
  kept_slots <- slots[keep]
  flat <- unlist(kept_slots)
  mapped <- follow_merged(flat, merged)
  n_moved <- length(unique(flat[flat != mapped]))
  mapped_slots <- relist(mapped, kept_slots)
  leaf <- map_chr(mapped_slots, last)
  unresolved_leaf <- !leaf %in% nodes$taxid
  inner <- unlist(map(mapped_slots, head, -1))
  dead_inner <- unique(inner[inner != "" & !inner %in% nodes$taxid])
  if (length(dead_inner) > 0) {
    stop(
      sample_id,
      ": taxids inside a path that the dump does not hold: ",
      str_c(head(dead_inner, 5), collapse = ", "),
      call. = FALSE
    )
  }

  kept <- tibble(
    node = leaf,
    rank = unname(rank_of_prefix[prefix[keep]]),
    taxpath = map_chr(mapped_slots, str_c, collapse = "|"),
    taxpathsn = map2_chr(clades[keep], mapped_slots, function(clade, slot) {
      name <- str_replace_all(str_remove(clade, "^[a-z]__"), "_", " ")
      name[slot == ""] <- ""
      str_c(name, collapse = "|")
    }),
    percentage = rows$abundance[keep],
    unresolved = unresolved_leaf
  )

  species_lost <- function(flag) sum(rows$abundance[flag & prefix == "s"])
  unresolved_species <- sum(kept$percentage[kept$unresolved & kept$rank == "species"])

  # A merged taxid can land two rows of one rank on the same node.
  out <- kept |>
    filter(!unresolved) |>
    group_by(rank, node) |>
    summarise(
      taxpath = first(taxpath),
      taxpathsn = first(taxpathsn),
      percentage = sum(percentage),
      n = n(),
      .groups = "drop"
    ) |>
    mutate(rank = factor(rank, levels = cami_ranks)) |>
    arrange(rank, desc(percentage)) |>
    mutate(rank = as.character(rank))

  expected_rank <- ncbi_ranks_for[out$rank]
  actual_rank <- nodes$rank[match(out$node, nodes$taxid)]
  rank_mismatch <- sum(!map2_lgl(actual_rank, expected_rank, `%in%`))

  list(
    sample_id = sample_id,
    database = input$database,
    rows = select(out, node, rank, taxpath, taxpathsn, percentage),
    n_in = nrow(rows),
    n_sgb = sum(is_sgb),
    n_no_taxid = sum(no_taxid),
    no_taxid_species = species_lost(no_taxid),
    n_moved = n_moved,
    n_unresolved = sum(kept$unresolved),
    unresolved_species = unresolved_species,
    n_collapsed = sum(out$n - 1L),
    n_rank_mismatch = rank_mismatch
  )
}

# The two CAMI invariants: no rank above 100, and no parent below the sum of its
# children. Both within MetaPhlAn's own rounding.
check_profile <- function(p) {
  rows <- p$rows
  totals <- rows |>
    group_by(rank) |>
    summarise(total = sum(percentage), .groups = "drop")
  over <- filter(totals, total > 100 + input_tol)
  if (nrow(over) > 0) {
    stop(
      p$sample_id,
      ": rank total above 100%: ",
      str_c(over$rank, " = ", signif(over$total, 8), collapse = "; "),
      call. = FALSE
    )
  }
  parent_node <- map_chr(str_split(rows$taxpath, fixed("|")), function(parts) {
    above <- head(parts, -1)
    above <- above[above != ""]
    if (length(above) == 0) NA_character_ else tail(above, 1)
  })
  bad <- tibble(node = parent_node, percentage = rows$percentage) |>
    filter(!is.na(node)) |>
    group_by(node) |>
    summarise(children = sum(percentage), .groups = "drop") |>
    left_join(select(rows, node, own = percentage), by = "node") |>
    filter(children > own + input_tol)
  if (nrow(bad) > 0) {
    stop(
      p$sample_id,
      ": children exceed parent at taxid(s) ",
      str_c(bad$node, collapse = ", "),
      call. = FALSE
    )
  }
  total_of <- function(r) sum(rows$percentage[rows$rank == r])
  message(sprintf(
    paste(
      "%s: %d rows, %d SGB rows dropped (each species = sum of its SGBs),",
      "%d rows with no taxid (%.4f%% at species), %d merged taxids,",
      "%d unresolved (%.4f%% at species), %d rows collapsed, %d rank mismatches,",
      "%d taxa, superkingdom %.6f%%, species %.6f%%, strain rows %d"
    ),
    p$sample_id,
    p$n_in,
    p$n_sgb,
    p$n_no_taxid,
    p$no_taxid_species,
    p$n_moved,
    p$n_unresolved,
    p$unresolved_species,
    p$n_collapsed,
    p$n_rank_mismatch,
    nrow(rows),
    total_of("superkingdom"),
    total_of("species"),
    sum(rows$rank == "strain")
  ))
  invisible(TRUE)
}

format_block <- function(p) {
  c(
    str_c("@SampleID:", p$sample_id),
    str_c("@Version:", profile_version),
    str_c("@Ranks:", ranks_header),
    "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE",
    str_c(
      p$rows$node,
      p$rows$rank,
      p$rows$taxpath,
      p$rows$taxpathsn,
      formatC(p$rows$percentage, format = "f", digits = percent_digits),
      sep = "\t"
    )
  )
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  message("Reading taxonomy from ", opts$taxdump)
  nodes <- read_dmp(
    file.path(opts$taxdump, "nodes.dmp"),
    c(1, 5),
    c("taxid", "rank")
  )
  merged <- read_dmp(
    file.path(opts$taxdump, "merged.dmp"),
    c(1, 3),
    c("old", "new")
  )
  message(sprintf("  %d nodes, %d merged taxids", nrow(nodes), nrow(merged)))

  profiles <- map(opts$inputs, profile_for_sample, nodes = nodes, merged = merged)

  databases <- unique(map_chr(profiles, "database"))
  if (length(databases) != 1) {
    stop("Inputs come from different databases: ", str_c(databases, collapse = ", "), call. = FALSE)
  }
  message("MetaPhlAn database: ", str_remove(databases, "^#"))

  ids <- map_chr(profiles, "sample_id")
  if (anyDuplicated(ids) > 0) {
    stop("Sample ID used twice: ", str_c(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  }

  walk(profiles, check_profile)

  dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
  writeLines(flatten_chr(map(profiles, format_block)), opts$out)
  message("Wrote ", opts$out, " with ", length(profiles), " sample block(s)")
}

main()
