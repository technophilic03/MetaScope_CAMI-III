#!/usr/bin/env Rscript

# Convert MetaScope species tables into a CAMI III taxonomic profile.
#
# CAMI defines PERCENTAGE as genome (cell) abundance, not read fraction:
# "the PERCENTAGE field should be 50 for both species and not 75 for A and 25
# for B, which would reflect the amount of sequence data (or number of reads)
# rather than the genome (or species) abundance."
#   https://cami-challenge.org/file-formats/#taxonomic-profiling
# MetaScope reports a read fraction, so every taxon is divided by its genome
# length and then renormalised. Lengths come from the reference build so that
# they describe the genomes the reads could actually align to.
#
# Usage:
#   Rscript code/metascope_to_cami.R \
#     --taxdump DIR --genome-lengths TSV --out FILE in1.csv [in2.csv ...]
#
#   --taxdump         directory holding nodes.dmp, names.dmp and merged.dmp
#   --genome-lengths  headerless TSV, taxid <TAB> genome length in bases
#   --out             output profile; one @SampleID block per input file
#   --sample-map      optional headerless TSV, input basename <TAB> sample ID
#
# Without --sample-map the sample ID is the input basename with
# ".metascope_id.csv" removed. The CAMI III toy gold standard names its samples
# "0" to "19" while MetaScope writes "sample_0.metascope_id.csv", and OPAL pairs
# a prediction to a gold standard by sample ID, so the map is needed there.
#
# MetaScope writes the literal "NA" as the taxid when accessionTaxa.sql has no
# taxid for an accession. Those rows leave the profile and the remaining taxa
# are rescaled to 100. The run reports how much read fraction that removed.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
})

# One entry per "|"-separated level of the CAMI rank string. Ranks listed
# together sit at the same depth and are alternatives to each other.
cami_levels <- list(
  c("acellular root", "cellular root", "other entries"),
  c("realm", "domain"),
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species",
  "strain"
)
ranks_header <- str_c(
  map_chr(cami_levels, str_c, collapse = ","),
  collapse = "|"
)
profile_version <- "0.9.2"
level_of_rank <- set_names(
  rep(seq_along(cami_levels), lengths(cami_levels)),
  unlist(cami_levels)
)
# Percentages are floored, not rounded, to six decimals. Flooring is what keeps
# the written file legal: a sum of floors never exceeds the floor of the sum, so
# every rank total stays at or below 100 and every parent stays at or above the
# sum of its children. Rounding can break both rules in the last digit.
percent_digits <- 6L
tol <- 1e-9

parse_args <- function(argv) {
  flags <- c(
    "--taxdump" = "taxdump",
    "--genome-lengths" = "genome_lengths",
    "--out" = "out",
    "--sample-map" = "sample_map"
  )
  opts <- list(
    taxdump = NA_character_,
    genome_lengths = NA_character_,
    out = NA_character_,
    sample_map = NA_character_,
    inputs = character(0)
  )
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
  for (need in c("taxdump", "genome_lengths", "out")) {
    if (is.na(opts[[need]])) {
      stop("Missing --", str_replace_all(need, "_", "-"), call. = FALSE)
    }
  }
  if (length(opts$inputs) == 0) {
    stop("No input .metascope_id.csv given", call. = FALSE)
  }
  absent <- opts$inputs[!file.exists(opts$inputs)]
  if (length(absent) > 0) {
    stop(
      "Input file not found: ",
      str_c(absent, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.na(opts$sample_map) && !file.exists(opts$sample_map)) {
    stop("Sample map not found: ", opts$sample_map, call. = FALSE)
  }
  opts
}

# OPAL pairs a prediction to its gold standard by sample ID, so the IDs written
# here must be the IDs the gold standard uses.
sample_ids_for <- function(inputs, map_path) {
  ids <- str_remove(basename(inputs), "\\.metascope_id\\.csv$")
  if (is.na(map_path)) {
    return(ids)
  }
  map_tbl <- read_tsv(
    map_path,
    col_names = c("file", "sample_id"),
    col_types = cols(.default = col_character()),
    progress = FALSE
  )
  mapped <- map_tbl$sample_id[match(basename(inputs), map_tbl$file)]
  if (anyNA(mapped)) {
    stop(
      "Sample map has no row for: ",
      str_c(basename(inputs)[is.na(mapped)], collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(mapped) > 0) {
    stop(
      "Sample map gives one ID to two inputs: ",
      str_c(unique(mapped[duplicated(mapped)]), collapse = ", "),
      call. = FALSE
    )
  }
  mapped
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

# Ancestor chains, root first, with taxid 1 left out. One vectorised lookup per
# depth step: indexing a three-million-name vector rescans it on every call.
lineages_of <- function(taxids, nodes) {
  chains <- as.list(taxids)
  frontier <- taxids
  for (depth in seq_len(100L)) {
    todo <- which(frontier != "1")
    if (length(todo) == 0) {
      return(chains)
    }
    parents <- nodes$parent[match(frontier[todo], nodes$taxid)]
    if (anyNA(parents)) {
      stop(
        "Broken taxonomy path at taxid ",
        str_c(frontier[todo][is.na(parents)], collapse = ", "),
        call. = FALSE
      )
    }
    frontier[todo] <- parents
    grew <- parents != "1"
    chains[todo[grew]] <- map2(parents[grew], chains[todo[grew]], ~ c(.x, .y))
  }
  stop("Taxonomy cycle: root not reached in 100 steps", call. = FALSE)
}

# Place each ancestor at its CAMI level. Ranks outside the CAMI list, such as
# "clade" or "subspecies", carry no level and drop out of the path.
path_slots <- function(chain, rank_of) {
  slots <- rep(NA_character_, length(cami_levels))
  lv <- unname(level_of_rank[unname(rank_of[chain])])
  keep <- !is.na(lv)
  if (anyDuplicated(lv[keep]) > 0) {
    stop(
      "Two ancestors share one CAMI level: ",
      str_c(chain, collapse = "|"),
      call. = FALSE
    )
  }
  slots[lv[keep]] <- chain[keep]
  slots
}

profile_for_sample <- function(
  csv_path,
  sample_id,
  lengths_tbl,
  nodes,
  sci,
  merged
) {
  ms <- read_csv(
    csv_path,
    col_types = cols(
      TaxonomyID = col_character(),
      EMProportion = col_double(),
      .default = col_guess()
    ),
    progress = FALSE
  )

  moved <- merged$new[match(ms$TaxonomyID, merged$old)]
  taxid <- if_else(is.na(moved), ms$TaxonomyID, moved)

  # MetaScope writes the literal "NA" when accessionTaxa.sql has no taxid for an
  # accession, and a taxdump newer than the reference can retire a taxid. Both
  # kinds of row carry no lineage, so they leave the profile. Their share cannot
  # be kept, because cell abundance needs a genome length they do not have.
  unnamed <- is.na(taxid) | taxid %in% c("NA", "")
  stale <- !unnamed & is.na(match(taxid, nodes$taxid))
  dropped <- unnamed | stale
  dropped_pct <- 100 * sum(ms$EMProportion[dropped])
  if (any(stale)) {
    message(sprintf(
      "%s: dropping %d taxid(s) absent from this taxdump: %s",
      sample_id,
      sum(stale),
      str_c(head(unique(taxid[stale]), 10), collapse = ", ")
    ))
  }
  taxid <- taxid[!dropped]
  em_kept <- ms$EMProportion[!dropped]
  if (length(taxid) == 0) {
    stop(sample_id, ": every row was dropped", call. = FALSE)
  }

  no_length <- taxid[is.na(match(taxid, lengths_tbl$taxid))]
  if (length(no_length) > 0) {
    stop(
      sample_id,
      ": ",
      length(no_length),
      " taxid(s) absent from the genome ",
      "length table: ",
      str_c(head(unique(no_length), 20), collapse = ", "),
      call. = FALSE
    )
  }

  # A merged taxid can collapse two MetaScope rows onto one node.
  counts <- tibble(taxid = taxid, em = em_kept) |>
    group_by(taxid) |>
    summarise(em = sum(em), .groups = "drop")

  # Read fraction to cell abundance, then rescaled to sum to 100.
  genome_length <- lengths_tbl$length[match(counts$taxid, lengths_tbl$taxid)]
  weight <- counts$em / genome_length
  counts$weight <- 100 * weight / sum(weight)

  chains <- lineages_of(counts$taxid, nodes)
  seen <- unique(flatten_chr(chains))
  rank_of <- set_names(nodes$rank[match(seen, nodes$taxid)], seen)
  slots <- map(chains, path_slots, rank_of = rank_of)

  # A node's percentage is the total weight of every input taxon beneath it.
  rows <- map2(slots, counts$weight, function(s, w) {
    lv <- which(!is.na(s))
    tibble(
      level = lv,
      node = s[lv],
      weight = w,
      taxpath = map_chr(
        lv,
        ~ str_c(replace_na(s[seq_len(.x)], ""), collapse = "|")
      )
    )
  }) |>
    list_rbind() |>
    group_by(node, level, taxpath) |>
    summarise(percentage = sum(weight), .groups = "drop") |>
    mutate(
      percentage = floor(percentage * 10^percent_digits) / 10^percent_digits
    )

  split_path <- count(rows, node) |> filter(n > 1)
  if (nrow(split_path) > 0) {
    stop(
      sample_id,
      ": one taxid reached by two paths: ",
      str_c(split_path$node, collapse = ", "),
      call. = FALSE
    )
  }

  used <- unique(flatten_chr(str_split(rows$taxpath, fixed("|"))))
  used <- used[used != ""]
  name_of <- set_names(sci$name[match(used, sci$taxid)], used)
  rows <- rows |>
    mutate(
      rank = unname(rank_of[node]),
      taxpathsn = map_chr(str_split(taxpath, fixed("|")), function(ids) {
        str_c(replace_na(unname(name_of[ids]), ""), collapse = "|")
      })
    ) |>
    arrange(level, desc(percentage))

  list(
    sample_id = sample_id,
    rows = rows,
    n_in = nrow(ms),
    n_moved = sum(!is.na(moved)),
    n_dropped = sum(dropped),
    dropped_pct = dropped_pct,
    n_taxa = nrow(counts)
  )
}

# The two numeric rules the CAMI format states, checked on what we will write.
check_profile <- function(p) {
  rows <- p$rows
  over <- rows |>
    group_by(rank) |>
    summarise(total = sum(percentage), .groups = "drop") |>
    filter(total > 100 + tol)
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
    filter(children > own + tol)
  if (nrow(bad) > 0) {
    stop(
      p$sample_id,
      ": children exceed parent at taxid(s) ",
      str_c(bad$node, collapse = ", "),
      call. = FALSE
    )
  }
  message(sprintf(
    paste(
      "%s: %d input rows, %d merged taxids, %d rows dropped (%.4f%% of reads),",
      "%d taxa, %d species rows, species total %.6f%%"
    ),
    p$sample_id,
    p$n_in,
    p$n_moved,
    p$n_dropped,
    p$dropped_pct,
    p$n_taxa,
    sum(rows$rank == "species"),
    sum(rows$percentage[rows$rank == "species"])
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
    c(1, 3, 5),
    c("taxid", "parent", "rank")
  )
  merged <- read_dmp(
    file.path(opts$taxdump, "merged.dmp"),
    c(1, 3),
    c("old", "new")
  )
  sci <- read_dmp(
    file.path(opts$taxdump, "names.dmp"),
    c(1, 3, 7),
    c("taxid", "name", "class")
  ) |>
    filter(class == "scientific name")
  message(sprintf(
    "  %d nodes, %d scientific names, %d merged taxids",
    nrow(nodes),
    nrow(sci),
    nrow(merged)
  ))

  lengths_tbl <- read_tsv(
    opts$genome_lengths,
    col_names = c("taxid", "length"),
    col_types = cols(taxid = col_character(), length = col_double()),
    progress = FALSE
  )
  if (any(lengths_tbl$length <= 0)) {
    stop("Genome length table holds a non-positive length", call. = FALSE)
  }
  # The taxonomy may have merged a taxid since the reference was indexed, so the
  # length table needs the same remap the sample tables get. A merge can put two
  # genomes on one taxid, and the median of them represents that taxon.
  lengths_tbl <- lengths_tbl |>
    mutate(taxid = coalesce(merged$new[match(taxid, merged$old)], taxid)) |>
    group_by(taxid) |>
    summarise(length = median(length), .groups = "drop")
  message(sprintf("Genome lengths for %d taxids", nrow(lengths_tbl)))

  profiles <- map2(
    opts$inputs,
    sample_ids_for(opts$inputs, opts$sample_map),
    profile_for_sample,
    lengths_tbl = lengths_tbl,
    nodes = nodes,
    sci = sci,
    merged = merged
  )
  walk(profiles, check_profile)

  dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
  writeLines(flatten_chr(map(profiles, format_block)), opts$out)
  message("Wrote ", opts$out, " with ", length(profiles), " sample block(s)")
}

main()
