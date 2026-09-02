#!/usr/bin/env Rscript

# Convert mOTUs 4 profiles into a CAMI III taxonomic profile.
#
# mOTUs 4 dropped the CAMI writer that version 3 had, and it names taxa with
# GTDB lineages instead of NCBI taxids. OPAL scores against a gold standard
# keyed on NCBI taxids, so every GTDB lineage has to be resolved to a taxid
# before the profile can be compared.
#
# Unlike metascope_to_cami.R there is no genome length correction here. CAMI
# defines PERCENTAGE as genome (cell) abundance, and mOTUs already reports that:
# its abundances come from marker gene coverage, which is per genome copy, not
# per read. Dividing by genome length again would double-correct.
#
# Usage:
#   Rscript code/motus4_to_cami.R --taxdump DIR --out FILE in1.motus4.relab ...
#
#   --taxdump   directory holding nodes.dmp and names.dmp
#   --out       output profile; one @SampleID block per input file
#
# The inputs are the .relab files, not the .motus4 files. One "motus profile -o"
# writes both: the named file holds integer counts, and "<name>.relab" holds the
# same rows as relative abundances. A CAMI PERCENTAGE is a proportion, so the
# .relab file is what this reads, and its column must sum to 1.
#
# merged.dmp is not needed here. metascope_to_cami.R reads it because MetaScope
# hands it taxids that the taxonomy may since have retired. A name lookup can
# only ever return a taxid the dump still holds.
#
# The sample ID is the abundance column header inside each mOTUs file, which
# run_motus4.sh sets with "motus profile -n". OPAL pairs a prediction to a gold
# standard by sample ID, so that value must be the ID the gold standard uses.
#
# How a GTDB lineage becomes a taxid:
#   1. GTDB suffixes that mark a split genus or phylum are removed, so
#      "g__Clostridium_AQ" is looked up as "Clostridium".
#   2. A species written "Unknown <genus> mOTUv4.0_nnnnnn" is a placeholder for
#      an unnamed organism, not a species, so it is skipped.
#   3. Names are matched to NCBI scientific names at the matching rank, deepest
#      rank first. The first rank that matches wins, and the mOTU's abundance is
#      assigned there. A mOTU that resolves only to its genus therefore appears
#      at genus, never at species.
#   4. When one name matches several taxids, the candidates are filtered to
#      those descended from a higher GTDB name in the same lineage.
#
# Anything that resolves to nothing leaves the profile, and the run reports how
# much abundance that removed.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tidyr)
})

# One entry per "|"-separated level of the CAMI rank string. Ranks listed
# together sit at the same depth and are alternatives to each other. This must
# stay identical to metascope_to_cami.R, or OPAL sees two rank vocabularies.
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

# GTDB prefix to the NCBI rank it should match. A domain is spelled "domain" in
# a taxdump from 2024 or later and "superkingdom" before that, so accept both.
gtdb_ranks <- c(
  d = "domain",
  p = "phylum",
  c = "class",
  o = "order",
  f = "family",
  g = "genus",
  s = "species"
)
ncbi_ranks_for <- list(
  d = c("domain", "superkingdom"),
  p = "phylum",
  c = "class",
  o = "order",
  f = "family",
  g = "genus",
  s = "species"
)

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
      stop("Missing --", str_replace_all(need, "_", "-"), call. = FALSE)
    }
  }
  if (length(opts$inputs) == 0) {
    stop("No input .motus4.relab file given", call. = FALSE)
  }
  absent <- opts$inputs[!file.exists(opts$inputs)]
  if (length(absent) > 0) {
    stop("Input file not found: ", str_c(absent, collapse = ", "), call. = FALSE)
  }
  opts
}

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

# GTDB marks a split taxon by appending an uppercase tag to the name, as in
# "Bacillota_A" or "Clostridium_AQ innocuum". NCBI carries the unsplit name.
clean_gtdb_name <- function(x) {
  x |>
    str_remove_all("_[A-Z]+(?=$| )") |>
    str_trim()
}

# One row per rank of one GTDB lineage string, deepest rank last.
parse_lineage <- function(lineage) {
  parts <- str_split(lineage, fixed(";"))[[1]] |> str_trim()
  prefix <- str_sub(parts, 1L, 1L)
  name <- clean_gtdb_name(str_remove(parts, "^[a-z]__"))
  keep <- prefix %in% names(gtdb_ranks) & name != ""
  # "Unknown <genus> mOTUv4.0_nnnnnn" names an organism with no described
  # species, so it must not be looked up as one.
  keep <- keep & !(prefix == "s" & str_detect(name, "^Unknown "))
  tibble(gtdb_rank = prefix[keep], name = name[keep])
}

# Resolve every distinct GTDB lineage to one NCBI taxid, deepest rank first.
resolve_lineages <- function(lineages, nodes, sci) {
  parsed <- tibble(lineage = lineages) |>
    mutate(part = map(lineage, parse_lineage)) |>
    unnest(part)

  # Candidate taxids for each (name, rank) pair, in one join rather than one
  # lookup per lineage.
  wanted <- distinct(parsed, gtdb_rank, name)
  named_nodes <- sci |>
    inner_join(select(nodes, taxid, rank), by = "taxid")
  cand <- wanted |>
    mutate(ncbi_rank = ncbi_ranks_for[gtdb_rank]) |>
    unnest(ncbi_rank) |>
    inner_join(named_nodes, by = c("name", "ncbi_rank" = "rank")) |>
    select(gtdb_rank, name, taxid)

  parsed <- parsed |>
    left_join(cand, by = c("gtdb_rank", "name"), relationship = "many-to-many")

  # A candidate is only usable if a higher GTDB name in the same lineage is one
  # of its ancestors. Compute every candidate's chain once.
  hits <- filter(parsed, !is.na(taxid))
  chains <- lineages_of(unique(hits$taxid), nodes)
  names(chains) <- unique(hits$taxid)
  name_of <- set_names(sci$name, sci$taxid)
  ancestor_names <- map(chains, ~ unname(name_of[.x]))

  depth_of <- set_names(seq_along(gtdb_ranks), names(gtdb_ranks))

  resolve_one <- function(df) {
    df <- arrange(df, desc(depth_of[gtdb_rank]))
    for (rk in unique(df$gtdb_rank[!is.na(df$taxid)])) {
      here <- filter(df, gtdb_rank == rk, !is.na(taxid))
      ids <- unique(here$taxid)
      if (length(ids) > 1) {
        # Test the closest GTDB name above this rank first. Two homonyms often
        # share a family but rarely a genus, so a shallow name separates
        # nothing and would leave the choice to row order.
        above <- df |>
          filter(depth_of[gtdb_rank] < depth_of[[rk]]) |>
          arrange(desc(depth_of[gtdb_rank]))
        for (nm in unique(above$name)) {
          ok <- map_lgl(ids, ~ nm %in% ancestor_names[[.x]])
          if (any(ok)) {
            ids <- ids[ok]
            if (length(ids) == 1) {
              break
            }
          }
        }
      }
      return(tibble(
        taxid = ids[1],
        resolved_rank = gtdb_ranks[[rk]],
        ambiguous = length(ids) > 1
      ))
    }
    tibble(taxid = NA_character_, resolved_rank = NA_character_, ambiguous = FALSE)
  }

  # The unassigned mOTU carries "d__;p__;...;s__" and so parses to no rows at
  # all. Rejoining onto the full list keeps it visible as unresolved instead of
  # letting it vanish before it can be counted.
  parsed |>
    group_by(lineage) |>
    group_modify(~ resolve_one(.x)) |>
    ungroup() |>
    right_join(tibble(lineage = lineages), by = "lineage") |>
    mutate(ambiguous = replace_na(ambiguous, FALSE))
}

profile_for_sample <- function(path, nodes, sci, resolved) {
  header <- read_lines(path, n_max = 2L)
  if (length(header) < 2L || !str_starts(header[[1]], "#")) {
    stop(path, ": not a mOTUs 4 profile", call. = FALSE)
  }
  cols <- str_split(header[[2]], fixed("\t"))[[1]]
  if (length(cols) != 3L) {
    stop(path, ": expected 3 columns, found ", length(cols), call. = FALSE)
  }
  sample_id <- str_trim(cols[[3]])
  if (sample_id == "") {
    stop(path, ": the abundance column has no sample name", call. = FALSE)
  }

  tbl <- read_tsv(
    path,
    skip = 2L,
    col_names = c("motu", "lineage", "abundance"),
    col_types = cols(
      motu = col_character(),
      lineage = col_character(),
      abundance = col_double()
    ),
    progress = FALSE
  )

  total <- sum(tbl$abundance)
  if (abs(total - 1) > 1e-4) {
    stop(sample_id, ": abundances sum to ", signif(total, 8), ", not 1", call. = FALSE)
  }

  tbl <- tbl |>
    left_join(resolved, by = "lineage") |>
    filter(abundance > 0)

  dropped <- is.na(tbl$taxid)
  dropped_pct <- 100 * sum(tbl$abundance[dropped])
  kept <- tbl[!dropped, ]
  if (nrow(kept) == 0) {
    stop(sample_id, ": no mOTU resolved to a taxid", call. = FALSE)
  }

  # Several mOTUs can land on one node, either as strains of one species or as
  # unnamed organisms that both fall back to their genus.
  counts <- kept |>
    group_by(taxid) |>
    summarise(weight = sum(abundance), .groups = "drop") |>
    mutate(weight = 100 * weight / sum(weight))

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
    n_in = nrow(tbl),
    n_dropped = sum(dropped),
    dropped_pct = dropped_pct,
    n_ambiguous = sum(kept$ambiguous, na.rm = TRUE),
    by_rank = table(kept$resolved_rank),
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
      "%s: %d mOTUs, %d unresolved (%.4f%% of abundance), %d ambiguous names,",
      "%d taxa, %d species rows, species total %.6f%%"
    ),
    p$sample_id,
    p$n_in,
    p$n_dropped,
    p$dropped_pct,
    p$n_ambiguous,
    p$n_taxa,
    sum(rows$rank == "species"),
    sum(rows$percentage[rows$rank == "species"])
  ))
  message("  resolved at: ", str_c(names(p$by_rank), " ", p$by_rank, collapse = ", "))
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
  sci <- read_dmp(
    file.path(opts$taxdump, "names.dmp"),
    c(1, 3, 7),
    c("taxid", "name", "class")
  ) |>
    filter(class == "scientific name") |>
    select(taxid, name)
  message(sprintf("  %d nodes, %d scientific names", nrow(nodes), nrow(sci)))

  lineages <- opts$inputs |>
    map(~ read_tsv(
      .x,
      skip = 2L,
      col_names = c("motu", "lineage", "abundance"),
      col_types = cols(.default = col_character()),
      progress = FALSE
    )$lineage) |>
    flatten_chr() |>
    unique()
  message(sprintf("%d distinct GTDB lineages across %d file(s)", length(lineages), length(opts$inputs)))

  resolved <- resolve_lineages(lineages, nodes, sci)
  message(sprintf(
    "  resolved %d, unresolved %d",
    sum(!is.na(resolved$taxid)),
    sum(is.na(resolved$taxid))
  ))
  unresolved <- resolved$lineage[is.na(resolved$taxid)]
  if (length(unresolved) > 0) {
    message("  first unresolved: ", str_c(head(unresolved, 5), collapse = " | "))
  }

  profiles <- map(opts$inputs, profile_for_sample, nodes = nodes, sci = sci, resolved = resolved)
  walk(profiles, check_profile)

  dir.create(dirname(opts$out), recursive = TRUE, showWarnings = FALSE)
  writeLines(flatten_chr(map(profiles, format_block)), opts$out)
  message("Wrote ", opts$out, " with ", length(profiles), " sample block(s)")
}

main()
