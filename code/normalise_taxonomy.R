#!/usr/bin/env Rscript

# Put CAMI profiles into one NCBI taxonomy by following merged.dmp, and keep
# only the sample blocks being compared.
#
# The CAMI II paper (Meyer et al. 2022, Nat Methods) did the same during its
# assessments: "NCBI's merged.dmp file was used to map synonymous taxa". Five
# taxonomy vintages meet in this comparison -- the gold standards, the MetaScope
# reference, the mOTUs 3 marker database (2019), mOTUs 4's name lookup and the
# MetaPhlAn index -- and a profiler that reports a taxid NCBI has since merged
# into a gold-standard taxid scores a false positive and a false negative for
# what is one organism.
#
# Usage:
#   Rscript code/normalise_taxonomy.R --taxdump DIR --out-dir DIR \
#       [--samples ID,ID,...] in1.profile in2.profile ...
#
#   --taxdump   directory holding nodes.dmp, names.dmp and merged.dmp
#   --out-dir   each input is written here under its own file name
#   --samples   keep only these @SampleID blocks, in input order; every ID
#               listed must be present in every input
#
# What changes, and nothing else does:
#   - every integer taxid in TAXID and in each TAXPATH slot is followed through
#     merged.dmp to the end of its chain
#   - a TAXPATHSN slot whose taxid moved takes the dump's scientific name, so a
#     new taxid is never paired with a stale name
#   - rows that land on the same RANK and TAXID are summed into the first
#
# What passes through untouched: headers, comments, column order, extra columns
# such as the gold standard's _CAMI_genomeID, empty TAXPATH slots, and taxids
# that are not integers -- the gold standard's strains are synthetic
# "<parent>.<N>" IDs that no dump holds. A row with nothing to change is written
# back byte for byte.
#
# Checked on every block: each rank's total is unchanged, and the number of
# CAMI rule violations (a rank above 100, children above their parent) is no
# higher than in the input. Some inputs break those rules already through
# rounding, so the check is that normalising adds none, not that there are none.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(tibble)
})

tol <- 1e-6

parse_args <- function(argv) {
  flags <- c("--taxdump" = "taxdump", "--out-dir" = "out_dir", "--samples" = "samples")
  opts <- list(
    taxdump = NA_character_,
    out_dir = NA_character_,
    samples = NA_character_,
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
  for (need in c("taxdump", "out_dir")) {
    if (is.na(opts[[need]])) {
      stop("Missing --", str_replace_all(need, "_", "-"), call. = FALSE)
    }
  }
  if (length(opts$inputs) == 0) {
    stop("No input profile given", call. = FALSE)
  }
  absent <- opts$inputs[!file.exists(opts$inputs)]
  if (length(absent) > 0) {
    stop("Input file not found: ", str_c(absent, collapse = ", "), call. = FALSE)
  }
  out <- file.path(opts$out_dir, basename(opts$inputs))
  if (any(normalizePath(out, mustWork = FALSE) == normalizePath(opts$inputs))) {
    stop("--out-dir would overwrite an input", call. = FALSE)
  }
  opts$keep <- if (is.na(opts$samples)) NULL else str_trim(str_split_1(opts$samples, ","))
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
# its end. Non-integer IDs and empty slots match nothing and pass through.
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

# A block starts at its @SampleID line. Comment and blank lines directly above
# it belong to it too: mOTUs 3 writes three "#" lines and a blank line before
# each @SampleID, and they should travel with the block they describe.
split_blocks <- function(lines) {
  starts <- which(str_starts(lines, "@SampleID:"))
  if (length(starts) == 0) {
    stop("No @SampleID line", call. = FALSE)
  }
  loose <- str_starts(lines, "#") | str_trim(lines) == ""
  first <- starts
  floor_at <- c(0L, head(starts, -1))
  for (k in seq_along(starts)) {
    s <- starts[[k]]
    while (s - 1L > floor_at[[k]] && loose[[s - 1L]]) {
      s <- s - 1L
    }
    first[[k]] <- s
  }
  last <- c(tail(first, -1) - 1L, length(lines))
  list(
    preamble = if (first[[1]] > 1L) lines[seq_len(first[[1]] - 1L)] else character(0),
    blocks = map2(first, last, ~ lines[.x:.y]),
    ids = str_trim(str_remove(lines[starts], "^@SampleID:"))
  )
}

# Rule violations, counted rather than raised: some inputs carry rounding
# violations of their own, and the test is only that normalising adds none.
violations <- function(rank, taxid, taxpath, pct) {
  over <- sum(tapply(pct, rank, sum) > 100 + tol)
  parent <- map_chr(str_split(taxpath, fixed("|")), function(parts) {
    above <- head(parts, -1)
    above <- above[above != ""]
    if (length(above) == 0) NA_character_ else tail(above, 1)
  })
  own <- tapply(pct, taxid, sum)
  children <- tapply(pct[!is.na(parent)], parent[!is.na(parent)], sum)
  shared <- intersect(names(children), names(own))
  over + sum(children[shared] > own[shared] + tol)
}

normalise_block <- function(block, merged, sci) {
  header_at <- which(str_starts(block, "@@"))
  if (length(header_at) != 1) {
    stop("Block without exactly one @@ column line: ", block[str_starts(block, "@SampleID")], call. = FALSE)
  }
  columns <- str_split_1(str_remove(block[[header_at]], "^@@"), fixed("\t"))
  col <- match(c("TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"), columns)
  if (anyNA(col)) {
    stop("Column line lacks TAXID, RANK, TAXPATH, TAXPATHSN or PERCENTAGE: ", block[[header_at]], call. = FALSE)
  }
  data_at <- which(!str_starts(block, "[@#]") & str_trim(block) != "")
  fields <- str_split(block[data_at], fixed("\t"))
  if (any(lengths(fields) < max(col))) {
    stop("Data row with fewer columns than the header", call. = FALSE)
  }
  get <- function(j) map_chr(fields, j)
  taxid <- get(col[[1]])
  rank <- get(col[[2]])
  taxpath <- get(col[[3]])
  names_sn <- get(col[[4]])
  pct <- as.numeric(get(col[[5]]))
  if (anyNA(pct)) {
    stop("PERCENTAGE is not a number in every row", call. = FALSE)
  }

  new_taxid <- follow_merged(taxid, merged)
  slots <- str_split(taxpath, fixed("|"))
  new_slots <- relist(follow_merged(unlist(slots), merged), slots)
  sn <- str_split(names_sn, fixed("|"))
  moved_row <- new_taxid != taxid | map2_lgl(slots, new_slots, ~ any(.x != .y))

  rebuilt <- block[data_at]
  for (r in which(moved_row)) {
    f <- fields[[r]]
    f[[col[[1]]]] <- new_taxid[[r]]
    f[[col[[3]]]] <- str_c(new_slots[[r]], collapse = "|")
    if (length(sn[[r]]) == length(slots[[r]])) {
      renamed <- slots[[r]] != new_slots[[r]]
      sn[[r]][renamed] <- sci$name[match(new_slots[[r]][renamed], sci$taxid)]
      f[[col[[4]]]] <- str_c(coalesce(sn[[r]], ""), collapse = "|")
    }
    rebuilt[[r]] <- str_c(f, collapse = "\t")
  }

  # Rows that now share a RANK and TAXID are one taxon. The row that already
  # carried that taxid is kept, because its path agrees with it -- a renumbered
  # mOTUs 3 recall row can bring an ambiguous path such as "(31977/541000)" --
  # and the others are summed into it and removed. If every row moved, the
  # first is kept.
  key <- str_c(rank, "\t", new_taxid)
  new_pct <- pct
  drop <- logical(length(key))
  path_conflicts <- 0L
  for (k in unique(key[duplicated(key)])) {
    rows <- which(key == k)
    unmoved <- rows[taxid[rows] == new_taxid[rows]]
    keep <- if (length(unmoved) > 0) unmoved[[1]] else rows[[1]]
    new_pct[[keep]] <- sum(pct[rows])
    path_conflicts <- path_conflicts +
      sum(map_chr(new_slots[rows], str_c, collapse = "|") != str_c(new_slots[[keep]], collapse = "|"))
    f <- str_split_1(rebuilt[[keep]], fixed("\t"))
    f[[col[[5]]]] <- as.character(new_pct[[keep]])
    rebuilt[[keep]] <- str_c(f, collapse = "\t")
    drop[setdiff(rows, keep)] <- TRUE
  }
  n_collapsed <- sum(drop)
  kept <- !drop

  before <- tapply(pct, rank, sum)
  after <- tapply(new_pct[kept], rank[kept], sum)
  if (!setequal(names(before), names(after)) || any(abs(before - after[names(before)]) > tol)) {
    stop("A rank total changed", call. = FALSE)
  }
  v_before <- violations(rank, taxid, taxpath, pct)
  v_after <- violations(
    rank[kept], new_taxid[kept],
    map_chr(new_slots[kept], str_c, collapse = "|"), new_pct[kept]
  )

  out <- block
  out[data_at] <- rebuilt
  # Guarded because x[-integer(0)] is empty, not x.
  if (n_collapsed > 0) {
    out <- out[-data_at[drop]]
  }
  list(
    lines = out,
    moved = tibble(old = taxid, new = new_taxid, rank = rank, pct = pct) |>
      filter(old != new),
    n_rows = length(data_at),
    n_rewritten = sum(moved_row),
    n_collapsed = n_collapsed,
    path_conflicts = path_conflicts,
    integer_ids = new_taxid[str_detect(new_taxid, "^[0-9]+$")],
    integer_pct = new_pct[str_detect(new_taxid, "^[0-9]+$")],
    v_before = v_before,
    v_after = v_after
  )
}

normalise_file <- function(path, out_path, keep, merged, nodes, sci) {
  parsed <- split_blocks(read_lines(path, progress = FALSE))
  chosen <- seq_along(parsed$blocks)
  if (!is.null(keep)) {
    missing <- setdiff(keep, parsed$ids)
    if (length(missing) > 0) {
      stop(basename(path), ": sample(s) not present: ", str_c(missing, collapse = ", "), call. = FALSE)
    }
    chosen <- which(parsed$ids %in% keep)
  }
  done <- map(parsed$blocks[chosen], normalise_block, merged = merged, sci = sci)

  v_before <- sum(map_int(done, "v_before"))
  v_after <- sum(map_int(done, "v_after"))
  if (v_after > v_before) {
    stop(basename(path), ": normalising added CAMI rule violations (", v_before, " -> ", v_after, ")", call. = FALSE)
  }

  writeLines(c(parsed$preamble, flatten_chr(map(done, "lines"))), out_path)

  moved <- bind_rows(map(done, "moved"))
  ids <- flatten_chr(map(done, "integer_ids"))
  pcts <- flatten_dbl(map(done, "integer_pct"))
  dead <- !ids %in% nodes$taxid
  message(sprintf(
    paste(
      "%s: %d of %d blocks kept, %d rows; %d distinct taxids merged, %d rows rewritten,",
      "%d rows collapsed (%d with differing paths); %d taxids absent from nodes.dmp",
      "(largest %.4f%%); rule violations %d before, %d after"
    ),
    basename(path), length(chosen), length(parsed$blocks),
    sum(map_int(done, "n_rows")),
    n_distinct(moved$old), sum(map_int(done, "n_rewritten")),
    sum(map_int(done, "n_collapsed")), sum(map_int(done, "path_conflicts")),
    n_distinct(ids[dead]), if (any(dead)) max(pcts[dead]) else 0,
    v_before, v_after
  ))
  if (nrow(moved) > 0) {
    detail <- moved |>
      group_by(old, new, rank) |>
      summarise(blocks = n(), mean_pct = mean(pct), .groups = "drop") |>
      arrange(desc(mean_pct))
    message(str_c(
      sprintf("    %s -> %s  %-12s in %d block(s), mean %.4f%%", detail$old, detail$new, detail$rank, detail$blocks, detail$mean_pct),
      collapse = "\n"
    ))
  }
  invisible(TRUE)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  message("Reading taxonomy from ", opts$taxdump)
  merged <- read_dmp(file.path(opts$taxdump, "merged.dmp"), c(1, 3), c("old", "new"))
  nodes <- read_dmp(file.path(opts$taxdump, "nodes.dmp"), 1, "taxid")
  sci <- read_dmp(
    file.path(opts$taxdump, "names.dmp"),
    c(1, 3, 7),
    c("taxid", "name", "class")
  ) |>
    filter(class == "scientific name") |>
    select(taxid, name)
  message(sprintf("  %d nodes, %d scientific names, %d merged taxids", nrow(nodes), nrow(sci), nrow(merged)))
  if (!is.null(opts$keep)) {
    message("Keeping ", length(opts$keep), " sample block(s): ", str_c(opts$keep, collapse = ","))
  }

  dir.create(opts$out_dir, recursive = TRUE, showWarnings = FALSE)
  walk(opts$inputs, function(path) {
    normalise_file(path, file.path(opts$out_dir, basename(path)), opts$keep, merged, nodes, sci)
  })
  message("Wrote ", length(opts$inputs), " profile(s) to ", opts$out_dir)
}

main()
