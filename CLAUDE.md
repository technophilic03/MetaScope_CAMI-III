# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

Benchmarking the **MetaScope** shotgun metagenomic profiler against the **CAMI III human-gut toy dataset**
(20 simulated samples with gold-standard taxonomic profiles). Not a software package — a
data + scripts working directory. There is no git repo, no test suite, and no build step.

Three top-level pieces:

| Path | Contents |
|---|---|
| `code/` | MetaScope driver (`run_metascope_metag.R`) + its job wrapper |
| `db/` | NCBI taxonomy dumps and the `taxonomizr` SQLite DB (`accessionTaxa.sql`, ~28 GB) |
| `toy_db/` | CAMI III `human-gut-toy` download (~176 GB) + extracted `taxonomic_profiles/` gold standard |

## Filesystem access constraint

Only these trees may be read or written:

- `/home/yl2800`
- `/home/yl2800/wejlab/work/Yaoan` (`wejlab` → `/projects/f_wj183_1` → `/projectsp/f_wj183_1`)
- `/home/yl2800/wejlab/reflib`

`/home/yl2800/wejlab/work/Yaoan/amagi` is not readable — skip it. The `/home/yl2800/wejlab/...`
and `/projectsp/f_wj183_1/...` spellings are the same files; R's `.libPaths()` and the shell
disagree on which one they print.

## Cluster environment (Rutgers Amarel, Slurm)

- Lab partition `p_wj183_1` (14-day limit, 6 nodes), account `general`. Also `main`/`mem`/`gpu` (3-day).
  Past jobs from this user all ran `--partition=p_wj183_1 --account=general`.
- Never run alignment or DB builds on the login node — `salloc`/`sbatch` first. The `amarel-setup`
  skill covers connection, staging, and submission conventions.
- Node-local scratch: use `$TMPDIR`; shared scratch is `/scratch/yl2800`.
- `singularity`/`apptainer` and `nextflow` are on `PATH`. Prebuilt MetaScope containers live one
  level up in `../metascope_16s*.sif`.

### R: pair the interpreter with its library (read this before running any R)

Packages here are compiled against one R version and break under another. The failure is
`undefined symbol: R_duplicateAsResizable`, which does not name the real cause. So always choose the
interpreter and the library together. There are **two** working pairings, and they are not interchangeable.

**Pairing A — R 4.5.2, used by the pipeline.** `code/run_metascope_metag.sbatch` sets it:

```bash
export MODULEPATH=$MODULEPATH:/projects/community/modulefiles
module load R/4.5.2
export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"
```

**Pairing B — R 4.6.1.** Reached by full path, with `R_LIBS_USER=/home/yl2800/wejlab/work/Yaoan/R`,
which is what `~/.Renviron` sets. Use it for scripts run outside the sbatch:

```bash
/cache/home/yl2800/r-envs/r-base-4.6.1/bin/Rscript --vanilla script.R
```

`MetaScope` must be installed in whichever library runs `code/run_metascope_metag.R` (or point at
`../MetaScope_fork`, a fork of wejlab/MetaScope at v1.99.10). Confirm before a batch run, because the
array fails 20 times otherwise:

```bash
module load R/4.5.2
R_LIBS_USER=/home/yl2800/wejlab/work/Yaoan/R-4.5 Rscript -e 'packageVersion("MetaScope")'
```

## Pipeline architecture

`code/run_metascope_metag.R` is a thin positional-argument CLI (`args[1..10]`:
read1, read2, indexDir, expTag, outDir, tmpDir, threads, comma-joined targets, comma-joined filters, taxDB)
that chains MetaScope's three stages:

1. `align_target_bowtie()` — align reads to the **target** libraries (microbial genomes) → BAM in `tmpDir`
2. `filter_host_bowtie()` — remove reads matching **filter** libraries (host/plasmid); emits `csv.gz`, not BAM
3. `metascope_id()` — EM reassignment of multi-mapping reads (`maxitsEM = 100`) against `accession_path`
   (the `accessionTaxa.sql` taxonomy DB) → final abundance table in `outDir`

Both alignment stages share one hand-tuned bowtie2 string:
`--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 9 --score-min L,0,1.7`. Keep target and filter identical —
the EM step assumes consistent alignment sensitivity.

### Reference libraries

`libs` are bowtie2 index **basenames** inside `lib_dir`, not paths. On Amarel they are in
`/home/yl2800/wejlab/reflib/2025_reference/bowtie_indices/`:

- `target_reference` — bacteria, archaea, fungi, protozoa, viral
- `filter_reference` — human, mouse, plasmid

Note these are `.bt2l` (large index) files, and the names differ from the legacy script's
`target_reference_2` / `filter_reference`.

### `code/run_metascope_metag_new.qsub` is legacy and does not run here

It is an **SGE array job from BU SCC**, not Slurm. Everything in it needs porting before use:
`SGE_TASK_ID`/`NSLOTS` → `SLURM_ARRAY_TASK_ID`/`SLURM_CPUS_PER_TASK`; `/restricted/projectnb/pathoscope/...`
paths → the Amarel paths above; `module load R/4.4.3` → the 4.6.1 interpreter; and the input glob
(`MeSS_out/fastq/*_R1.fq.gz`) → CAMI toy reads. Treat it as a specification of the intended job shape,
not as a runnable script.

## Taxonomy database (`db/`)

Built by `db/build_accession_tax.R` via `taxonomizr` from `names.dmp`, `nodes.dmp`, and the two
`accession2taxid` dumps. Source dumps are mirrored in `/home/yl2800/wejlab/reflib/2026_accession_taxa/`.

```bash
/cache/home/yl2800/r-envs/r-base-4.6.1/bin/Rscript --vanilla db/build_accession_tax.R
```

Both `accession2taxid` dumps go to one `read.accession2taxid()` call with `overwrite = TRUE`. That is
deliberate, and it must stay that way. The default `overwrite = FALSE` refuses to append once the
`accessionTaxa` table exists, so two sequential calls silently load only `nucl_gb` and skip `nucl_wgs`
(see `db/build_accession_tax.log`: *"already contains table accessionTaxa"*). WGS accessions then
resolve to `NA` and a large fraction of hits vanish at the `metascope_id` step.

The database was rebuilt that way and works as of 2026-08-31. Verify before trusting a run — neither
accession may come back `NA`:

```r
accessionToTaxa(c("Z17240.1", "AAAA02000001.1"), "db/accessionTaxa.sql")  # -> 9606 and a taxid
```

Pass a bare accession. A `>` left on from a FASTA header makes the lookup miss and return `NA`, which
looks exactly like a missing accession.

A rebuild is multi-hour and produces a ~28 GB file, so submit it as a batch job.

## CAMI III data and evaluation

`toy_db/` holds the raw downloads listed in `CAMI3_toy_dataset_download.list` (86 URLs from
`s3.bi.denbi.de/.../cami3_toydata/human-gut-toy/`): one set of `{reads,contigs,bam,gsa}.tar.gz` per
`sample_0..19`, plus the pooled gold-standard assembly and `source_genomes.tar.gz`. **All 86 are already
downloaded**; only `taxonomic_profiles.tar.gz` has been extracted, so per-sample data still needs unpacking:

```bash
tar xzf toy_db/sample_0_reads.tar.gz -C toy_db/
```

The `.list` file has **CRLF line endings** — `wget -i` tolerates them, but any `while read` loop over it
yields paths with a trailing `\r` and every file looks missing. Strip with `${u%$'\r'}`.

`toy_db/taxonomic_profiles/taxonomic_profile_<N>.txt` is the ground truth, in **CAMI bioboxes profiling
format** (v0.9.2): `@`-prefixed headers then columns `TAXID RANK TAXPATH TAXPATHSN PERCENTAGE`, with
per-rank percentages that each sum to 100 within a rank. MetaScope output is a flat per-taxon abundance
table, so any comparison needs a conversion step to this format (the CAMI `OPAL` tool consumes it)
— none exists in this repo yet.

## Related trees outside this project

- `../MetaScope_fork/` — fork of `wejlab/MetaScope` (`origin` = technophilic03, `upstream` = wejlab);
  where MetaScope source changes belong
- `../metascope-skill/nf-core-metascopeprolifer/` — Nextflow port of the pipeline, with its own CLAUDE.md
- `/home/yl2800/wejlab/reflib/` — shared lab reference data (bowtie2 indices, taxonomy, BLAST DBs)

The `metascope-slurm` skill generates Amarel submission scripts for the **Nextflow** pipeline from SRA
accessions; it does not target the R driver in `code/`.
