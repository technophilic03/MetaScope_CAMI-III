# Restart point for samples whose alignment and EM already succeeded.
# Job 61290578 lost seven tasks to the OOM killer inside metascope_blast(); the
# .updated.bam and metascope_id.csv they produced are intact, so only the BLAST
# stage needs repeating.
now <- Sys.time()

suppressPackageStartupMessages({
    library(MetaScope)
    library(tidyverse)
    })

args <- commandArgs(trailingOnly = TRUE)

output_id_csv <- args[1]
bam_path <- args[2]
expTag <- args[3]
outDir <- args[4]
tmpDir <- args[5]
# metascope_blast() does `if (!is.numeric(num_threads)) num_threads <- 1`, so a
# character here silently drops BLAST to one core.
threads <- as.integer(args[6])
stopifnot(!is.na(threads), threads >= 1)
taxDB <- args[7]
blast_db <- args[8]

stopifnot(file.exists(output_id_csv), file.exists(bam_path))

# num_threads sets BiocParallel::MulticoreParam(workers = num_threads) in
# metascope_blast(). Each forked worker rebuilds names(all_fastas) — 16.1M
# strings for these samples — so peak memory scales with the worker count, not
# with the sample. At 32 workers seven of nineteen tasks passed 320G and died.
message("running MetaBlast")
metascope_blast(metascope_id_path = output_id_csv,
                bam_file_path = bam_path,
                tmp_dir = tmpDir,
                out_dir = outDir,
                sample_name = expTag,
                fasta_dir = tmpDir,
                num_results = 100,
                num_reads = 50,
                hit_list = 10,
                num_threads = threads,
                db_path = blast_db,
                quiet = FALSE,
                db = "ncbi",
                accession_path = taxDB)
metascope_blast_path <- file.path(outDir, paste0(expTag, ".metascope_blast.csv"))
message("MetaBlast COMPLETE")

message("Running blast reassignment")
blast_reassignment(metascope_blast_path,
                   species_threshold = 0.5,
                   num_hits = 100,
                   blast_tmp_dir = paste0(tmpDir, "/blast"),
                   out_dir = outDir,
                   sample_name = expTag)
message("BLAST REASSIGNMENT COMPLETE")

message(capture.output(Sys.time() - now))
