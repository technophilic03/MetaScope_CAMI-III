library(MetaScope)
setwd("/projectsp/f_wj183_1/work/Yaoan/CAMI-II/db")
dir.create("bt2_index", showWarnings = FALSE)
mk_bowtie_index(ref_dir   = "ref_cat",
                lib_dir   = "bt2_index",
                lib_name  = "cami2_refseq",
                threads   = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
                overwrite = FALSE)
