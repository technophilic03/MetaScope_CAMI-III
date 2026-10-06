# 0. Package and paths ----
library(data.table)
library(taxonomizr)

ref_fasta <- 
old_gb <- "db/ncbi_taxonomy_accession2taxid/nucl_gb.accession2taxid.gz"
cur_wgs <- "db/ncbi_taxonomy_09292026/nucl_wgs.accession2taxid.gz"
cur_dead <- "db/ncbi_taxonomy_09292026/dead_wgs.accession2taxid.gz"
