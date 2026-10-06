library(taxonomizr)

db <- "/home/yl2800/wejlab/work/Yaoan/CAMI-II/db/ncbi_taxonomy"
out <- "/home/yl2800/wejlab/work/Yaoan/CAMI-II/refseq/accessionTaxa_2019_gb.sql"
accession <- "/home/yl2800/wejlab/work/Yaoan/CAMI-II/db/ncbi_taxonomy_accession2taxid"

# getAccession2taxid(outDir = accession, types='nucl_wgs')
read.names.sql(file.path(db, "names.dmp"), out, overwrite = TRUE)
read.nodes.sql(file.path(db, "nodes.dmp"), out, overwrite = TRUE)

read.accession2taxid(file.path(accession, "nucl_gb.accession2taxid.gz"), out, vocal=TRUE, overwrite=TRUE)
message("Done building accession taxa!")
