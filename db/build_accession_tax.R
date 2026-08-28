library(taxonomizr)

db <- "/home/yl2800/wejlab/work/Yaoan/CAMI-III/db"
out <- "/home/yl2800/wejlab/work/Yaoan/CAMI-III/db/accessionTaxa.sql"

read.names.sql(file.path(db, "names.dmp"), out)
read.nodes.sql(file.path(db, "nodes.dmp"), out)
read.accession2taxid(file.path(db, c("nucl_gb.accession2taxid.gz", "nucl_wgs.accession2taxid.gz")), out, vocal=TRUE, overwrite=TRUE)
message("Done building accession taxa!")
