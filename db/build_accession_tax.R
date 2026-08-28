library(taxonomizr)

db <- "/home/yl2800/wejlab/work/Yaoan/CAMI-III/db"
out <- "/home/yl2800/wejlab/work/Yaoan/CAMI-III/db/accessionTaxa.sql"

read.names.sql(file.path(db, "names.dmp"), out)
read.nodes.sql(file.path(db, "nodes.dmp"), out)
read.accession2taxid(file.path(db, "nucl_gb.accession2taxid.gz"), out, vocal = TRUE)
read.accession2taxid(file.path(db, "nucl_wgs.accession2taxid.gz"), out, vocal = TRUE)

message("Done building accession taxa!")
