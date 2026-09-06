#!/bin/bash -l
#SBATCH --job-name=blast_build
#SBATCH --partition=p_wj183_1
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=yl2800@rutgers.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/logs/%x_%A_%a.out
#SBATCH --error=/projectsp/f_wj183_1/work/Yaoan/CAMI-III/output/logs/%x_%A_%a.err

# Build one BLAST nucleotide database per library in db/manifests/libraries.txt,
# from the same FASTA files db/build_bowtie_indices.sbatch turned into bowtie2
# indices.
#
#   sbatch --array=1-$(wc -l < db/manifests/libraries.txt)%4 db/build_blast_db.sbatch
#
# Then join the per-library databases into one name that blastn can open:
#
#   cd /projectsp/f_wj183_1/work/Yaoan/CAMI-III/db/blast_db
#   blastdb_aliastool -dbtype nucl -title "CAMI III RefSeq" -out cami_refseq \
#     -dblist "$(paste -sd' ' ../manifests/libraries.txt)"
#
# and point args[12] of code/run_metascope_b.R at that alias:
#
#   blast_db="/projectsp/f_wj183_1/work/Yaoan/CAMI-III/db/blast_db/cami_refseq"
#
# WHY BUILD THIS AT ALL
# The run so far searched /home/yl2800/wejlab/reflib/2024_blast_nt/nt. That is
# NCBI nt from 2024, and it holds taxa the CAMI III reference does not. BLAST can
# therefore move a read onto a taxon MetaScope could never report, and CAMI then
# scores that as a false positive of the profiler rather than of the database.
# Searching the same sequences the reads were aligned to removes that.
# The cost is that BLAST stops being an outside opinion: it can no longer show
# that a read came from something absent from the target index.
#
# WHY ONE DATABASE PER LIBRARY
# The whole reference is 652.76 Gbp and makeblastdb is single threaded, so one
# database would be a single serial job of several days with no restart point.
# One per library builds them in parallel, lets a failed library be rebuilt
# alone, and matches the bowtie2 layout, so a difference between the two steps
# points at one library. blastdb_aliastool then presents them as one database.
#
# WHY A TAXID MAP IS NOT OPTIONAL
# metascope_blast() asks blastn for the staxid column. blastn reads staxid out of
# the database, never out of accessionTaxa.sql. A database built without
# -taxid_map returns staxid 0 for every hit and the taxonomy joins downstream all
# fail. db/blast_taxid_map.R writes that map, and says why it uses
# accessionTaxa.sql rather than assembly_summary.txt.
#
# DISK
# Expect roughly 200-250 GB for the finished databases: 652.76 Gbp at two bits a
# base is 163 GB, and -parse_seqids adds accession indexes on top. Nothing large
# is staged on the way, because makeblastdb reads the FASTA from a pipe.

set -euo pipefail

stem=/projectsp/f_wj183_1/work/Yaoan/CAMI-III

manifest_dir="${stem}/db/manifests"
db_dir="${stem}/db/blast_db"
map_dir="${db_dir}/taxid_maps"
taxa_db="${stem}/db/accessionTaxa.sql"

export MODULEPATH=$MODULEPATH:/projects/community/modulefiles
module purge
module load blast
module load R/4.5.2
module list

export R_LIBS_USER="/home/yl2800/wejlab/work/Yaoan/R-4.5"

for tool in makeblastdb blastdbcmd Rscript; do
    command -v "$tool" >/dev/null || {
        echo "$tool not on PATH" >&2; exit 1; }
done

lib=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${manifest_dir}/libraries.txt")
[[ -n "$lib" ]] || { echo "no library on line ${SLURM_ARRAY_TASK_ID}" >&2; exit 1; }
manifest="${manifest_dir}/${lib}.manifest"

for f in "$manifest" "$taxa_db" "${stem}/db/blast_taxid_map.R"; do
    [[ -e "$f" ]] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done

mkdir -p "$db_dir" "$map_dir"

echo "Library: $lib"
echo "Files in manifest: $(wc -l < "$manifest")"
df -h "$db_dir"

# One decompression of the whole library. Gzip members concatenate legally, so
# the manifest can be fed through a single gzip. A 21458-name command line would
# overrun ARG_MAX, which is why the names go through the pipe and not through
# makeblastdb's -in.
stream_fasta() {
    while read -r f; do cat "$f"; done < "$manifest" | gzip -dc
}


### STEP 1  LIST THE ACCESSIONS ##############################################
#
# The accession is the first token of the header with the ">" removed, which is
# also the seq-id -parse_seqids will store. The two must be the same string or
# makeblastdb matches nothing in the map.

acc_file="${TMPDIR:-/tmp}/${lib}.acc"

# pipefail turns both a truncated gzip and a library with no header at all into a
# bare exit 1, so name the two possibilities here rather than leaving the log
# blank. Letting the failure through matters: a short accession list would set a
# low n_in, and step 4 would then compare it against an equally short database
# and call the pair correct.
echo "Step 1: reading FASTA headers"
stream_fasta | LC_ALL=C grep '^>' | LC_ALL=C sed 's/^>//; s/[[:space:]].*//' > "$acc_file" || {
    echo "READING HEADERS FAILED for $manifest: a truncated gzip, or no '>' line at all" >&2
    exit 1; }

n_in=$(wc -l < "$acc_file")
echo "  $n_in sequences"

# -parse_seqids refuses a repeated accession within one database: makeblastdb
# stops with "Duplicate seq_ids are found" after it has already read the whole
# library. Find them now instead, and name them. A repeat across two libraries is
# fine and expected, because each library is its own database; RefSeq lists some
# plasmids both in their parent assembly and in the plasmid collection.
#
# There is no right copy to drop, so this stops rather than deduplicating.
dupes=$(LC_ALL=C sort "$acc_file" | LC_ALL=C uniq -d)
if [[ -n "$dupes" ]]; then
    echo "DUPLICATE ACCESSIONS IN ${lib}: $(printf '%s\n' "$dupes" | wc -l); first three:" >&2
    printf '%s\n' "$dupes" | head -3 >&2
    exit 1
fi


### STEP 2  RESOLVE THEM TO TAXIDS ###########################################

echo "Step 2: resolving accessions to taxids"
Rscript --vanilla "${stem}/db/blast_taxid_map.R" \
    "$acc_file" "$taxa_db" "${map_dir}/${lib}.taxid_map"


### STEP 3  BUILD THE DATABASE ###############################################
#
# makeblastdb takes its FASTA from standard input, which is its documented
# default for -in. It splits the output into volumes on its own and writes a
# .nal alias over them, so ${lib} is one name to blastn however many volumes it
# ends up with.

echo "Step 3: building the database"
stream_fasta | makeblastdb \
    -dbtype nucl \
    -blastdb_version 5 \
    -parse_seqids \
    -taxid_map "${map_dir}/${lib}.taxid_map" \
    -title "CAMI III RefSeq ${lib}" \
    -out "${db_dir}/${lib}"


### STEP 4  COUNT THE SEQUENCES BACK OUT #####################################
#
# A database cut short by a full disk still opens and still answers queries, and
# the missing sequences would read as a MetaScope result rather than as a broken
# reference. Count them.

n_out=$(blastdbcmd -db "${db_dir}/${lib}" -info \
        | LC_ALL=C sed -n 's/^[[:space:]]*\([0-9,]*\) sequences;.*/\1/p' \
        | tr -d ,)
echo "Sequences in: $n_in"
echo "Sequences out: $n_out"
[[ "$n_out" -eq "$n_in" ]] || {
    echo "SEQUENCE COUNT CHANGED: $n_in in, $n_out out" >&2; exit 1; }


### STEP 5  CHECK ONE TAXID CAME THROUGH #####################################
#
# The count above says nothing about taxids, and a database whose taxids are all
# 0 fails silently: blastn still returns hits, and metascope_blast() drops them
# at the first join. Ask the database for one accession and compare.

read -r probe_acc probe_taxid < "${map_dir}/${lib}.taxid_map"
db_taxid=$(blastdbcmd -db "${db_dir}/${lib}" -entry "$probe_acc" -outfmt '%T')
echo "Taxid check: $probe_acc map=$probe_taxid db=$db_taxid"
[[ "$db_taxid" == "$probe_taxid" ]] || {
    echo "TAXID NOT STORED: $probe_acc is $probe_taxid in the map, $db_taxid in the database" >&2
    exit 1; }

rm -f "$acc_file"
echo "DONE $lib"
