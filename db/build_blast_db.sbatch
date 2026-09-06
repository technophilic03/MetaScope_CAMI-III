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


# Check that the tools we need are on PATH. 
for tool in makeblastdb blastdbcmd Rscript; do
    command -v "$tool" >/dev/null || {
        echo "$tool not on PATH" >&2; exit 1; }
done

# Check that the manifest, the accession-to-taxid database, and the R script exist. 
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

# Stream the FASTA files listed in the manifest, decompressing them on the fly. 
stream_fasta() {
    while read -r f; do cat "$f"; done < "$manifest" | gzip -dc
}


### STEP 1  LIST THE ACCESSIONS ##############################################

acc_file="${TMPDIR:-/tmp}/${lib}.acc"

echo "Step 1: reading FASTA headers"
stream_fasta | LC_ALL=C grep '^>' | LC_ALL=C sed 's/^>//; s/[[:space:]].*//' > "$acc_file" || {
    echo "READING HEADERS FAILED for $manifest: a truncated gzip, or no '>' line at all" >&2
    exit 1; }

n_in=$(wc -l < "$acc_file")
echo "  $n_in sequences"

# detect duplicates if any, because makeblastdb will silently drop duplicates and we want to know about them.
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

echo "Step 3: building the database"
stream_fasta | makeblastdb \
    -dbtype nucl \
    -blastdb_version 5 \
    -parse_seqids \
    -taxid_map "${map_dir}/${lib}.taxid_map" \
    -title "CAMI III RefSeq ${lib}" \
    -out "${db_dir}/${lib}"


### STEP 4  COUNT THE SEQUENCES BACK OUT #####################################

# The count of sequences in the database should match the count of sequences in the manifest. 
# If not, something went wrong (e.g., a truncated gzip file, or a duplicate accession that makeblastdb silently dropped).

n_out=$(blastdbcmd -db "${db_dir}/${lib}" -info \
        | LC_ALL=C sed -n 's/^[[:space:]]*\([0-9,]*\) sequences;.*/\1/p' \
        | tr -d ,)
echo "Sequences in: $n_in"
echo "Sequences out: $n_out"
[[ "$n_out" -eq "$n_in" ]] || {
    echo "SEQUENCE COUNT CHANGED: $n_in in, $n_out out" >&2; exit 1; }


### STEP 5  CHECK ONE TAXID CAME THROUGH #####################################

# The taxid of one accession in the map should match the taxid of that accession in the database.
read -r probe_acc probe_taxid < "${map_dir}/${lib}.taxid_map"
db_taxid=$(blastdbcmd -db "${db_dir}/${lib}" -entry "$probe_acc" -outfmt '%T')
echo "Taxid check: $probe_acc map=$probe_taxid db=$db_taxid"
[[ "$db_taxid" == "$probe_taxid" ]] || {
    echo "TAXID NOT STORED: $probe_acc is $probe_taxid in the map, $db_taxid in the database" >&2
    exit 1; }

rm -f "$acc_file"
echo "DONE $lib"
