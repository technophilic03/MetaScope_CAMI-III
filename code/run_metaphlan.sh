#!/bin/bash -l

# Profile one CAMI II strain-madness sample with MetaPhlAn 4.2.
#
# Usage:
#   run_metaphlan.sh readsTar sampleName sampleId outDir scratchDir threads dbDir index
#
#   sampleName   file stem for the outputs, e.g. strmgCAMI2_sample_0
#   sampleId     the @SampleID OPAL pairs on, e.g. strmgCAMI2_short_read_sample_0.
#                --sample_id writes it into the profile and metaphlan_to_cami.R
#                reads it back from there, so there is one source of truth.
#   dbDir/index  the database directory and the pinned index name, e.g.
#                mpa_vJan25_CHOCOPhlAnSGB_202503
#
# Outputs, in outDir:
#   <sampleName>.metaphlan.tsv   the relative-abundance profile
#   <sampleName>.mapout.bz2      the bowtie2 hits, so a re-run with different
#                                profiling options can skip the alignment
#
# Choices, each against the MetaPhlAn default only where the benchmark needs it:
#
#   No R1/R2 split. MetaPhlAn does not use pairing except when it subsamples,
#   and this does not subsample, so the interleaved CAMI file goes in as is.
#
#   --skip_unclassified_estimation. Without it 4.2 rescales every clade to an
#   estimated unclassified fraction and adds an "UNCLASSIFIED -1" row, which the
#   gold standard has no counterpart for.
#
#   --offline with a pinned index. The default index is "latest", which is
#   resolved against the server's mpa_latest file, and that file already names
#   mpa_vJan26_CHOCOPhlAnSGB_202605. Unpinned, the database would depend on the
#   day the job ran.

set -euo pipefail

if [[ $# -ne 8 ]]; then
    echo "usage: $0 readsTar sampleName sampleId outDir scratchDir threads dbDir index" >&2
    exit 1
fi

readsTar="$1"
sampleName="$2"
sampleId="$3"
outDir="$4"
scratchDir="$5"
threads="$6"
dbDir="$7"
index="$8"

[[ -s "$readsTar" ]] || { echo "MISSING OR EMPTY INPUT: $readsTar" >&2; exit 1; }
for ext in 1 2 3 4 rev.1 rev.2; do
    [[ -s "${dbDir}/${index}.${ext}.bt2l" ]] \
        || { echo "MISSING INDEX FILE: ${dbDir}/${index}.${ext}.bt2l" >&2; exit 1; }
done
[[ -s "${dbDir}/${index}.pkl" ]] || { echo "MISSING: ${dbDir}/${index}.pkl" >&2; exit 1; }
command -v metaphlan >/dev/null || { echo "metaphlan is not on PATH" >&2; exit 1; }

mkdir -p "$outDir"
rm -rf "$scratchDir"
mkdir -p "$scratchDir"

profileOut="${outDir}/${sampleName}.metaphlan.tsv"
mapOut="${outDir}/${sampleName}.mapout.bz2"
readsFq="${scratchDir}/${sampleName}.fq.gz"

now=$SECONDS
echo "PROFILING ${sampleName} AS SAMPLE ${sampleId}"
# Not "metaphlan --version": that also probes the database server, and on a
# compute node with no route out it hung for 7 minutes instead of printing.
# The interpreter beside the metaphlan script is the conda env's own; a bare
# "python" can resolve to another env earlier on PATH.
"$(dirname "$(command -v metaphlan)")/python" \
    -c 'import importlib.metadata as m; print("MetaPhlAn", m.version("metaphlan"))'

# The tarball holds one interleaved, gzipped FASTQ. Copy that member out as is;
# MetaPhlAn reads gzip itself.
tar xzOf "$readsTar" --wildcards '*/anonymous_reads.fq.gz' > "$readsFq"
[[ -s "$readsFq" ]] || { echo "NOTHING EXTRACTED FROM $readsTar" >&2; exit 1; }

# MetaPhlAn refuses to overwrite an existing map file, so a re-run starts clean.
rm -f "$mapOut" "$profileOut"

metaphlan "$readsFq" \
    --input_type fastq \
    --db_dir "$dbDir" \
    --index "$index" \
    --offline \
    --mapout "$mapOut" \
    --tmp_dir "$scratchDir" \
    --nproc "$threads" \
    --sample_id "$sampleId" \
    --skip_unclassified_estimation \
    -t rel_ab \
    -o "$profileOut"

### ASSERT THE OUTPUT
# The database line, the sample id and the column header are what the converter
# relies on, so a renamed column or a different database fails here rather than
# feeding the converter the wrong field.
[[ -s "$profileOut" ]] || { echo "EMPTY OUTPUT: $profileOut" >&2; exit 1; }

gotDb=$(head -1 "$profileOut")
[[ "$gotDb" == "#${index}" ]] \
    || { echo "WRONG DATABASE LINE IN $profileOut: got '${gotDb}', want '#${index}'" >&2; exit 1; }

gotId=$(awk -F'\t' '/^#SampleID/ { print $2; exit }' "$profileOut")
[[ "$gotId" == "$sampleId" ]] \
    || { echo "WRONG SAMPLE ID IN $profileOut: got '${gotId}', want '${sampleId}'" >&2; exit 1; }

wantHeader=$'#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species'
grep -qxF "$wantHeader" "$profileOut" \
    || { echo "COLUMN HEADER CHANGED IN $profileOut" >&2; exit 1; }

# Kingdom rows are the whole classified community, so they sum to 100.
awk -F'\t' '
    !/^#/ && $1 ~ /^k__[^|]*$/ { s += $3; n++ }
    END {
        if (n == 0) { print "NO KINGDOM ROWS" > "/dev/stderr"; exit 1 }
        d = s - 100; if (d < 0) d = -d
        if (d > 1e-3) { printf "KINGDOM ROWS SUM TO %.6f, NOT 100\n", s > "/dev/stderr"; exit 1 }
        printf "%d kingdom rows, summing to %.6f\n", n, s
    }' "$profileOut"

awk -F'\t' '!/^#/ && NF { n++ } END { printf "%d clade rows\n", n }' "$profileOut"

# The FASTQ copy is large and fully reproducible from the tarball.
rm -rf "$scratchDir"

echo "PROFILE STEP COMPLETE after $((SECONDS - now))s"
