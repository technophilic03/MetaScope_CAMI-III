#!/bin/bash
#
# Fetch the CAMI II strain-madness short-read samples.
#
# 100 tarballs of 2,274,274,241 bytes each, 227 GB in total, from the CAMI II
# data repository at frl.publisso.de (DOI 10.4126/FRL01-006425521). No
# registration and no special client: the cami-challenge.org site still shows a
# request form for CAMI 2 data, but the repository serves these over plain
# anonymous HTTPS.
#
# Each tarball holds two members:
#   short_read/<timestamp>_sample_<N>/reads/anonymous_reads.fq.gz   the reads
#   short_read/<timestamp>_sample_<N>/reads/reads_mapping.tsv.gz    CAMISIM truth
# The reads are ONE interleaved gzipped FASTQ, 150 bp, with /1 and /2 header
# suffixes -- the same convention as CAMI III, so the splitter in
# run_metascope_b.sbatch works on them unchanged.
#
# Downloads are sequential on purpose. The host's tolerance for parallel
# connections is unknown, and this is a one-off cost; wget -c makes it resumable.
#
# Usage:
#   bash code/fetch_strmg_reads.sh            # all 100
#   bash code/fetch_strmg_reads.sh 1 99       # an inclusive index range
#
set -euo pipefail

ROOT=https://frl.publisso.de/data/frl:6425521
BASE=strain/short_read
PFX=strmgCAMI2
DEST=/projectsp/f_wj183_1/work/Yaoan/CAMI-II/strmg_db
FIRST=${1:-0}
LAST=${2:-99}

mkdir -p "$DEST"
cd "$DEST"

# The md5 manifest covers the whole repository; keep only our 100 rows, with the
# directory prefix stripped so "md5sum -c" can run from this directory.
if [[ ! -s md5s.strmg_reads ]]; then
    wget -q -O md5sums.txt "${ROOT}/md5sums.txt"
    grep "^.*  ${BASE}/${PFX}_sample_[0-9]*_reads\.tar\.gz$" md5sums.txt \
      | sed "s#  ${BASE}/#  #" > md5s.strmg_reads
fi
want=$(wc -l < md5s.strmg_reads)
[[ "$want" -eq 100 ]] || { echo "EXPECTED 100 MD5 ROWS, GOT ${want}" >&2; exit 1; }

for i in $(seq "$FIRST" "$LAST"); do
    f="${PFX}_sample_${i}_reads.tar.gz"

    # Skip anything already present and already correct, so a re-run is cheap.
    if [[ -s "$f" ]] && grep " ${f}\$" md5s.strmg_reads | md5sum -c --status - 2>/dev/null; then
        echo "HAVE ${f}"
        continue
    fi

    wget -c -nv -O "$f" "${ROOT}/${BASE}/${f}"

    grep " ${f}\$" md5s.strmg_reads | md5sum -c - \
      || { echo "MD5 MISMATCH: ${f}" >&2; exit 1; }
done

echo
have=$(for i in $(seq 0 99); do [[ -s "${PFX}_sample_${i}_reads.tar.gz" ]] && echo x; done | wc -l)
echo "TARBALLS PRESENT: ${have} / 100"
