#!/bin/bash

## Count co-occurring motifs
## Usage: count_co-occuring_motifs.sh bedfile.bed maxDist

# This script will calculate the number of co-occurring motifs in a bed file,
# constrained by the maximum distance allowed. Does not work with gzipped files.

f="$1"
maxDist="$2"

midpointBed="${BMO_DIR}/src/midpointBed.sh"

tmp_file="`date +%F_%N`.tmp"

# 1) center on motif midpoint
# 2) add extentions
# 3) remove negative values (in case we extend to < 0 coordinates)
# 4) regular bedtools intersection with -c flag
cat ${f} | ${midpointBed} |\
    awk "{OFS=\"\t\"; print \$1,\$2 - ${maxDist},\$3 + ${maxDist}}" |\
    sed 's/\t-[0-9]*\t/\t0\t/g' |\
    intersectBed -F 1.0 -c -a stdin -b ${f} | cut -f 4 > ${tmp_file}

paste ${f} ${tmp_file} && rm ${tmp_file}