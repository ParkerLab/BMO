#!/bin/bash

# midpointBed
# Returns the midpoints from bed file
#
# Usage `z(cat) file.bed(.gz) | midpointBed
#
#
# Author: Ricardo D'Oliveira Albanus

cat "-" | perl -e 'while (<>) {chomp; @d=split/\t/, $_; $mid = int ($d[1] + ($d[2] - $d[1])/2); print "$d[0]\t$mid\t", $mid+1, "\t", join("\t", splice @d,3), "\n";}'
#$d[3]\n";}'