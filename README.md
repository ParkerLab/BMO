# BMO
## Disclaimer
This is BMO!

Football is my best friend :)

## Usage
To run BMO:
```
Rscript bmo.R --f1 negativeBinomialAtac.bed --f2 motifCo-occurrences.bed -o outputDir -n motif_name
```
The output will be
1. Chromosome name
2. Start position
3. End position
4. Motif name
5. PWM score (by FIMO)
6. Strand
7. Number of ATAC-seq fragments
8. Number of co-occurring motifs
9. ATAC-seq fragments p-value
10. Co-occurring motifs p-value
11. Combined p-value
12. Adjusted p-value
13. BMO score (-log10 adjusted p-value)
14. Bound 0/1
