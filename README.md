# ParkerLab BMO

# Overview
BMO (pronounced *beemo* - yes, that BMO) is an algorithm to predict TF binding from ATAC-seq data without using footprints. BMO uses negative binomial models of ATAC-seq fragments and number of co-occurring motifs to determine the likelihood of a given motif instance being bound.

**Disclaimer**: As of now, BMO has only been tested with paired-end ATAC-seq data.

# Required input
1. ATAC-seq experiment pruned BAM file\*
2. MACS2 peak calls\*\*
3. Motif BED6\*\*\* files. One per PWM scan (we use [FIMO](http://meme-suite.org/doc/fimo.html)).

\* High quality, no duplicates, and properly mapped reads only.

\** We call them with [MACS2](https://github.com/taoliu/MACS) using flags `--broad --nomodel --shift -100 --extsize 200 --keep-dup all` and then intersect out [blacklisted regions](https://sites.google.com/site/anshulkundaje/projects/blacklists).

\*\*\* We haven't tested with BED4, but it probably works.

# Running BMO
## Using Snake
We strongly recommend you run BMO using [Snakemake]( https://snakemake.readthedocs.io), an extremely powerful tool for running reproducible pipelines. In order to run BMO using Snakemake, you just need to copy both `src/Snakefile` and `input/config.yaml` to wherever you want to run the analysis in your machine/cluster and then update the config file to point to your motif BED files directory and your sample's BAM and MACS2 peak calls files. The example config file we included is (hopefully) self-explanatory\*. After that, just run the code below and Snakemake will take care of the rest. Simple as that.
```
snakemake [-j {threads}] [--resources io_limit={concurrency}] --configfile config.yaml
```
* Flags in [brackets] are optional.
* `{threads}` is an integer value for the number of cores. 
* `{concurrency}` is also an integer, and determines the number of maximum I/O intensive jobs to run simultaneously (we recommend starting with 1-3 and keeping an eye on the `%iowait` column of `sar` to see how much your machine/cluster can handle). Alternatively, the high I/O jobs can be pushed to RAM by changing `use_dev_shm: True` in the config file, but we only advise that if you are running on a big cluster with up to 15Gb spare RAM and `/dev/shm` is writable.

Because BMO is under active development, we strongly encourage you to run it this way, as we will modify the pipeline steps as the code matures. Running each step manually can potentially require extensive code rewriting.

\* But feel free to e-mail me at albanus@umich.edu in case it's not.

## Manually (standalone version)
*Coming soon*

# Interpreting results
Once the pipeline finishes running, you can find results in the `work/bmo/{sample}/bound` folder. Each BED file there will correspond to the predicted bound motif instances for each PWM scan. The format is a standard BED6 format, with the -log<sub>10</sub> adjusted *p*-value in the 5<sup>th</sup> column. If you want to use different filtering than our defaults, then look for the files at `work/bmo/{sample}/all`. The output is as described below:
1. Chromosome name
2. Start position
3. End position
4. Feature name
5. Feature score (*e.g.* FIMO score)
6. Strand
7. Number of ATAC-seq fragments in feature vicinity
8. Number of co-occurring motifs in feature vicinity
9. ATAC-seq fragments negative bionomial *p*-value
10. Co-occurring negative bionomial motifs *p*-value
11. Combined *p*-value
12. Adjusted combined *p*-value (using Benjamini-Yekutieli)
13. BMO score (-log<sub>10</sub> adjusted *p*-value)
14. Predicted bound 0/1 (1 is *yes*)
