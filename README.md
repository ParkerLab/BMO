# ParkerLab BMO

# Overview
BMO (pronounced *beemo* - yes, that BMO) is an algorithm to predict TF binding from ATAC-seq data without using footprints. BMO uses negative binomial models of ATAC-seq fragments and number of co-occurring motifs to determine the likelihood of a given motif instance being bound.

**Disclaimer**: BMO has only been tested with paired-end ATAC-seq data.
**Platforms tested**: Linux, macOS

# Required input
1. ATAC-seq experiment pruned <sup>1</sup> and indexed BAM file
2. ATAC-seq peak calls <sup>2</sup>
3. Motif BED6 <sup>3</sup> files. One file per PWM scan <sup>4</sup>

<sup>1</sup> High quality, no duplicates, and properly mapped reads only.
<sup>2</sup> We call them with [MACS2](https://github.com/taoliu/MACS), using flags `--broad --nomodel --shift -100 --extsize 200 --keep-dup all` and then intersect out [blacklisted regions](https://sites.google.com/site/anshulkundaje/projects/blacklists).
<sup>3</sup> Mandatory. Pipeline will crash otherwise.
<sup>4</sup> We use [FIMO](http://meme-suite.org/doc/fimo.html).

# Running BMO
## Using Snakemake
We strongly recommend you run BMO using [Snakemake](https://snakemake.readthedocs.io), an extremely powerful tool for running reproducible pipelines. In order to run BMO this way, you need to copy both the `Snakefile` and `input/config.yaml` to wherever you want to run the analysis in your machine/cluster and then update the config file as described in the next section. After that, just run the code below and Snakemake will take care of the rest. Simple as that.
```
snakemake [-j {threads}] [--resources io_limit={concurrency}] --configfile config.yaml
```
* Flags in [brackets] are optional.
* `{threads}` is an integer value for the total number of cores Snakemake is allowed to use. if running on a HPC cluster, then set to 999 or some arbitrarily high value.
* `{concurrency}` is also an integer, and determines the number of maximum I/O intensive jobs to run simultaneously (we recommend starting with 1-3 and keeping an eye on the `%iowait` column of `sar` to see how much your machine can handle). Alternatively, the high I/O jobs can be pushed to RAM or to an SSD partition by changing `use_shm: True` in the config file (see details in the next section). If using the latter option, then also add `shm_limit={shm_concurrency}` to the `--resources` call above, where `{shm_concurrency}` is an integer for the maximum number of concurring jobs when using the RAM/SSD partition. If either `io_limit` or `shm_limit` are not set, then all jobs will be submitted with no regards to maximum concurrency (which should not be an issue if running in a cluster).

Reminder: Because BMO is still under active development, we again strongly encourage you to run it using Snakemake. Running each step manually may potentially require extensive code rewriting as the code matures.


### Setting up the configuration file
You can find the configuration file skeleton in `config/config.yaml`. Most fields should be self-explanatory, but let's go through all of them just in case:
* `bmo_dir`: this is the folder where you cloned this repository.
* `motif_dir`: path to your motif BED6 files.
* `motif_ext`: extension associated with the motif files. Do not add a dot in the beginning (`bed.gz` instead of `.bed.gz`).
* `motif_file`: plain text file containing the list of motif names (no paths or extensions), one per line.
* `results`: directory where the output will be saved.
* `ionice`: a string to append before each command call to set disk/CPU priority (*e.g.* `ionice -c2 -n7`). For Mac users: `ionice` is not supported in macOS, so the default value should be set to blank or `""`.
* `use_shm`: *True* or *False*. For some of the I/O heavy steps, BMO can use a shared memory partition (*i.e.* RAM that behave as disk space, like `/dev/shm`). This can increase BMO's concurrency when submitting jobs on a single machine with >10Gb of RAM to spare. Can also be used to send jobs to a SSD partition, instead.
* `shm_dir`: path for shared memory operations if `use_shm: True`.
* `samples`: information about each ATAC-seq experiment to be processed. The `bamfile` and `peakfile` fields under each sample handle are mandatory and point to the experiment's BAM and peak calls, respectively. Attention: make sure that the sample handles are unique.

IMPORTANT: Because the config is a [yaml](http://yaml.org) file, **make sure to maintain all the indentations as they appear** or unexpected behaviors may occurr.


## Manually (standalone version)
*Coming soon*



# Interpreting BMO results
## Processed BED files
Once the pipeline finishes running, you can find results in the `work/bmo/{sample}/bound` folder. The BED6 files there will correspond to the **predicted bound motif instances** for each motif.

1. Chromosome name
2. Start position
3. End position
4. Feature name
5. BMO score (-log<sub>10</sub> adjusted *p*-value)
6. Strand


## Full output
If you want to use different filtering thresholds than our defaults, then look for the files at `work/bmo/{sample}/all`. The output is as described below:

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