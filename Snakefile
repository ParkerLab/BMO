#! /usr/bin/env python
#
# The Parker Lab (theparkerlab.org)
# University of Michigan, Ann Arbor
#

###
## HELPER FUNCTIONS
#

def get_motifs():
    """List all motifs"""
    for motif in open(config['motif_file'],'r').readlines():
        yield motif.strip()

def get_samples():
    """List all samples"""
    for sample in sorted(config['samples'].keys()):
        yield sample

def get_peaks(sample):
    """Get peak file associated with sample from the config"""
    peakfile = config["samples"][sample]["peakfile"]
    return peakfile

if config["ionice"] is None:
    IONICE = ""
else:
    IONICE = config["ionice"]

###
## PIPELINE
#

## Desired output
#

rule all:
    input:
        expand(os.path.join(config['results'], "bmo", "{sample}", "bound", "{motif}.bound.bed"), sample=get_samples(), motif=get_motifs())


if config["use_shm"] is True:
    print("Motifs will be processed using /dev/shm...")
    include: os.path.join(config["bmo_dir"], "rules", "measure_raw_signal_shm.smk")
else:
    print("Motifs will not be processed using /dev/shm...")
    include: os.path.join(config["bmo_dir"], "rules", "measure_raw_signal.smk")


rule motifs_outside_peaks:
    input:
        motif = rules.measure_raw_signal.output,
        peaks = lambda wildcards: get_peaks(wildcards.sample)
    output:
       os.path.join(config['results'], "raw_signals", "{sample}", "outside_peaks", "{motif}.bed")
    shell:
        "{IONICE} intersectBed -v -a {input.motif} -b {input.peaks} > {output}"

rule count_overlapping_motifs:
    input:
        rules.measure_raw_signal.output
    output:
          os.path.join(config['results'], "raw_signals", "{sample}", "co_occurring", "{motif}.bed")
    params: 
        com_src = os.path.join(config["bmo_dir"], "src", "count_co-occuring_motifs.sh"),
        bmo_dir = config["bmo_dir"],
        vicinity = 100
    shell:
        "{IONICE} {params.com_src} {input} {params.vicinity} {params.bmo_dir} > {output}"

rule fit_nbinoms:
    input:
        raw_signal   = rules.measure_raw_signal.output,
        motif_counts = rules.count_overlapping_motifs.output,
        outside_peaks = rules.motifs_outside_peaks.output
    output:
        os.path.join(config['results'], "negative_binomials", "{sample}", "{motif}.bed.gz")
    params:
        nb_src = os.path.join(config["bmo_dir"], "src", "rawSigNBModel.R"),
        in_handle = "{motif}.bed",
        out_handle = config['results'] + "/negative_binomials/{sample}/{motif}",
        d1 = os.path.join(config['results'], "raw_signals", "{sample}", "all_regions/"), # note trailing slashes on both
        d2 = os.path.join(config['results'], "raw_signals", "{sample}", "outside_peaks/"),
    shell:
        """
        {IONICE} Rscript {params.nb_src} -f {params.in_handle} --dir1 {params.d1} --dir2 {params.d2} \
        -o {params.out_handle} -c 7 --writeBed
        """

rule BMO:
    input:
        atac_nb = rules.fit_nbinoms.output,
        motif_counts = rules.count_overlapping_motifs.output
    output:
        os.path.join(config['results'], "bmo", "{sample}", "bound", "{motif}.bound.bed")
    params:
        bmo_src = os.path.join(config["bmo_dir"], "src", "bmo.R"),
        bmo_output_dir = os.path.join(config['results'], "bmo", "{sample}")
    threads: 1
    shell:
        """
        {IONICE} Rscript {params.bmo_src} --f1 {input.atac_nb} --f2 {input.motif_counts} \
        -o {params.bmo_output_dir} -n {wildcards.motif} -p {threads}
        """