#! /usr/bin/env python
#
# The Parker Lab (theparkerlab.org)
# University of Michigan, Ann Arbor
#

from os.path import join

###
# HELPER FUNCTIONS
#


def get_motifs():
    """List all motifs"""
    for motif in open(config['motif_file'], 'r').readlines():
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
# PIPELINE
#

# Desired output
#

rule all:
    input:
        expand(
            join(config['results'], "bmo", "{sample}",
                         "bound", "{motif}.bound.bed"),
            sample=get_samples(), motif=get_motifs()
        ),
        # expand(
        #     join(config['results'], "fvice", "{sample}", 
        #          "vplots", "{motif}.png"),
        #     sample=get_samples(), motif=get_motifs()
        # ),


if config["use_shm"] is True:
    print("Motifs will be processed using /dev/shm...")
    include:
        join(config["bmo_dir"], "rules",
                     "measure_raw_signal_shm.smk")
else:
    print("Motifs will not be processed using /dev/shm...")
    include:
        join(config["bmo_dir"], "rules", "measure_raw_signal.smk")


rule motifs_outside_peaks:
    input:
        motif = rules.measure_raw_signal.output,
        peaks = lambda wildcards: get_peaks(wildcards.sample)
    output:
        join(config['results'], "raw_signals",
                     "{sample}", "outside_peaks", "{motif}.bed")
    conda:
        "envs/bmo.yml"
    shell:
        "{IONICE} intersectBed -v -a {input.motif} -b {input.peaks} > {output}"

rule count_overlapping_motifs:
    input:
        rules.measure_raw_signal.output
    output:
        join(config['results'], "raw_signals",
                     "{sample}", "co_occurring", "{motif}.bed")
    params:
        com_src = join(
            config["bmo_dir"], "src", "count_co-occuring_motifs.sh"),
        bmo_dir = config["bmo_dir"],
        vicinity = 100,
    conda:
        "envs/bmo.yml"
    shell:
        """
        {IONICE} {params.com_src} {input} {params.vicinity} {params.bmo_dir}\
            > {output}
        """

rule fit_nbinoms:
    input:
        raw_signal = rules.measure_raw_signal.output,
        motif_counts = rules.count_overlapping_motifs.output,
        outside_peaks = rules.motifs_outside_peaks.output
    output:
        join(config['results'], "negative_binomials",
                     "{sample}", "{motif}.bed.gz")
    params:
        nb_src = join(config["bmo_dir"], "src", "rawSigNBModel.R"),
        in_handle = "{motif}.bed",
        out_handle = config['results'] + \
            "/negative_binomials/{sample}/{motif}",
        # note trailing slashes on both
        d1 = join(config['results'],
                          "raw_signals", "{sample}", "all_regions/"),
        d2 = join(config['results'], "raw_signals",
                          "{sample}", "outside_peaks/"),
    conda:
        "envs/bmo.yml"
    shell:
        """
        {IONICE} Rscript {params.nb_src} -f {params.in_handle} \
            --dir1 {params.d1} --dir2 {params.d2} \
            -o {params.out_handle} -c 7 --writeBed
        """

rule BMO:
    input:
        atac_nb = rules.fit_nbinoms.output,
        motif_counts = rules.count_overlapping_motifs.output
    output:
        join(config['results'], "bmo",
                     "{sample}", "bound", "{motif}.bound.bed")
    params:
        bmo_src = join(config["bmo_dir"], "src", "bmo.R"),
        bmo_output_dir = join(config['results'], "bmo", "{sample}")
    conda:
        "envs/bmo.yml"
    threads:
        1
    shell:
        """
        {IONICE} Rscript {params.bmo_src} --f1 {input.atac_nb} \
            --f2 {input.motif_counts} \
            -o {params.bmo_output_dir} -n {wildcards.motif} -p {threads}
        """

# rule filter_co_occurring:
#     input:
#         rules.BMO.output,
#     output:
#         join(config['results'], "fvice", "{sample}", "bound_no_overlap", "{motif}.bed"),
#     params:
#         py = join(config["bmo_dir"], "src", "filter_bed_co-occurring.py"),
#     conda:
#         "envs/bmo.yml"
#     shell:
#         """
#         {IONICE} {params.py} -i {input} -d 500 -c 5 -t strict > {output}
#         """

# rule vsignal:
#     input:
#         bam = lambda wildcards: config["samples"][wildcards.sample]["bamfile"],
#         bai = lambda wildcards: config["samples"][wildcards.sample]["bamfile"] + ".bai",
#         motif = rules.filter_co_occurring.output
#     output:
#         join(config['results'], "fvice", "{sample}", "vplots", "{motif}.vsignal.gz")
#     params:
#         py = join(config["bmo_dir"], "src", "measure_signal"),
#         flags = "-r 500 -f 3 -F 4 -F 8 -q 30",
#     conda:
#         "envs/bmo.yml"
#     resources:
#         io_limit = 1
#     threads: 4
#     shell:
#         """
#         {IONICE} {params.py} -p {threads} {params.flags} {input.bam} \
#             {input.motif} | gzip -c > {output}
#         """

# rule vplot:
#     input:
#         rules.vsignal.output
#     output:
#         join(config['results'], "fvice", "{sample}", "vplots", "{motif}.png"),
#     params:
#         R = join(config["bmo_dir"], "src", "makeVplots.R"),
#         handle = join(config['results'], "fvice", "{sample}", "vplots", "{motif}"),
#         name = "{sample}" + "_{motif}"
#     conda:
#         "envs/bmo.yml"
#     shell:
#         """
#         {IONICE} Rscript {params.R} -f {input} -o {params.handle} 
#             -n {params.name} --ylim 1.0
#         """