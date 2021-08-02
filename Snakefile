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


# Config
IONICE = config["ionice"]

BMO_SRC = config["bmo_dir"]
RESULTS = config['results']
BAM_DIR = config["bamdir"]
PEAK_DIR = config["summits"]
MOTIF_DIR = config["motif_dir"]

# Wildcards
motifs = get_motifs()
samples = config["clusters"]



rule all:
    input:
        expand(
            join(RESULTS, "fvice", "{sample}", "vplots", "{motif}.png"),
            sample=samples, motif=motifs
        ),



rule measure_raw_signal:
    input:
        motif = join(MOTIF_DIR, "{motif}.bed.gz"),
        bam = join(BAM_DIR, "atac__{sample}.bam"),
        bai = join(BAM_DIR, "atac__{sample}.bam.bai"),
    output:
        join(RESULTS, "raw_signals", "{sample}", "all_regions", "{motif}.bed")
    params:
        mrs_src = join(BMO_SRC, "src", "measureRawSignal.py"),
    conda:
        "envs/bmo.yml"
    threads: 5
    resources:
        io_limit = 1
    shell:
        """
        {IONICE} {params.mrs_src} -p {threads} -b {input.bam} \
            -m {input.motif} > {output}
        """

rule motifs_outside_peaks:
    input:
        motif = rules.measure_raw_signal.output,
        peaks = join(PEAK_DIR, "{sample}_peaks.min2.bed"),
    output:
        join(RESULTS, "raw_signals", "{sample}", "outside_peaks", "{motif}.bed")
    conda:
        "envs/bmo.yml"
    shell:
        "{IONICE} intersectBed -v -a {input.motif} -b {input.peaks} > {output}"

rule count_overlapping_motifs:
    input:
        rules.measure_raw_signal.output
    output:
        join(RESULTS, "raw_signals", "{sample}", "co_occurring", "{motif}.bed")
    params:
        com_src = join(
            BMO_SRC, "src", "count_co-occuring_motifs.sh"),
        bmo_dir = BMO_SRC,
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
        join(RESULTS, "negative_binomials",
                     "{sample}", "{motif}.bed.gz")
    params:
        nb_src = join(BMO_SRC, "src", "rawSigNBModel.R"),
        in_handle = "{motif}.bed",
        out_handle = RESULTS + \
            "/negative_binomials/{sample}/{motif}",
        # note trailing slashes on both
        d1 = join(RESULTS,
                          "raw_signals", "{sample}", "all_regions/"),
        d2 = join(RESULTS, "raw_signals",
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
        join(RESULTS, "bmo",
                     "{sample}", "bound", "{motif}.bound.bed")
    params:
        bmo_src = join(BMO_SRC, "src", "bmo.R"),
        bmo_output_dir = join(RESULTS, "bmo", "{sample}")
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

rule filter_co_occurring:
    input:
        rules.BMO.output,
    output:
        join(RESULTS, "fvice", "{sample}", "bound_no_overlap", "{motif}.bed"),
    conda:
        "envs/bmo.yml"
    shell:
        """
        {IONICE} ~/scripts/filter_bed_co-occurring.py -i {input} \
            -d 500 -c 5 -t strict > {output}
        """

rule vsignal:
    input:
        bam = rules.measure_raw_signal.input.bam,
        motif = rules.filter_co_occurring.output
    output:
        join(RESULTS, "fvice", "{sample}", "vplots", "{motif}.vsignal.gz")
    params:
        py = join(BMO_SRC, "src", "measure_signal"),
        flags = "-r 500 -f 3 -F 4 -F 8 -q 30",
    conda:
        "envs/bmo.yml"
    resources:
        io_limit = 1
    threads: 4
    shell:
        """
        {IONICE} {params.py} -p {threads} {params.flags} {input.bam} \
            {input.motif} | gzip -c > {output}
        """

rule vplot:
    input:
        rules.vsignal.output
    output:
        join(RESULTS, "fvice", "{sample}", "vplots", "{motif}.png"),
    params:
        R = join(BMO_SRC, "src", "makeVplots.R"),
        handle = join(RESULTS, "fvice", "{sample}", "vplots", "{motif}"),
        name = "{sample}" + "_{motif}"
    conda:
        "envs/bmo.yml"
    shell:
        """
        {IONICE} Rscript {params.R} -f {input} -o {params.handle} 
            -n {params.name} --ylim 1.0
        """