# /dev/shm version

rule push_bam_to_shm:
    input:
        bam = lambda wildcards: config["samples"][wildcards.sample]["bamfile"],
        bai = lambda wildcards: config["samples"][
            wildcards.sample]["bamfile"] + ".bai"
    output:
        bam = temporary(
            join(config["shm_dir"], "{sample}.pruned.bam")
        ),
        bai = temporary(
            join(config["shm_dir"], "{sample}.pruned.bam.bai")
        ),
    conda:
        "../envs/bmo.yml"
    shell:
        """
        {IONICE} cp {input.bam} {output.bam} && samtools index {output.bam}
        """

rule push_motif_to_shm:
    input:
        lambda wildcards: join(
            config['motif_dir'], "{}.{}".format(
                wildcards.motif, config["motif_ext"])
        )
    output:
        temporary(join(config["shm_dir"], "motifs", "{motif}.bed.gz"))
    shell:
        """
        {IONICE} cp {input} {output}
        """

rule measure_raw_signal:
    input:
        motif = rules.push_motif_to_shm.output,
        bam = rules.push_bam_to_shm.output.bam,
        bai = rules.push_bam_to_shm.output.bai,
    output:
        join(config['results'], "raw_signals", "{sample}",
                     "all_regions", "{motif}.bed")
    params:
        mrs_script = join(
            config["bmo_dir"], "src", "measureRawSignal.py"),
        tmp = join(config["shm_dir"],
                           "{sample}_{motif}.raw_signal.bed"),
    conda:
        "../envs/bmo.yml"
    threads:
        15
    resources:
        shm_limit = 1
    shell:
        """
        {params.mrs_script} -p {threads} -b {input.bam} -m {input.motif} > \
            {params.tmp} && mv {params.tmp} {output}
        """
