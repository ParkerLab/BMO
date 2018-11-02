# /dev/shm version

rule push_bam_to_shm:
    input: 
        bam = lambda wildcards: config["samples"][wildcards.sample]["bamfile"],
        bai = lambda wildcards: config["samples"][wildcards.sample]["bamfile"] + ".bai"
    output:
        bam = temporary(os.path.join(config["shm_dir"], "{sample}.pruned.bam")),
        bai = temporary(os.path.join(config["shm_dir"], "{sample}.pruned.bam.bai"))
    shell:
        "{IONICE} cp {input.bam} {output.bam} && samtools index {output.bam}"

rule push_motif_to_shm:
    input:
        lambda wildcards: os.path.join(config['motif_dir'], "{}.{}".format(wildcards.motif, config["motif_ext"]))
    output:
        temporary(os.path.join(config["shm_dir"], "motifs", "{motif}.bed.gz"))
    shell:
        "{IONICE} cp {input} {output}"

rule measure_raw_signal:
    input:
        motif = rules.push_motif_to_shm.output,
        bam = rules.push_bam_to_shm.output.bam,
        bai = rules.push_bam_to_shm.output.bai,
    output:
        os.path.join(config['results'], "raw_signals", "{sample}", "all_regions", "{motif}.bed")
    params:
        mrs_src = os.path.join(config["bmo_dir"], "src", "measureRawSignal.py"),
    threads: 15
    resources:
        shm_limit = 1
    params:
        tmp = os.path.join(config["shm_dir"], "{sample}_{motif}.raw_signal.bed")
    shell:
        "{params.mrs_src} -p {threads} -b {input.bam} -m {input.motif} > {params} && mv {params} {output}"