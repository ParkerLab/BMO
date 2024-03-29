# Regular version
rule measure_raw_signal:
    input:
        motif = lambda wildcards: os.path.join(config['motif_dir'], "{}.{}".format(wildcards.motif, config["motif_ext"])),
        bam = lambda wildcards: config["samples"][wildcards.sample]["bamfile"],
        bai = lambda wildcards: config["samples"][wildcards.sample]["bamfile"] + ".bai",
    output:
        os.path.join(config['results'], "raw_signals", "{sample}", "all_regions", "{motif}.bed")
    params:
        mrs_src = os.path.join(config["bmo_dir"], "src", "measureRawSignal.py"),
    threads: 15
    resources:
        io_limit = 1
    shell:
        "{IONICE} {params.mrs_src} -p {threads} -b {input.bam} -m {input.motif} > {output}"
