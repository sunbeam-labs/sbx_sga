ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
try:
    SBX_SGA_VERSION = get_ext_version("sbx_sga")
except NameError:
    SBX_SGA_VERSION = "0.0.0"


localrules:
    all_sga_snippy,


rule all_sga_snippy:
    input:
        ISOLATE_FP / "reports" / "snippy.report",


rule sga_snippy:
    input:
        rp1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        rp2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        vcf=ISOLATE_FP / "snippy" / "{sample}" / "snps.vcf",
        tab=ISOLATE_FP / "snippy" / "{sample}" / "snps.tab",
    params:
        ref=Cfg["sbx_sga"]["snippy_ref"],
    log:
        LOG_FP / "sga_snippy_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_snippy_{sample}.tsv"
    threads: 8
    conda:
        "envs/snippy.yml"
    shell:
        """
        snippy --cpus {threads} --outdir $(dirname {output.vcf}) --force \
            --ref {params.ref} --R1 {input.rp1} --R2 {input.rp2} > {log} 2>&1
        """


rule snippy_summary:
    input:
        expand(ISOLATE_FP / "snippy" / "{sample}" / "snps.tab", sample=Samples),
    output:
        ISOLATE_FP / "reports" / "snippy.report",
    script:
        "scripts/snippy.py"
