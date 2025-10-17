ISOLATE_FP = Cfg["all"]["output_fp"] / "isolate"
TOOLS = {
    "shovill": ["number of contigs", "min coverage", "max coverage", "mean coverage"],
    "checkm": ["Completeness", "Contamination"],
    "sylph": ["Taxonomic_abundance", "Contig_name"],
    "mlst": ["Schema", "ST", "Alleles"],
    "bakta": [
        "Length",
        "GC",
        "N50",
        "CDSs",
        "tRNAs",
        "tmRNAs",
        "rRNAs",
        "hypotheticals",
        "CRISPR arrays",
    ],
    "mash": ["Mash_Contamination", "Contaminated_Spp"],
}

try:
    SBX_SGA_VERSION = get_ext_version("sbx_sga")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_SGA_VERSION = "0.0.0"


localrules:
    all_sga,


rule all_sga:
    input:
        expand(ISOLATE_FP / "quast" / "{sample}" / "report.tsv", sample=Samples),
        f"{ISOLATE_FP}/final_summary.tsv",
        f"{ISOLATE_FP}/reports/amr.report",


## Assembly
rule sga_shovill:
    input:
        rp1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        rp2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    log:
        LOG_FP / "sga_shovill_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_shovill_{sample}.tsv"
    conda:
        "envs/shovill.yml"
    resources:
        mem_mb=32000,
    shell:
        """
        (shovill --force --assembler skesa --outdir $(dirname {output.contigs}) --R1 {input.rp1}  --R2 {input.rp2} &> {log} &&
        mv $(dirname {output.contigs})/contigs.fa {output.contigs}) ||
        touch {output.contigs}
        """


## Parse shovill output
rule shovill_summary:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        statistics=ISOLATE_FP / "shovill" / "{sample}" / "parsed_summary.tsv",
    log:
        LOG_FP / "sga_shovill_summary_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_shovill_summary_{sample}.tsv"
    script:
        "scripts/shovill.py"


# Combine shovill outputs
rule combine_shovill_summary:
    input:
        reports=expand(
            ISOLATE_FP / "shovill" / "{sample}" / "parsed_summary.tsv", sample=Samples
        ),
    output:
        shovill_report=ISOLATE_FP / "reports" / "shovill.report",
    params:
        suffix="",
        header=True,
    log:
        LOG_FP / "sga_combine_shovill_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_combine_shovill_summary.tsv"
    script:
        "scripts/concat_files.py"


# Taxonomic classification
rule sga_sylph:
    input:
        rp1=QC_FP / "decontam" / "{sample}_1.fastq.gz",
        rp2=QC_FP / "decontam" / "{sample}_2.fastq.gz",
    output:
        report=ISOLATE_FP / "sylph" / "{sample}" / "{sample}.tsv",
    threads: 8
    params:
        ref=Cfg["sbx_sga"]["sylph_ref"],
    log:
        LOG_FP / "sga_sylph_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_sylph_{sample}.tsv"
    conda:
        "envs/sylph.yml"
    resources:
        mem_mb=32000,
        runtime=120,
    shell:
        """
    if [ $(zcat {input.rp1} | wc -l) -ge 4 ] && [ $(zcat {input.rp2} | wc -l) -ge 4 ]; then
        sylph sketch -1 "{input.rp1}" -2 "{input.rp2}" -t 1 -d "$(dirname "{output.report}")" > "{log}" 2>&1
        sylph profile "{params.ref}" "$(dirname "{output.report}")"/*.sylsp -t {threads} -o "{output.report}" >> "{log}" 2>&1
    else
        touch {output.report}
    fi
        """


rule combine_sylph_summary:
    input:
        summaries=expand(
            ISOLATE_FP / "sylph" / "{sample}" / "{sample}.tsv", sample=Samples
        ),
    output:
        all_summary=ISOLATE_FP / "reports" / "sylph.report",
    log:
        LOG_FP / "sga_combine_sylph_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_combine_sylph_summary.tsv"
    params:
        suffix="",
        header=True,
    script:
        "scripts/concat_files.py"


### Assembly QC
rule sga_checkm:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        quality_report=ISOLATE_FP / "checkm" / "{sample}" / "quality_report.tsv",
    params:
        ref=Cfg["sbx_sga"]["checkm_ref"],
    log:
        LOG_FP / "sga_checkm_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_checkm_{sample}.tsv"
    conda:
        "envs/checkm2.yml"
    resources:
        mem_mb=16000,
    shell:
        """
        if [ -s {input.contigs} ]; then
            checkm2 predict \\
            --force \\
            -x fa \\
            -i {input.contigs} \\
            -o $(dirname {output.quality_report}) \\
            --database_path {params.ref} \\
            &> {log} || touch {output.quality_report}
        else
            touch {output.quality_report}
        fi
        """


rule combine_checkm_summary:
    input:
        summaries=expand(
            ISOLATE_FP / "checkm" / "{sample}" / "quality_report.tsv", sample=Samples
        ),
    output:
        all_summary=ISOLATE_FP / "reports" / "checkm.report",
    log:
        LOG_FP / "sga_combine_checkm_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_combine_checkm_summary.tsv"
    params:
        suffix="",
        header=True,
    script:
        "scripts/concat_files.py"


rule sga_quast:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        quast_dir=ISOLATE_FP / "quast" / "{sample}" / "report.tsv",
    log:
        LOG_FP / "sga_quast_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_quast_{sample}.tsv"
    conda:
        "envs/quast.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            quast.py \\
                -o $(dirname {output.quast_dir}) \\
                {input.contigs} \\
                &> {log} || touch {output.quast_dir}
        else
            touch {output.quast_dir}
        fi
        """


rule sga_mash:
    input:
        reads=expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        agg=temp(ISOLATE_FP / "mash" / "{sample}" / "{sample}.fastq"),
        win=temp(ISOLATE_FP / "mash" / "{sample}" / "{sample}_winning.tab"),
        sort=ISOLATE_FP / "mash" / "{sample}" / "{sample}_sorted_winning.tab",
    params:
        ref=Cfg["sbx_sga"]["mash_ref"],
    log:
        LOG_FP / "sga_mash_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mash_{sample}.tsv"
    conda:
        "envs/mash.yml"
    resources:
        mem_mb=16000,
    shell:
        """
        zcat {input.reads} > {output.agg}

        if [ -s {output.agg} ]; then
            mash screen -w -p 8 {params.ref} {output.agg} > {output.win} 2> {log}
            sort -gr {output.win} > {output.sort} 2>> {log}
        else
            touch {output.win} {output.sort}
        fi
        """


rule mash_summary:
    input:
        sorted_reports=ISOLATE_FP / "mash" / "{sample}" / "{sample}_sorted_winning.tab",
    output:
        summary=ISOLATE_FP / "mash" / "{sample}" / "{sample}_summary.tsv",
    log:
        LOG_FP / "sga_mash_summary_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mash_summary_{sample}.tsv"
    script:
        "scripts/mash.py"


rule combine_mash_summary:
    input:
        summary_reports=expand(
            ISOLATE_FP / "mash" / "{sample}" / "{sample}_summary.tsv", sample=Samples
        ),
    output:
        mash_report=ISOLATE_FP / "reports" / "mash.report",
    shell:
        """
        echo -e "SampleID\\tMash_Contamination\\tContaminated_Spp" > {output.mash_report}
        cat {input.summary_reports} >> {output.mash_report}
        """


# Typing
rule sga_mlst:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        mlst=ISOLATE_FP / "mlst" / "{sample}" / "{sample}.mlst",
    log:
        LOG_FP / "sga_mlst_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mlst_{sample}.tsv"
    conda:
        "envs/mlst.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            mlst --nopath {input.contigs} > {output.mlst} 2> {log}
        else
            mkdir -p $(dirname {output.mlst})
            touch {output.mlst}
        fi
        """


rule mlst_parse:
    input:
        reports=ISOLATE_FP / "mlst" / "{sample}" / "{sample}.mlst",
    output:
        mlst_report=ISOLATE_FP / "mlst" / "{sample}" / "parsed_mlst.txt",
    log:
        LOG_FP / "sga_mlst_parse_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_mlst_parse_{sample}.tsv"
    script:
        "scripts/mlst.py"


rule mlst_summary:
    input:
        reports=expand(
            ISOLATE_FP / "mlst" / "{sample}" / "parsed_mlst.txt", sample=Samples
        ),
    output:
        amr_report=ISOLATE_FP / "reports" / "mlst.report",
    params:
        suffix="",
        header=True,
    log:
        LOG_FP / "sga_mlst_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_mlst_summary.tsv"
    script:
        "scripts/concat_files.py"


### Annotation
rule sga_bakta:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        bakta=ISOLATE_FP / "bakta" / "{sample}" / "{sample}.txt",
    params:
        ref=Cfg["sbx_sga"]["bakta_ref"],
    log:
        LOG_FP / "sga_bakta_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_bakta_{sample}.tsv"
    conda:
        "envs/bakta.yml"
    resources:
        mem_mb=32000,
        runtime=180,
    shell:
        """
        if [ -s {input.contigs} ]; then
            bakta --force --db {params.ref} \\
            --output $(dirname {output.bakta}) \\
            --prefix {wildcards.sample} \\
            --skip-plot {input.contigs} \\
            &> {log}
        else
            touch {output.bakta}
        fi
        """


rule parse_bakta_report:
    input:
        bakta=ISOLATE_FP / "bakta" / "{sample}" / "{sample}.txt",
    output:
        parsed_report=ISOLATE_FP / "bakta" / "{sample}" / "parsed_summary.tsv",
    log:
        LOG_FP / "sga_bakta_parse_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_bakta_parse_{sample}.tsv"
    script:
        "scripts/bakta.py"


rule combine_bakta_summary:
    input:
        reports=expand(
            ISOLATE_FP / "bakta" / "{sample}" / "parsed_summary.tsv", sample=Samples
        ),
    output:
        bakta_report=ISOLATE_FP / "reports" / "bakta.report",
    params:
        suffix="",
        header=True,
    log:
        LOG_FP / "sga_combine_bakta_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_combine_bakta_summary.tsv"
    script:
        "scripts/concat_files.py"


### AMR Profiling
rule sga_abritamr:
    input:
        contigs=ISOLATE_FP / "shovill" / "{sample}" / "{sample}.fa",
    output:
        abritamr=ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out",
    log:
        LOG_FP / "sga_abritamr_{sample}.log",
    benchmark:
        BENCHMARK_FP / "sga_abritamr_{sample}.tsv"
    conda:
        "envs/abritamr.yml"
    shell:
        """
        if [ -s {input.contigs} ]; then
            abritamr run \\
                --contigs {input.contigs} \\
                --prefix {wildcards.sample} \\
                &> {log}
            mv {wildcards.sample} $(dirname $(dirname {output.abritamr}))
        else
            mkdir -p $(dirname {output.abritamr})
            touch {output.abritamr}
        fi     
    """


rule abritamr_summary:
    input:
        reports=expand(
            ISOLATE_FP / "abritamr" / "{sample}" / "amrfinder.out", sample=Samples
        ),
    output:
        amr_report=ISOLATE_FP / "reports" / "amr.report",
    params:
        suffix="",
        header=True,
    log:
        LOG_FP / "sga_abritamr_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_abritamr_summary.tsv"
    script:
        "scripts/concat_files.py"


# Final Summary Report
rule all_summary:
    input:
        reports=expand(ISOLATE_FP / "reports" / "{tool}.report", tool=TOOLS.keys()),
    output:
        final_report=ISOLATE_FP / "final_summary.tsv",
    params:
        tools=TOOLS,
    log:
        LOG_FP / "sga_all_summary.log",
    benchmark:
        BENCHMARK_FP / "sga_all_summary.tsv"
    script:
        "scripts/summarize_all.py"


rule clean_sga:
    input:
        QC_FP / ".decontam_cleaned",
        rules.all_sga.input,
    output:
        touch(ISOLATE_FP / ".sga_cleaned"),
    params:
        log_dir=LOG_FP,
        benchmark_dir=BENCHMARK_FP,
    shell:
        """
        isolate_dir=$(dirname {output[0]})
        log_dir={params.log_dir}
        benchmark_dir={params.benchmark_dir}

        rm -rf $log_dir
        rm -rf $benchmark_dir
        find "$isolate_dir/mash" -type f -name "*.fastq" -delete
        """


rule test_sga:
    input:
        f"{ISOLATE_FP}/reports/shovill.report",
        f"{ISOLATE_FP}/reports/mlst.report",
        f"{ISOLATE_FP}/reports/checkm.report",
        f"{ISOLATE_FP}/reports/amr.report",
        #f"{ISOLATE_FP}/reports/bakta.report", # Missing test database
        f"{ISOLATE_FP}/reports/mash.report",
