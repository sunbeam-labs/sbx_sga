import subprocess as sp
import sys
from pathlib import Path

CPUS = 8

# Run script like: `python script_name.py /path/to/reads/ /path/to/reference_genome.fa`
# Will probably want to submit this to the cluster for real datasets
# Use an sbatch script and include the above as the main script body
if __name__ == "__main__":
    # Gather reads
    # For now, manually excluding mock/control samples
    reads_fp = Path(sys.argv[1])
    reads_files = list(reads_fp.glob("*.fastq.gz"))
    # Pair reads
    samples = {
        f.name.split("_1.fastq.gz")[0].split("_2.fastq.gz")[0] for f in reads_files
    }
    paired_reads = [
        (sample, reads_fp / f"{sample}_1.fastq.gz", reads_fp / f"{sample}_2.fastq.gz")
        for sample in samples
    ]

    # Select reference genome
    # Again, for now, doing this manually
    ref_fp = Path(sys.argv[2])

    # Run snippy per-sample
    for sample, r1, r2 in paired_reads:
        outdir = Path("snippy_outputs") / sample
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = [
            "snippy",
            "--cpus",
            str(CPUS),
            "--outdir",
            str(outdir),
            "--ref",
            ref_fp,
            "--R1",
            str(r1),
            "--R2",
            str(r2),
        ]
        print("Running ", " ".join(cmd))
        sp.run(cmd)

    # Run snippy-core to create core SNP alignment
    snippy_outdirs = [str(p) for p in Path("snippy_outputs").iterdir() if p.is_dir()]
    core_outdir = Path("snippy_core_output")
    core_outdir.mkdir(parents=True, exist_ok=True)
    cmd = ["snippy-core", "--prefix", "core", *snippy_outdirs, "--ref", ref_fp]
    print("Running ", " ".join(cmd))
    sp.run(cmd, cwd=str(core_outdir), shell=True)

    # Run iqtree on core SNP alignment
    core_snp_alignment = core_outdir / "core.full.aln"
    iqtree_outdir = core_outdir / "iqtree" / "core.full.aln.treefile"
    iqtree_outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "iqtree",
        "-s",
        str(core_snp_alignment),
        "-m",
        "MFP",
        "-B",
        "1000",
        "-T",
        "AUTO",
        "-pre",
        str(core_outdir),
    ]
    print("Running ", " ".join(cmd))
    sp.run(cmd)
