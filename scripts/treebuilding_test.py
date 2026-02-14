#!/usr/bin/env python3

from pathlib import Path
import subprocess as sp
import sys
import argparse
import shutil
from datetime import datetime

CPUS_DEFAULT = 8
SNIPPY_OUTDIR = Path("snippy_outputs").resolve()
CORE_OUTDIR = Path("snippy_core_output").resolve()

def run_cmd(cmd, cwd=None, capture_to=None):
    print(f"\n> Running: {cmd}")
    completed = sp.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)

    if capture_to:
        with open(capture_to, "w") as fh:
            fh.write(f"## CMD: {cmd}\n")
            fh.write(f"## Returncode: {completed.returncode}\n\n")
            fh.write("## STDOUT\n")
            fh.write(completed.stdout or "")
            fh.write("\n\n## STDERR\n")
            fh.write(completed.stderr or "")

    return completed.returncode

#### this feature will likely be deprecated later on with more advanced sample name handling
def detect_simple_pairs(reads_dir):
    """
    Assumes sample names are in this format for now:
        <sample>_R1.fastq.gz
        <sample>_R2.fastq.gz
    """

    r1_files = sorted(reads_dir.glob("*_R1.fastq.gz"))
    samples = {}

    for r1 in r1_files:
        sample = r1.name.replace("_R1.fastq.gz", "")
        r2 = reads_dir / f"{sample}_R2.fastq.gz"
        if r2.exists():
            samples[sample] = (r1, r2)
        else:
            print(f"WARNING: R2 missing for sample {sample}, skipping.")

    return samples
####

def main(argv):
    parser = argparse.ArgumentParser(description="treebuilding pipeline using snippy snippy-core and iqtree")
    parser.add_argument("reads_dir", type=Path)
    parser.add_argument("reference", type=Path)
    parser.add_argument("--cpus", type=int, default=CPUS_DEFAULT)
    parser.add_argument("-g", "--graft-tree", type=Path)
    args = parser.parse_args(argv)

    reads_dir = args.reads_dir
    ref = args.reference
    cpus = args.cpus
    print(args)
    if not reads_dir.is_dir():
        sys.exit(f"ERROR: Reads directory {reads_dir} does not exist")

    if not ref.exists():
        sys.exit(f"ERROR: Reference file {ref} does not exist")

    
    # Make snippy outdir
    SNIPPY_OUTDIR.mkdir(exist_ok=True)
    

    # Detect samples
    pairs = detect_simple_pairs(reads_dir)
    if not pairs:
        sys.exit("ERROR: No valid samples detected")

    print(f"Detected {len(pairs)} paired samples.")

    successful = []

    # Run snippy on new samples only, skip already processed samples
    for sample, (r1, r2) in pairs.items():
        outdir = SNIPPY_OUTDIR / sample
        outdir.mkdir(parents=True, exist_ok=True)

        runlog = outdir / f"{sample}.snippy.log"


        if (outdir / "snps.vcf").exists():
            print(f"Snippy already ran for: {sample}")
            successful.append(str(outdir))
        else:
            print(f"Running Snippy on {sample}...")
            cmd = (
		  f"snippy --cpus {cpus} "
		  f"--outdir {outdir} "
		  f"--ref {ref} "
		  f"--R1 {r1} "
		  f"--R2 {r2} "
		  f"--force"
		  )
            rc = run_cmd(cmd, capture_to=str(runlog))
            
            if (outdir / "snps.vcf").exists():
            	print(f"Snippy OK: {sample}")
            	successful.append(str(outdir))
            else:
            	print(f"Snippy failed: {sample}, see log: {runlog}")
            	
            

    if not successful:
        sys.exit("ERROR: No successful snippy outputs...stopping pipeline.")

    # Run snippy-core
    CORE_OUTDIR.mkdir(exist_ok=True)

    core_log = CORE_OUTDIR / "snippy_core.log"
    sample_dirs = " ".join(successful)

    cmd = f"snippy-core --prefix core {sample_dirs} --ref {ref}"
    rc = run_cmd(cmd, cwd=str(CORE_OUTDIR), capture_to=str(core_log))

    aln = CORE_OUTDIR / "core.full.aln"
    if not aln.exists() or aln.stat().st_size == 0:
        sys.exit("ERROR: snippy-core failed, core.full.aln not generated.")

    print("snippy-core OK.")

    # Run iqtree
    iqtree_dir = CORE_OUTDIR / "iqtree"
    iqtree_dir.mkdir(exist_ok=True)

    tree_prefix = iqtree_dir / "core.full.aln"
    
    if args.graft_tree is not None:
        cmd = (
        f"iqtree -s {aln} -m GTR+G -B 1000 -T AUTO -g {args.graft_tree} "
        f"-pre {tree_prefix}"
    ) 
    else:
    	cmd = (
        f"iqtree -s {aln} -m GTR+G -B 1000 -T AUTO "
        f"-pre {tree_prefix}"
    	)

    iqlog = iqtree_dir / "iqtree.log"

    rc = run_cmd(cmd, cwd=str(iqtree_dir), capture_to=str(iqlog))

    treefile = iqtree_dir / "core.full.aln.treefile"
    if not treefile.exists():
        sys.exit("ERROR: IQ-TREE failed, no treefile produced.")

    print("\nPIPELINE COMPLETED SUCCESSFULLY")
    print(f"snippy outputs: {SNIPPY_OUTDIR.resolve()}")
    print(f"snippy-core outputs: {CORE_OUTDIR.resolve()}")
    print(f"iqtree outputs: {treefile.resolve()}")


if __name__ == "__main__":
    main(sys.argv[1:])

