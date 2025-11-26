#!/usr/bin/env python3
"""
Local test wrapper to simulate SLURM array execution for RNA-seq pipeline.
It loops over all _R1 fastq files in the input directory and runs
run_pipeline.sh for each sample, setting SLURM_ARRAY_TASK_ID.
"""

import subprocess
import os
from pathlib import Path
import argparse
import shlex

def parse_args():

    parser = argparse.ArgumentParser(description="Local test of SLURM array pipeline")

    parser.add_argument(
        "--root",
        required=True,
        help="Path to pipeline root"
        )
    
    parser.add_argument(
        "--indir",
        required=True,
        help="Directory containing input fastq files"
        )
    
    parser.add_argument(
        "--steps",
        nargs="+",
        default=["trim", "align", "count"],
        help="Pipeline steps to run"
    )

    parser.add_argument(
        "--conda-env",
        default="rnaseq",
        help="Conda environment to use"
    )

    parser.add_argument(
        "--runScript",
        default="run_pipeline.sh",
        help="Path to the SLURM wrapper script"
    )

    return parser.parse_args()

def main():
    args = parse_args()
    indir = Path(args.indir).resolve()
    run_script = Path(args.root).resolve() / "scripts" / args.runScript

    r1_files = sorted(indir.glob("*_R1*.fastq*"))
    if not r1_files:
        raise FileNotFoundError(f"No _R1 fastq files found in {indir}")

    print(f"Found {len(r1_files)} samples in {indir}. Running locally...\n")

    for idx, r1 in enumerate(r1_files):
        r2 = Path(str(r1).replace("_R1", "_R2"))
        if not r2.exists():
            print(f"Skipping {r1.name}: missing paired file {r2.name}")
            continue

        # set environment variable to simulate SLURM task ID
        os.environ["SLURM_ARRAY_TASK_ID"] = str(idx)

        # build command
        cmd = [
            "bash",
            str(run_script),
            "--root", str(args.root),
            "--indir", str(indir),
            "--steps", *args.steps,
            "--conda-env", args.conda_env
        ]

        print(f"\n=== Running sample {idx}: {r1.name} ===")
        print(f"Command: {' '.join(shlex.quote(c) for c in cmd)}\n")

        # run the pipeline
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)

if __name__ == "__main__":
    main()
