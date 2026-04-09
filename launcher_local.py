# region Imports

from pathlib import Path
<<<<<<< HEAD
import sys, argparse, subprocess, time, logging, os
=======
import sys, argparse, subprocess, time, logging
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.config_loader import ConfigLoader
<<<<<<< HEAD
from src.counts import Counts
=======
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6

# endregion

def parse_args():
    """
    Accepts CLI arguments for launcher_local.py
    """
    parser = argparse.ArgumentParser(description="Local launcher for RNA-seq pipeline, runs samples in parallel batches")

    parser.add_argument(
        "--root",
        required=True,
        help="Path to root directory of pipeline (RNAseq_BULK)"
    )

    parser.add_argument(
        "--indir",
        required=True,
        help="Path to directory containing input fastq files"
    )

    parser.add_argument(
        "--steps",
        nargs="+",
        default=["trim", "align", "count"],
        help="Pipeline steps to run (trim, align, count)"
    )

    parser.add_argument(
        "--max_parallel",
        type=int,
        default=3,
        help="Maximum number of samples to run in parallel (default: 3)"
    )

<<<<<<< HEAD
    parser.add_argument(
        "--python",
        default="/mnt/rds/genetics02/DrummLab/jpf85/.conda/envs/rnaseq/bin/python",
        help="Path to conda environment python (default: rnaseq env)"
    )

    return parser.parse_args()

=======
    return parser.parse_args()


>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
def find_sample_pairs(indir: Path):
    """
    Finds all R1/R2 fastq pairs in the input directory
    Params:
        indir               Path to directory containing fastq files
    Returns:
        list of (r1, r2) Path tuples
    """
    # find all R1 files, ignore non-fastq files
    r1_files = sorted([f for f in indir.glob("*_R1*.fastq*") if f.is_file()])

    if not r1_files:
        raise FileNotFoundError(f"No R1 fastq files found in {indir}")

    pairs = []
    for r1 in r1_files:
        # find matching R2
        r2 = indir / r1.name.replace("_R1", "_R2")
        if not r2.exists():
            print(f"WARNING: No R2 found for {r1.name}, skipping")
            continue
        pairs.append((r1, r2))

    return pairs

<<<<<<< HEAD
=======

>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
def setup_logging(root: Path):
    """
    Sets up logging to both console and a launcher log file
    """
    log_file = root / "launcher.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

<<<<<<< HEAD
def run_sample(root: Path, indir: Path, r1: Path, r2: Path, steps: list, log_dir: Path, python_path: Path):
=======

def run_sample(root: Path, indir: Path, r1: Path, r2: Path, steps: list, log_dir: Path):
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
    """
    Launches main.py for a single sample pair as a background subprocess
    Params:
        root            Path to pipeline root
        indir           Path to input directory
        r1              Path to R1 fastq file
        r2              Path to R2 fastq file
        steps           list of steps to run
        log_dir         Path to directory to write sample log to
<<<<<<< HEAD
        python_path     Path to the conda env's python
=======
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
    Returns:
        subprocess.Popen object
    """
    # get sample name for log file
    sample_name = r1.name.split("_R1")[0]
    log_file = log_dir / f"{sample_name}_run.log"

    # build command
    cmd = [
<<<<<<< HEAD
        str(python_path), str(root / "scripts" / "main.py"),
=======
        "python", str(root / "scripts" / "main.py"),
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
        "--root", str(root),
        "--indir", str(indir),
        "--sample1", r1.name,
        "--sample2", r2.name,
        "--steps", *steps
    ]

<<<<<<< HEAD
    # pass conda bin dir in PATH so fastp, STAR etc... are found in the proper location
    conda_bin = str(python_path.parent)
    env = os.environ.copy()
    env["PATH"] = conda_bin + ":" + env.get("PATH","")

=======
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
    # open log file for this sample
    log_handle = open(log_file, "w")

    # launch subprocess in background
    process = subprocess.Popen(
        cmd,
        stdout=log_handle,
<<<<<<< HEAD
        stderr=log_handle,
        env=env
=======
        stderr=log_handle
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
    )

    return process, log_handle, sample_name

<<<<<<< HEAD
=======

>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
def launcher():
    """
    Main launcher function - finds all sample pairs and runs them in parallel batches
    """

    # ---------------------------------------------------
    # Parse args and load config
    # ---------------------------------------------------

    args = parse_args()
    root = Path(args.root)
    indir = Path(args.indir)
    cfg = ConfigLoader(root / "config.yaml")
<<<<<<< HEAD
    python_path = Path(args.python)
=======
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6

    # ---------------------------------------------------
    # Set up logging
    # ---------------------------------------------------

    logger = setup_logging(root)
    logger.info("=" * 60)
    logger.info("RNAseq Bulk Local Launcher")
    logger.info("=" * 60)

    # ---------------------------------------------------
    # Find all sample pairs
    # ---------------------------------------------------

    pairs = find_sample_pairs(indir)
    total = len(pairs)
    logger.info(f"Found {total} sample pairs to process")
    logger.info(f"Steps: {args.steps}")
    logger.info(f"Max parallel jobs: {args.max_parallel}")
    logger.info("-" * 60)

    # ---------------------------------------------------
    # Set up log directory for sample logs
    # ---------------------------------------------------

    run_name = cfg.get("project", "name")
    run_dir = root / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    sample_logs_dir = run_dir / "sample_logs"
    sample_logs_dir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------
    # Run samples in parallel batches
    # ---------------------------------------------------

    # track active jobs as list of (process, log_handle, sample_name)
    active_jobs = []
    completed = 0
    failed = 0
    pair_queue = list(pairs)

    while pair_queue or active_jobs:

        # fill up slots with new jobs if we have capacity
        while pair_queue and len(active_jobs) < args.max_parallel:
            r1, r2 = pair_queue.pop(0)
            sample_name = r1.name.split("_R1")[0]
            logger.info(f"Starting sample {completed + failed + len(active_jobs) + 1}/{total}: {sample_name}")

            process, log_handle, name = run_sample(
<<<<<<< HEAD
                root, indir, r1, r2, args.steps, sample_logs_dir, python_path
=======
                root, indir, r1, r2, args.steps, sample_logs_dir
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
            )
            active_jobs.append((process, log_handle, name))

        # check status of active jobs
        still_running = []
        for process, log_handle, name in active_jobs:
            # poll to see if process is done (returns None if still running)
            retcode = process.poll()

            if retcode is None:
                # still running
                still_running.append((process, log_handle, name))
            elif retcode == 0:
                # finished successfully
                log_handle.close()
                completed += 1
                logger.info(f"✓ Completed: {name} ({completed}/{total})")
            else:
                # finished with error
                log_handle.close()
                failed += 1
                logger.error(f"✗ FAILED: {name} (return code {retcode}) - check {sample_logs_dir}/{name}_run.log")

        active_jobs = still_running

        # wait a bit before checking again to avoid busy-waiting
        if active_jobs or pair_queue:
            time.sleep(30)

    # ---------------------------------------------------
    # Summary
    # ---------------------------------------------------

    logger.info("=" * 60)
    logger.info(f"All samples processed!")
    logger.info(f"  Completed: {completed}/{total}")
    logger.info(f"  Failed:    {failed}/{total}")
    logger.info("=" * 60)

<<<<<<< HEAD
    # check how many samples failed
=======
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6
    if failed > 0:
        logger.warning(f"{failed} samples failed - check logs in {sample_logs_dir}")
    else:
        logger.info("All samples completed successfully!")

<<<<<<< HEAD
    # plot pca if successful
    if completed > 0:
        logger.info("Running count summarization and PCA...")
        counter = Counts(root, cfg)
        counter.preprocess_pipeline()
    else:
        logger.error("No samples completed - skipping summarization")
=======
>>>>>>> 4b05d6bb8a238ecc81539a6154592e6f800d1ff6

if __name__ == "__main__":
    launcher()
