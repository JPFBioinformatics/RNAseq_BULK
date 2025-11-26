# region Imports

from pathlib import Path
import sys, argparse,subprocess,shlex

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.config_loader import ConfigLoader
from src.utils import log_subprocess, get_max_threads, get_total_memory

# endregion

def parse_args():
    """
    accepts CLI arguments for launcher.py 
    """
    parser = argparse.ArgumentParser(description="Launcher to submit RNA-seq HPC pipeline as a SLURM array job")

    parser.add_argument(
        "--root",
        required=True,
        help = "path to root dir of pipeline"
    )

    parser.add_argument(
        "--indir",
        required=True,
        help="directory containing input fastq files"
    )

    parser.add_argument(
        "--steps",
        nargs ="+",
        default=["trim","align","count"],
        help="pipeline steps to execute"
    )

    parser.add_argument(
        "--runScript",
        default="run_pipeline.sh",
        help="SLURM wrappr script to run"
    )

    return parser.parse_args()

def launcher():
    """
    launches the run_pipeline.sh file with dynamically updated values for SLURM jobs
    """

    # ---------------------------------------------------
    # Parse CLI arguments and load config
    # ---------------------------------------------------

    # set up args/config and make sure config is fine
    args = parse_args()
    cfg = ConfigLoader(Path(args.root) / "config.yaml")
    cfg.check_bools()

    # get run name and handle run dir
    name = str(cfg.get("project","name"))
    jobname = shlex.quote(name)
    project_dir =  cfg.get_path("project","name",base_path=args.root)
    project_dir.mkdir(parents=True,exist_ok=True)

    # count samples
    r1_files = sorted(Path(args.indir).glob("*_R1*.fastq*"))
    num_samples = len(r1_files)
    if num_samples == 0:
        raise FileNotFoundError(f"No R1 fastq files found in {args.indir}")

    # get cfg values
    max_threads = get_max_threads(cfg)
    total_mem = get_total_memory(cfg,max_threads)
    conda_env = cfg.get("project","conda_env")

    # run script path
    run_script_path = shlex.quote(str(Path(args.root) / args.runScript))

    # get steps
    steps = shlex.quote(' '.join(args.steps))

    # ---------------------------------------------------
    # build sbatch command
    # ---------------------------------------------------
    
    cmd = (
        f"sbatch "
        f"--array=0-{num_samples - 1} "
        f"--cpus-per-task={max_threads} "
        f"--mem={total_mem} "
        f"--job-name={jobname} "
        f"{run_script_path} "
        f"--root {shlex.quote(args.root)} "
        f"--indir {shlex.quote(args.indir)} "
        f"--steps {steps} "
        f"--conda-env {shlex.quote(conda_env)}"
    )

    # run subprocess and log
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    log_subprocess(result, project_dir, "launcher")

if __name__ == "__main__":
    launcher()