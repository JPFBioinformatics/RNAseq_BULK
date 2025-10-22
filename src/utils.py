from pathlib import Path
import subprocess
import os

def log_subprocess(result: subprocess.CompletedProcess, log_dir: Path, header: str, prefix: str):
    """
    Collects logs produced by python subprocess and puts them in a text file for easy viewing
    Params:
        result:                     subprocess.CompletedProcess object that can be logged, result of running subprocess
        log_dir:                    path to the file within log folder to write to, log folder for this run
        header:                     header for this section of the log file created
        prefix:                     prefix for the log file, for example "fastp" for the logs corrosponding to a fastp run
    """
    # make sure log_dir exists
    log_dir.mkdir(parents=True,exist_ok=True)

    # path to log file, creates file if it does not exist
    log_file = log_dir / f"{prefix}.log"

    with open(log_file, "a") as f:
        f.write(f"========================= {header} =========================\n")
        f.write("===== STDOUT =====")
        f.write(result.stdout or "No stdtout\n")
        f.write("===== STDERR =====")
        f.write(result.stderr or "No stderr\n\n")

