from pathlib import Path
import subprocess
import os
from datetime import datetime
import json
from src.config_loader import ConfigLoader

def log_subprocess(result: subprocess.CompletedProcess, log_dir: Path, step: str):
    """
    Collects logs produced by python subprocess and puts them in a text file for easy viewing
    Params:
        result:                     subprocess.CompletedProcess object that can be logged, result of running subprocess
        log_dir:                    path to the sample file where we will write log file
        step:                       step (fastp, starAlign etc...) that this subprocess was associated with
    """
    # get timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # make sure log_dir exists
    log_dir.mkdir(parents=True,exist_ok=True)

    # path to log file
    log_file = log_dir / "subprocess_log.jsonl"

    # dict of values to store
    data = {
        "step": step,
        "log_ts": timestamp,
        "stdout": result.stdout,
        "stedrr": result.stderr
    }

    with open(log_file, "a") as f:
        f.write(json.dumps(data) + "\n")

def find_name(r1: Path, r2: Path):
    """
    Generates sring names for two given input files under standard naming conventions for this pipeline
    Params:
        r1                          Path object pointing to the location of the forward read file including suffix _R1 before file extensions
        r2                          Same as r1 but for the reverse read file
    Returns:
        three strings for sample names
        name                        Name of the sample
        name1                       R1 forward read name
        name2                       R2 reverse read name
    """

    # make sure R1 and R2 are paths
    r1 = Path(r1)
    r2 = Path(r2)

    # get shared smaple name
    name1 = r1.name.split("_R1")[0]
    name2 = r2.name.split("_R2")[0]

    # if the user followed naming conventions in README.txt then get the name of the sample
    if name1 == name2:
        name = name1
    # else generate a new name
    else:
        name = f"{name1}_{name2}"

    return name

def remove_file_extensions(file: Path):
    """
    stems a path object until there are no more file extensiosn present
    Params:
        file                        Path object specifying a file that you want to find the name of without file extensinos (.gz, .fastq, .bam etc...)
    Returns:
        string file name without extensions
    """
    # ensure file is a Path object
    file = Path(file)

    # remove suffix while one exists
    while file.suffix:
        file = file.with_suffix("")

    return file.name

def check_name(file: Path):
    """
    Checks to make sure that a file does not have a '.' in an invalid location (not associated with an extension)
    used to check .fastq/.fastq.gz input files, not for later in pipeline
    """

    # make sure file is a Path
    file = Path(file)

    # get the file's name
    name = file.name

    # get the portion before the extensions
    if "_R1" in name:
        name = name.split("_R1")[0]
    if "_R2" in name:
        name = name.split("_R2")[0]

    # check for incorrect '.' placement
    if "." in name:
        raise ValueError(f"Filename {name} contains a '.' character not associated witha file extension, remove it and try again")

def check_bool(list: list):
    """
    checks if input list are all boolean
    Params:
        list                        list of values to make sure are booleans
    """
    # check entries to ensure they are booleans
    for entry in list:
        if not isinstance(entry,bool) or entry is None:
            raise ValueError(f"{entry} is not boolean, check config.yaml values must be true, True, false, or False")
