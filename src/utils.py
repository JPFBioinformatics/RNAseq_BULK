from pathlib import Path
import subprocess
from datetime import datetime
import json
from src.config_loader import ConfigLoader
import os

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
        "stderr": result.stderr
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
        name                        Name of the sample
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

def generate_paired_samples(indir: Path, sample1=None, sample2=None):
    """
    Generates a list of paired samples from input directory
    Params:
        indir                           Path to directory containing fastq files
        sample1                         Optional, name of/Path to forward read file if you want to run on a single pair (include extensions)
        sample2                         Optional, name of/Path to reverse read file if you want to run on a single pair (include extensions)
    Returns:
        paired_samples                  list of tuples [(r1_1,r2_1),(r1_2,r2_2)...] where each entry is a pair of Fastq files to process
    """
    in_path = Path(indir)
    paired_samples = []

    # helper method to ensure Path object is absolute
    def ensure_path(f):
        p = Path(f)
        if not p.is_absolute():
            p = in_path / p
        return p

    # if the user specifies a single pair to process
    if sample1 and sample2:
        r1 = ensure_path(sample1)
        r2 = ensure_path(sample2)
        # raise error if files do not exist
        if not r1.exists() or not r2.exists():
            raise FileNotFoundError(f"Could not find specified files: {r1} and {r2}")
        # append to paired_samples
        paired_samples.append((r1,r2))
        
    # if the user specifies just one sample then generate its pair
    elif sample1:
        r1 = ensure_path(sample1)
        r2 = ensure_path(str(r1.name).replace("_R1","_R2"))
        # raise error if fiels do not exist
        if not r1.exists() or not r2.exists():
            raise FileNotFoundError(f"Could not find specified files: {r1} and {r2}")
        # append to paired samples
        paired_samples.append((r1,r2))

    # if no specific samples are given
    else:
        # grab all r1 files
        r1_files = sorted(in_path.glob("*_R1*.fastq*"))
        # raise error if R1 files do not exist
        if not r1_files:
            raise FileNotFoundError(f"no R1 files found in {in_path}")
        # find the partner of each R1 file and append
        for r1 in r1_files:
            r2 = in_path / r1.name.replace("_R1","_R2")
            # check to make sure R2 exists
            if not r2.exists():
                raise FileNotFoundError(f"Missing R2 file for sample {r1.name}")
            # append to output list
            paired_samples.append((r1,r2))

    return paired_samples

def get_max_threads(cfg: ConfigLoader):
    """
    gets max threads specified in config.yaml
    Params:
        cfg                         ConfigLoader object
    Returns:
        max_threads                 maximum number of threads needed for the job
    """
    # get max threads needed to pass on to the sbatch command
    tools_dir = cfg.get("tools")
    tools = tools_dir.keys()

    # find max threads
    max_threads = 0
    for key in tools:
        threads = cfg.get("tools",str(key),"threads")
        if threads > max_threads:
            max_threads = threads

    return max_threads

def get_total_memory(cfg: ConfigLoader, max_threads: str):
    """
    calculates max memory required for the job
    Params:
        cfg                         ConfigLoader object
        max_threads                 max_threads needed for the SLURM job
    Returns:
        total_mem                   total memory needed for the job based on samtools sort mem specified
    """
    # get sort memory needed for samtools
    samtools_sort_mem = cfg.get("tools","samtools","sortMemory")

    # split value into number and units
    num_part = samtools_sort_mem[:-1]
    unit = samtools_sort_mem[-1].upper()

    # check unit value
    if unit not in ("G","M"):
        raise ValueError(f"Unknown memory unit: {unit} must be G (gigabytes) or M (megabytes)")
    # convert value to float
    try:
        mem_val = float(num_part)
    except ValueError:
        raise ValueError(f"could not convert memory value {mem_val} to number")
    # calculate total memory value needed

    mem_val *= int(max_threads)
    # convert to string total memory
    total_mem = f"{int(mem_val)}{unit}"

    return total_mem

def get_scratch(root: Path):
    """
    gets scratch/temp dir from HPC or as specified by the CLI command if running locally
    Params:
        root                        Path to root dir of project to put temp in if not running on cluster with Scratch
    """

    # choose node-local SSD if available
    tmp = os.environ.get("TMPDIR")
    if tmp:
        return Path(tmp)
    
    # fallback to glaobal if local not available
    pfs = os.environ.get("PFSDIR")
    if pfs:
        return Path(pfs)
    
    # if neither are available then make temp dir in root
    return Path(root) / "temp"

def get_STAR_suffix(cfg: ConfigLoader):
    """
    Gets the suffxi that STAR will add to the bam file generated by Star_wrapper STARAligner.align()
    Params:
        cfg                     ConfigLoader object with config.yaml loaded
    Returns:
        full_suffix             the full suffix STAR will add to bam file when it is aligned
    """
    # get the sample type configuration
    SAMtype = cfg.get("tools","STAR","outSAMtype")

    # split string at the space
    parts = str(SAMtype).split()
    file_type = parts[0]
    sorting = parts[1] if len(parts)>1 else None

    # find file type for suffix
    if file_type == "BAM":
        suffix_type = "bam"
    elif file_type == "SAM":
        suffix_type = "sam"
    else:
        raise ValueError(f"Unsupported STAR outSAMtype {SAMtype}")
    
    # check if sorting was specified
    if not sorting:
        full_suffix = f"_Aligned.out.{suffix_type}"

    # if it was then add requried parts to STAR suffix
    else:
        if sorting not in ["SortedByCoordinate","SortedByName"]:
            raise ValueError(f"Unsupported STAR outSAMtype {SAMtype}")

        elif sorting == "SortedByCoordinate":
            suffix_sort = "sortedByCoord"
        else:
            suffix_sort = "sortedByName"

        full_suffix = f"_Aligned.{suffix_sort}.out.{suffix_type}"

    return full_suffix

    