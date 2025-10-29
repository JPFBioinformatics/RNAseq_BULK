from src.config_loader import ConfigLoader
from pathlib import Path
import argparse

def parse_args():
    """
    Accepts CLI arguments passed after caling main.py
    ex: python main.py --config /path/to/config.yaml --sample Sample1 --steps trim
        will run the "trim" section of main.py on sample 1
        allows for easier parallelization on HPC cluster
    """
    parser = argparse.ArgumentParser(description="Rna-seq Bulk")

    parser.add_argument(
        "--root",
        required=True,
        help="Path to root directory of pipeline, typically /home/user/RNAseq_bulk"
    )

    parser.add_argument(
        "--indir",
        requried=True,
        help="Path to directory where input files are stored"
    )

    parser.add_argument(
        "--sample1",
        required=False,
        help="Path to the forward read file of single sample you want to be processed"
    )

    parser.add_argument(
        "--sample2",
        required=False,
        help="Path to the reverse read file of single sample you want to be processed"
    )

    parser.add_argument(
        "--steps",
        required=False,
        nargs="+",
        default=["trim","align","count"],
        help="Piepeline steps to run (trim, align, or count)"
    )

    parser.add_argument(
        "--tempDir",
        required=True,
        help="Temporary directory to store intermediate fiels in, scratch on HPC clusters"
    )

    return parser.parse_args()

def main():

    # ---------------------------------------------------
    # Parse CLI arguments and handle config.yaml
    # ---------------------------------------------------

    # parse args
    args = parse_args()
    # load config
    cfg =  ConfigLoader(Path(args.root / "config.yaml"))


    # ---------------------------------------------------
    # Set up directories
    # ---------------------------------------------------

    # get run name
    run_name = cfg.get("project","name")

    # create directory to store outputs from this pipeline run
    root_dir = Path(args.root)
    run_dir = root_dir / run_name
    run_dir.mkdir(parents=True,exist_ok=True)

    # create reference dir
    ref_dir = run_dir / "reference"
    ref_dir.mkdir(parents=True,exist_ok=True)
    
    # create all directories for outputs
    data_dirs = cfg.get("data_dirs")
    data_keys = list(data_dirs.keys())
    # create a dir for each key in data_dirs of config.yaml
    for key in data_keys:
        dir_path = cfg.get_path("data_dirs", key, base_path=run_dir)
        dir_path.mkdir(parents=True,exist_ok=True)


    # ---------------------------------------------------
    # Get sample(s) to process
    # ---------------------------------------------------

    # path to dir containing raw data
    in_path = Path(args.indir)

    # sample list
    samples = []

    # check if specific forward and reverse reads are given
    if args.sample1 and args.sample2:
        samples.append(in_path / Path(args.sample1))
        samples.append(in_path / Path(args.sample2))
        
    # if just one sample is given then use that one to find the second file (this is default action)
    elif args.sample1:
        samples.append(in_path / Path(args.sample1))
        samples.append(in_path / Path(str(args.sample1).replace("_R1","_R2")))

    # if samples not specified just grab all forward reads in the file
    else:

        # if indir is a diretory grab all valid files
        if in_path.is_dir():
            samples = list(in_path.glob("*_R1*.fastq*"))
        else:
            raise FileNotFoundError(f"No directory found at --indir {in_path}")
        
        # if there were no valid files in directory raise error
        if not samples:
            raise FileNotFoundError(f"No valid fastq (compressed or not) found at --indir {in_path}")




if __name__ == "__main__":
    main()