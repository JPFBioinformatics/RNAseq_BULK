from pathlib import Path
from src.utils import log_subprocess
import subprocess
from config_loader import ConfigLoader

# sort, index, filter in that order


class SamtoolsWrapper:
    """
    Wrapper class for common samtools operations on BAM files:
        -sort
        -index
        -filter
        -flatstat (QC for logs)
    """

    def __init__(self, log_dir: Path, config: Path):
        """
        Params:
            log_dir                         Path to the directory to store subprocess logs (/run_name/logs)
            config                          Path to config.yaml
        """
        self.log_dir = log_dir
        self.config = config

    def sort_bam(self, bam_file: Path):
        """
        Sorts given BAM file to prepare for indexing, first step of samtools processing
        Params:
            bam_file                        bam file to be sorted
        """
        cfg = ConfigLoader(self.config)


        # build output file name
        out_name = f"{bam_file.stem}_sorted.bam"

        # build command
        cmd = [
            "samtools",
            "sort",
            "-@", str(self.threads),
            "-o", str(out_name), str(bam_file)
        ]

        # execute command
        result = subprocess(cmd, capture_output=True, text=True)

        # log results
        log_subprocess(result, self.log_dir, f"{bam_file.stem} sort", "samtools")

    def index_bam(self, bam_file: Path):
        """
        Indexes given sorted bam file and produces a bam_file.bam.bai comapinon file in same directory as bam_file, second step of samtools processing
        Params:
            bam_file                        Path to sorted bam file to be indexed
        """

        # build command
        cmd = [
            "samtools",
            "index", str(bam_file)
        ]

        # run command
        result = subprocess(cmd, capture_output=True, text=True)

        # log results
        log_subprocess(result, self.log_dir, f"{bam_file.stem} indexing", "samtools")

    def filter_bam(self, min_mapq: int=10, remove_unmapped: bool=True):
        """
        Filters input bam file to remove maps of poor quality
        Params:
            min_mapq                        minimum quality of mapping to keep the read, at default (10) it only keeps reads >90% likely to 
        """
