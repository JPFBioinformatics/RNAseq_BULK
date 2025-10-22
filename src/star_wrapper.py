from pathlib import Path
from src.utils import log_subprocess
import subprocess

class STARIndexBuilder:
    """
    Class to build the reference idnex from genome FastA file and GTF annotation
    """

    def __init__(self,genome_fasta: Path, gtf_file: Path, index_dir: Path, log_dir: Path, threads: int=4):
        """
        Params:
            genome_fasta                        Path object pointing to the genome fasta
            gtf_file                            Path object pointing to the gtf file
            index_dir                           Path object pointing to the output directory to build star index files in (star_index dir in reference)
            log_dir                             Path object pointing to the directory to store subprocess logs in (/run_name/logs)
            threads                             Number of CPU threads to assign to task
        """
        self.genome_fasta = Path(genome_fasta)
        self.gtf_file = Path(gtf_file)
        self.index_dir = Path(index_dir)
        self.log_dir = Path(log_dir)
        self.threads = threads

        # build output directory (/root/run_name/reference/star_index/)
        self.index_dir.mkdir(parents=True,exist_ok=True)

    def build_index(self):
        """
        Creates a star index directory
        """

        # build command
        cmd = [
            "STAR",
            "--runThreadN", str(self.threads),
            "--runMode", "genomeGenerate",
            "--genomeDir", str(self.index_dir),
            "--genomeFastaFiles", str(self.genome_fasta),
            "--sjbdGTFfile", str(self.gtf_file)
        ]

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result, self.log_dir, "Star Reference Index Production", "STARIndexBuilder")

class STARAligner:
    """
    Class to align trimmed fastq files to reference index built by STARIndexBuilder, produces a single BAM file from forward and reverse reads
    """

    def __init__(self, index_dir: Path, out_dir: Path, log_dir: Path, threads: int=4):
        """
        Params:
            index_dir                           Path object pointing to the star_index diretctory to map reads to
            out_dir                             Path object pointing to where the BAM files are to be stored, temporary storage space due to large size of BAM files
            log_dir                             Path object pointing to where the subprocess logs will be stored (/run/logs)
            threads                             Number of CPU threads to assign to task
        """
        self.index_dir = index_dir
        self.out_dir = out_dir
        self.log_dir = log_dir
        self.threads = threads

    def align(self, r1: Path, r2: Path, sample_name: str):
        """
        Preforms alignment of file r1 and r2 to the Star index

        Params:
            r1                                  Path object pointing to the trimmed forward read to map
            r2                                  Path object pointing to the trimmed reverse read to map
            sample_name                         Name to save combined reads under
        """

        # output file location/name
        out_prefix = self.out_dir / f"{sample_name}_"

        # build command
        cmd = [
            "STAR",
            "--runThreadN", str(self.threads),
            "--genomeDir", str(self.index_dir),
            "--readFilesIn", str(r1), str(r2),
            "--readFilesCommand", "zcat",         # zcat for zipped fastq files
            "--outFileNamePrefix", str(out_prefix),
            "--outSAMtype", "BAM", "SotredByCoordinate"
        ]

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result, self.log_dir, f"{sample_name}", "STARAligner")