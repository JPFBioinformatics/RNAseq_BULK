# uses FastP to trim data and produce a pre-trim QC report and a post-trim QC report
import subprocess
from pathlib import Path
from src.utils import log_subprocess
from src.config_loader import ConfigLoader

class QCTrimmer:
    """
    Class to run FastQC on one or more FastQ files
    """

    def __init__(self, root: Path):
        """
        Params:
            root                        Path to the root file of the pipeline (RNAseq_bulk)
        """
        self.root = root

    def run_fastp(self, r1_in: Path, r2_in: Path, adapt_fwd: str=None, adapt_rev: str=None):
        """
        Runs FastP QC, trim, second QC and stores files as specified
        Params:
            r1_in:                      path to the .gz forward read file
            r2_in:                      path to the .gz reverse read file
        """
        # get root path
        root = self.root
        #load config
        cfg = ConfigLoader(root / "config.yaml")
        # get other dir paths
        project = self.root / cfg.get("project","name")
        qc_path= project / cfg.get("project","data_dirs","qc_report")
        trimmed = project/ cfg.get("project","temp_dirs","trimmed")
        
        # get threads for this operation
        threads = cfg.get("fastp", "threads")

        # get names for the input samples
        name1 = r1_in.stem
        name2 = r2_in.stem
        # get full run name

        # build fastp command, ensuring that any Path objects are converted to str
        cmd = [
            "fastp",
            "-i", str(r1_in),
            "-I", str(r2_in),
            "-o", str(trimmed / f"{name1}_trimmed.fastq.gz"),
            "-O", str(trimmed / f"{name2}_trimmed.fastq.gz"),
            "-h", str(qc_path / f"{name1}.html"),
            "-j", str(qc_path / f"{name2}.json"),
            "-w", str(threads)
        ]

        # include specified adapters, if not default to auto detect
        if adapt_fwd:
            cmd += ["--adapter_sequence", adapt_fwd]
        if adapt_rev:
            cmd += ["--adapter_sequence_r2", adapt_rev]

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # get location of log directory
        log_dir = self.log_dir
        log_dir.mkdir(parents=True,exist_ok=True)

        # get name of these samples
        r1_name = r1_in.name
        r2_name = r2_in.name
        # build header for log file
        header = f"Samples: {r1_name} and {r2_name}"

        # log subprocess
        log_subprocess(result, log_dir, header, "fastp")
