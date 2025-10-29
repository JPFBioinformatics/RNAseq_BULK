# uses FastP to trim data and produce a pre-trim QC report and a post-trim QC report
import subprocess
from pathlib import Path
from src.utils import log_subprocess, find_name, check_bool
from src.config_loader import ConfigLoader

class QCTrimmer:
    """
    Class to run FastQC on one or more FastQ files
    """

    def __init__(self, root: Path, temp_dir: Path):
        """
        Params:
            root                        Path to the root file of the pipeline (RNAseq_bulk)
            temp_dir                    Path to the temp dir for intermediate files
        """
        self.root = Path(root)
        self.config = root / "config.yaml"
        self.temp_dir = Path(temp_dir)

    def run_fastp(self, r1_in: Path, r2_in: Path):
        """
        Runs FastP QC, trim, second QC and stores files as specified
        Params:
            r1_in                       path to the .gz forward read file
            r2_in                       path to the .gz reverse read file
        Returns:
            r1_out, r2_out paths to the trimmed r1 and r2 fastq.gz files
        """
        
        # get names for the files/sample
        name = find_name(r1_in,r2_in)

        #load config
        cfg = ConfigLoader(self.root / "config.yaml")

        # get other dir paths
        project = cfg.get_path("project","name", base_path=self.root)
        sample_dir = project / name
        temp_dir = self.temp_dir / name

        # build the directories if they do not already exist
        for dir in [project,sample_dir,temp_dir]:
            dir.mkdir(parents=True,exist_ok=True)
            
        # get other values
        threads = cfg.get_threads("fastp")
        length_required = cfg.get("tools","fastp","length_required")
        qualified_quality_phred = cfg.get("tools","fastp","qualified_quality_phred")
        specifyAdapter = cfg.get("tools", "fastp", "specify_adapter")

        # build output files
        r1_out = temp_dir / "trimmed_R1.fastq.gz"
        r2_out = temp_dir / "trimmed_R2.fastq.gz"
        html_out = sample_dir / "fastP_QC.html"
        json_out = sample_dir / "fastP_QC.json"

        # build fastp command, ensuring that any Path objects are converted to str
        cmd = [
            "fastp",
            "-i", str(r1_in),
            "-I", str(r2_in),
            "-o", str(r1_out),
            "-O", str(r2_out),
            "-h", str(html_out),
            "-j", str(json_out),
            "-w", str(threads),
            "--length_required", length_required,
            "--qualified_quality_phred", qualified_quality_phred
        ]
        
        # check that bools are bools
        check_bool([specifyAdapter])

        # check if we want to specify adapters
        if specifyAdapter:
            # get adapter sequences
            adapter_sequence = cfg.get("tools","fastp","adapter_sequence")
            adapter_sequence_r2 = cfg.get("tools","fastp","adapter_sequence_r2")
            # include specified adapters in command
            if adapter_sequence:
                cmd.extend(["--adapter_sequence", adapter_sequence])
            if adapter_sequence_r2:
                cmd.extend(["--adapter_sequence_r2", adapter_sequence_r2])

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result, sample_dir, "fastP")
        
        # return location of temp trimmed files
        return r1_out, r2_out
