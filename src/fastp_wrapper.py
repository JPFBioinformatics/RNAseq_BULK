# region Imports

import subprocess, sys
from pathlib import Path

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.utils import log_subprocess, find_name
from src.config_loader import ConfigLoader

# endregion

class QCTrimmer:
    """
    Class to run FastQC on one or more FastQ files
    """

    def __init__(self, cfg: ConfigLoader, root: Path, temp_dir: Path, sample_dir: Path):
        """
        Params:
            cfg                         ConfigLoader object that has loaded the config.yaml file for this project
            root                        Path to the root file of the pipeline (RNAseq_bulk)
            temp_dir                    Path to the temp dir for intermediate files
            sample_dir                  dierctory where sample data is to be stored
        """
        self.root = Path(root)
        self.cfg = cfg
        self.temp_dir = Path(temp_dir)
        self.sample_dir = Path(sample_dir)

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
        cfg = self.cfg

        # get other dir paths
        project = cfg.get_path("project","name", base_path=self.root)
        sample_dir = self.sample_dir
        temp_dir = self.temp_dir / name

        # build the directories if they do not already exist
        for dir in [project,sample_dir,temp_dir]:
            dir.mkdir(parents=True,exist_ok=True)
            
        # get other values
        threads = cfg.get_threads("fastp")
        length_required = cfg.get("tools","fastp","length_required")
        qualified_quality_phred = cfg.get("tools","fastp","qualified_quality_phred")
        specifyAdapter = cfg.get("tools", "fastp", "specify_adapter")
        save_inputs = cfg.get("tools","fastp","save_inputs")

        # build output files
        r1_out = temp_dir / f"{name}_R1_trimmed.fastq.gz"
        r2_out = temp_dir / f"{name}_R2_trimmed.fastq.gz"
        html_out = sample_dir / "logs" / "fastP_QC.html"
        json_out = sample_dir / "logs" / "fastP_QC.json"

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
            "--length_required", str(length_required),
            "--qualified_quality_phred", str(qualified_quality_phred)
        ]

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

        if r1_out.exists() and r2_out.exists() and not save_inputs:
            try:
                r1_in.unlink()
                r2_in.unlink()
                print(f"FastP complete, deleted input files:\n{r1_in}\n{r2_in}\n")
            except Exception as e:
                print(f"Warning, could not delete FastP input files:\n{r1_in}\n{r2_in}\n")

        # return location of temp trimmed files
        return r1_out, r2_out
