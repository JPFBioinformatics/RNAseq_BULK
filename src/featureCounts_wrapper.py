# region Imports

from pathlib import Path
import sys, subprocess

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.config_loader import ConfigLoader
from src.utils import log_subprocess

# endregion

class FeatureCountsWrapper:
    """
    wrapper class for featureCounts that produce text files for downstream analysis
    """
    def __init__(self, cfg: ConfigLoader, root: Path, sample_dir: Path):
        """
        Params:
            cfg                         ConfigLoader object that has loaded config.yaml for this project
            root                        Path to root project folder
            sample_dir                  Directory to store sample data in
        """
        self.root = root
        self.cfg = cfg
        self.sample_dir = Path(sample_dir)

    def count_features(self, bam_file: Path):
        """
        Runs featurecounts to produce count files from bam files
        Params:
            bam_file                    path to the bam file to be counted
        """
        # sample name
        name = bam_file.stem.split("_Aligned")[0]

        # connect to config
        cfg = self.cfg

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)
        sample_dir = self.sample_dir
        # make sure dirs exist
        for dir in [project,ref_dir,sample_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get other values
        gtf_file = cfg.get_path("reference","gtf_file",base_path=ref_dir)
        threads = cfg.get_threads("featureCounts")
        strand = cfg.get("tools","featureCounts","strand_specific")
        feature_type = cfg.get("tools","featureCounts","feature_type")
        gtf_attr_type = cfg.get("tools","featureCounts","gtf_attr_type")
        ignoreDup = cfg.get("tools","featureCounts","ignoreDup")
        isPairedEnd = cfg.get("tools","featureCounts","isPairedEnd")
        save_files = cfg.get("project","save_files")


        # build output file
        out_file = sample_dir / f"{name}_counts.txt"

        # build command
        cmd = [
            "featureCounts",
            "-T", str(threads),
            "-s", str(strand),
            "-t", str(feature_type),
            "-g", str(gtf_attr_type),
            "-a", str(gtf_file),
            "-o", str(out_file)
        ]

        # append optional values
        if ignoreDup:
            cmd.append("--ignoreDup")
        if isPairedEnd:
            cmd.append("-p")

        # append the input file
        cmd.append(str(bam_file))

        # run subprocess
        result = subprocess.run(cmd,capture_output=True,text=True)

        # log subprocess
        log_subprocess(result,sample_dir,"featureCounts")
        
        # if we are not saving files then delete the input bam
        if out_file.exists() and not save_files:
            try:
                bam_file.unlink()
                print(f"FeatureCount successful, deleted input bam:\n{bam_file}\n")
            except Exception as e:
                print(f"Warning, could not delete file:\n{bam_file}\nError:\n{e}\n")
        

