from pathlib import Path
from src.config_loader import ConfigLoader
from src.utils import check_bool,log_subprocess
import subprocess

class FeatureCountsWrapper:
    """
    wrapper class for featureCounts that produce text files for downstream analysis
    """
    def __init__(self, root: Path):
        """
        Params:
            root                        Path to root project folder
        """
        self.root = root
        self.config = root / "config.yaml"

    def count_features(self, bam_file: Path):
        """
        Runs featurecounts to produce count files from bam files
        Params:
            bam_file                    path to the bam file to be counted
        """
        # sample name
        name = bam_file.stem.split("_Aligned")[0]

        # connect to config
        cfg = ConfigLoader(self.config)

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)
        sample_dir = project / name
        # make sure dirs exist
        for dir in [project,ref_dir,sample_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get other values
        gtf_file = cfg.get_path("reference","gtf_file",base_path=ref_dir)
        threads = cfg.get_threads("featureCounts")
        strand = cfg.get("tools","featureCounts","strand_specific")
        feature_type = cfg.get("tools","featureCounts","feature_type")
        gtf_attr_type = cfg.get("tools","featureCounts","gtf_attr_type")
        # bools
        ignoreDup = cfg.get("tools","featureCounts","ignoreDup")
        isPairedEnd = cfg.get("tools","featureCounts","isPairedEnd")
        # ensure booleans are booleans
        check_bool([ignoreDup,isPairedEnd])

        # build output file
        out_file = sample_dir / f"counts.txt"

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

