from pathlib import Path
from src.utils import log_subprocess, find_name
from src.config_loader import ConfigLoader
import subprocess

class STARIndexBuilder:
    """
    Class to build the reference idnex from genome FastA file and GTF annotation
    """

    def __init__(self, root: Path, temp_dir: Path):
        """
        Params:
            root                                Path to root folder (RNAseq_BULK)
            tempDir                             Path to temporary diretory for intermediate file storage
        """
        self.root = Path(root)
        self.config = self.root / "config.yaml"
        self.temp_dir = Path(temp_dir)

    def build_index(self):
        """
        Creates a star index directory
        Returns:
            path to star index directory
        """
        # load configs
        cfg = ConfigLoader(self.config)

        # get dirs
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)
        idx_dir = cfg.get_path("reference","star_index",base_path=ref_dir)

        # make sure directories exist
        for dir in [ref_dir, idx_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get other values
        genome_fasta = cfg.get_path("reference","genome_fasta",base_path=ref_dir)
        gtf_file = cfg.get_path("reference","gtf_file",base_path=ref_dir)
        threads = cfg.get_threads("STAR")

        # build command
        cmd = [
            "STAR",
            "--runThreadN", str(threads), 
            "--runMode", "genomeGenerate",
            "--genomeDir", str(idx_dir),
            "--genomeFastaFiles", str(genome_fasta),
            "--sjbdGTFfile", str(gtf_file)
        ]

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result, idx_dir, "STARIndexBuilder")

class STARAligner:
    """
    Class to align trimmed fastq files to reference index built by STARIndexBuilder, produces a single BAM file from forward and reverse reads
    """

    def __init__(self, root: Path, temp_dir: Path):
        """
        Params:
            root                                Path to root directory (RNAseq_BULK)
            temp_dir                            path to temporary directory for intermediate files
        """
        self.root = Path(root)
        self.config = self.root / "config.yaml"
        self.temp_dir = temp_dir

    def align(self, r1: Path, r2: Path):
        """
        Preforms alignment of file r1 and r2 to the Star index

        Params:
            r1                                  Path object pointing to the trimmed forward read to map
            r2                                  Path object pointing to the trimmed reverse read to map
            sample_name                         Name to save combined reads under
        """
        # get sample name
        name = find_name(r1,r2)

        # connect to config
        cfg = ConfigLoader(self.config)

        # get directories
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = project / name
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)
        star_index = cfg.get_path("reference","star_index",base_path=ref_dir)
        temp_dir = self.temp_dir / {name}

        # make sure directories exist
        for dir in [project,sample_dir,ref_dir,star_index,temp_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get threads
        threads = cfg.get_threads("STAR")

        # get additional parameters
        readFilesCommand = cfg.get("tools","STAR","file_type")
        outSAMtype = cfg.get("tools","STAR","outSAMtype")
        genomeload = cfg.get("tools","STAR","genomeload")
        outFilterMultimapMax = cfg.get("tools","STAR","outFilterMultimapMax")
        twopassMode = cfg.get("tools","STAR","twopassMode")
        sjdbOverhang = cfg.get("tools","STAR","sjdbOverhang")
        alignIntronMax = cfg.get("tools","STAR","alignIntronMax")
        outReadsUnmapped = cfg.get("tools","STAR","outReadsUnmapped")

        # specify output file
        out_file = sample_dir / f"{name}_"

        # build command
        cmd = [
            "STAR",
            "--runThreadN", str(threads),
            "--genomeDir", str(star_index),
            "--readFilesIn", str(r1), str(r2),
            "--readFilesCommand", str(readFilesCommand),
            "--outFileNamePrefix", str(out_file),
            "--outSAMtype", str(outSAMtype),
            "--gneomeLoad", str(genomeload),
            "--outFilterMultimapNmax", str(outFilterMultimapMax),
            "--twopassMode", str(twopassMode),
            "--sjdbOverhang", str(sjdbOverhang),
            "--alignIntronMax", str(alignIntronMax),
            "--outTmpDir", str(temp_dir)
        ]

        if outReadsUnmapped == "Within" or outReadsUnmapped == "Fastx":
            cmd.extend(["--outReadsUnmapped", outReadsUnmapped])

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result, sample_dir, "STARAligner")