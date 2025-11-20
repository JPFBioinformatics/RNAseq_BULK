# region Imports

from pathlib import Path
import subprocess, sys

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.utils import log_subprocess, find_name, get_STAR_suffix
from src.config_loader import ConfigLoader

# endregion

class STARIndexBuilder:
    """
    Class to build the reference idnex from genome FastA file and GTF annotation
    """

    def __init__(self, cfg: ConfigLoader, root: Path, temp_dir: Path):
        """
        Params:
            cfg                                 ConfigLoader object that has loaded config.yaml for this project
            root                                Path to root folder (RNAseq_BULK)
            tempDir                             Path to temporary diretory for intermediate file storage
        """
        self.root = Path(root)
        self.cfg = cfg
        self.temp_dir = Path(temp_dir)

    def build_index(self):
        """
        Creates a star index directory
        Returns:
            path to star index directory
        """
        # load configs
        cfg = self.cfg

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
            "--sjdbGTFfile", str(gtf_file)
        ]

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result, idx_dir, "STARIndexBuilder")

class STARAligner:
    """
    Class to align trimmed fastq files to reference index built by STARIndexBuilder, produces a single BAM file from forward and reverse reads
    """

    def __init__(self, cfg: ConfigLoader, root: Path, temp_dir: Path):
        """
        Params:
            cfg                                 ConfigLoader object that has loaded config.yaml for this project
            root                                Path to root directory (RNAseq_BULK)
            temp_dir                            path to temporary directory for intermediate files
        """
        self.root = Path(root)
        self.cfg = cfg
        self.temp_dir = temp_dir

    def align(self, r1: Path, r2: Path, cleanup=False):
        """
        Preforms alignment of file r1 and r2 to the Star index

        Params:
            r1                                  Path object pointing to the trimmed forward read to map
            r2                                  Path object pointing to the trimmed reverse read to map
            cleanup                             bool, if true then delete the input r1/r2 files

        Returns:
            path to aligned bam file
        """
        # get sample name
        print(f"\nStar Inputs:\n{r1}    {r2}\n")
        name = find_name(r1,r2)

        # connect to config
        cfg = self.cfg

        # get directories
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = project / name
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)
        star_index = cfg.get_path("reference","star_index",base_path=ref_dir)
        temp_dir = self.temp_dir / name

        # make sure directories exist
        for dir in [project,sample_dir,ref_dir,star_index,temp_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get threads
        threads = cfg.get_threads("STAR")

        # get additional parameters
        readFilesCommand = cfg.get("tools","STAR","file_type")
        outSAMtype = cfg.get("tools","STAR","outSAMtype")
        genomeload = cfg.get("tools","STAR","genomeLoad")
        outFilterMultimapMax = cfg.get("tools","STAR","outFilterMultimapMax")
        twopassMode = cfg.get("tools","STAR","twopassMode")
        sjdbOverhang = cfg.get("tools","STAR","sjdbOverhang")
        alignIntronMax = cfg.get("tools","STAR","alignIntronMax")
        outReadsUnmapped = cfg.get("tools","STAR","outReadsUnmapped")
        outFilterMismatchNoverLmax = cfg.get("tools","STAR","outFilterMismatchNoverLmax")

        # specify output file
        out_file = sample_dir / f"{name}"
        print(f"\naligner out_file:\n{out_file}")

        # build command
        cmd = [
            "STAR",
            "--runThreadN", str(threads),
            "--genomeDir", str(star_index),
            "--readFilesIn", str(r1), str(r2),
            "--readFilesCommand", str(readFilesCommand),
            "--outFileNamePrefix", str(out_file),
            "--outSAMtype", str(outSAMtype),
            "--genomeLoad", str(genomeload),
            "--outFilterMultimapNmax", str(outFilterMultimapMax),
            "--twopassMode", str(twopassMode),
            "--sjdbOverhang", str(sjdbOverhang),
            "--alignIntronMax", str(alignIntronMax),
            "--outTmpDir", str(temp_dir / "STAR")
        ]

        # add outreadsunmapped if specified
        if outReadsUnmapped == "Within" or outReadsUnmapped == "Fastx":
            cmd.extend(["--outReadsUnmapped", outReadsUnmapped])

        # add outFilterMismatchNoverLmax if specified
        if outFilterMismatchNoverLmax:
            cmd.extend(["--outFilterMismatchNoverLmax",outFilterMismatchNoverLmax])
            
        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result, sample_dir, "STARAligner")

        # get STAR suffix
        suffix = get_STAR_suffix(cfg)

        # build full path to output file with STAR suffixes attached
        output_full = out_file.with_name(out_file.name + suffix)

        # delete input trimmed files if successful
        if output_full.exists() and cleanup:
            try:
                r1.unlink()
                r2.unlink()
            except Exception as e:
                print(f"Warning: could not delete input FASTQ files {r1}, {r2}\n{e}")

        # return full path
        return output_full
    