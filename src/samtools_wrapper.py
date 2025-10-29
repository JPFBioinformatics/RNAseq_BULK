from pathlib import Path
from src.utils import log_subprocess, check_bool
import subprocess
from config_loader import ConfigLoader

# sort, index, filter in that order


class SamtoolsWrapper:
    """
    Wrapper class for common samtools operations on BAM files:
        -sort
        -filter
        -index
        -flatstat (QC for logs)
    pipeline will do these operations in above order by default
    """

    def __init__(self, root: Path, temp_dir: Path):
        """
        Params:
            root                            path to root folder of project
            temp_dir                        path to the temp dir for intermediate files
        """
        self.root = Path(root)
        self.config = self.root / "config.yaml"
        self.temp_dir = Path(temp_dir)

    def sort_file(self, file: Path):
        """
        Sorts given BAM file to prepare for indexing, first step of samtools processing
        Params:
            file                                bam or cram file to be sorted
        Returns:    
            path to the sorted file (in temp dir)
        """
        # get raw bam QC
        self.flagstat(file,"raw")

        # get sample name
        name = file.stem.split("_Aligned")[0]

        # load config
        cfg = ConfigLoader(self.config)
        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = project / name
        temp_dir = self.temp_dir / name
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)

        # make sur directories exist
        for dir in [project,sample_dir,ref_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get config values
        threads = cfg.get_threads("samtools")
        sortMemory = cfg.get("tools","samtools","sortMemory")
        genome = cfg.get_path("reference","genome_fasta",base_path=ref_dir)
        file_type = cfg.get("tools","samtools","file_type")

        # build output file
        out_file = temp_dir / f"{name}_Aligned_Sorted.{file_type}"

        # build command
        cmd = [
            "samtools",
            "sort",
            "-@", str(threads),
            "-m", str(sortMemory),
            "-o", str(out_file),
        ]

        # append ref if we are using cram
        if file_type == "cram":
            cmd.extend(["-T",str(genome)])

        # add input file
        cmd.append(str(file))

        # execute command
        result = subprocess(cmd, capture_output=True, text=True)

        # log results
        log_subprocess(result, sample_dir, "samtools_sort")

        # get sorted QC
        self.flagstat(out_file,"sorted")

        # return path to sorted file
        return out_file

    def index_file(self, file: Path):
        """
        Indexes given sorted bam file and produces a bam_file.bam.bai comapinon file in same directory as bam_file, second step of samtools processing
        Params:
            file                              Path to sorted bam/cram file to be indexed (we also filter before indexing by default but filtering is not required)
        """
        # get file name
        name = file.stem.split("_Aligned")[0]

        # load config
        cfg = ConfigLoader(self.config)

        # get config dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = project / name

        # make sure directories exist
        for dir in [project,sample_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get values from config
        threads = cfg.get_threads("samtools")
        file_type = cfg.get("tools","samtools","file_type")

        # name output file
        out_file = sample_dir / f"{name}.{file_type}"

        # build command
        cmd = [
            "samtools",
            "index",
            "-@", str(threads),
            "-o", str(out_file),
            str(file)
        ]

        # run command
        result = subprocess(cmd, capture_output=True, text=True)

        # log results
        log_subprocess(result, sample_dir, "samtools_index")

    def filter_file(self, file: Path):
        """
        Filters input bam file to remove maps of poor quality
        Params:
            file                            Path to bam file you want to filter, must be aligned and ideally sorted as well
        Returns:
            Path to the sorted filtered bam file (in the temp dir, we save .CRAM to keep later)
        """
        # get sample name
        name = file.stem.split("_Aligned")[0]

        # load config
        cfg = ConfigLoader(self.config)

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = project / name
        temp_dir = self.temp_dir / name
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)

        # make sure directories exist
        for dir in [project,sample_dir,temp_dir,ref_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get values from config
        threads = cfg.get_threads("samtools")
        minMapQuality = cfg.get("tools","samtools","minMapQuality")
        file_type = cfg.get("tools","samtools","outType")
        filter_f = cfg.get("tools","samtools","filter_f")
        filter_F = cfg.get("tools", "samtools", "filter_F")
        genome = cfg.get_path("reference","genome_fasta",base_path=ref_dir)

        # get file type
        if str(file_type) == "cram":
            outType = "-C"
            ext = ".cram"
        elif str(file_type) == "bam":
            outType = "-b"
            ext = ".bam"
        else:
            outType = ""
            ext = ".sam"

        # build output file
        out_file = temp_dir / f"{name}Aligned_Sorted_Filtered{ext}"

        # build base command
        cmd  = [
            "samtools",
            "view",
            str(outType),
            "-@", str(threads),
            "-q", str(minMapQuality),
            "-o", str(out_file)
        ]

        # append additional filters if specified
        if filter_f:
            cmd.extend(["-f",filter_f])
        if filter_F:
            cmd.extend(["-F", filter_F])

        # handle cram vs bam
        if file_type == "cram":
            cmd.extend(["-T",genome])

        # append input file
        cmd.append(str(file))
        
        # run subprocess
        result = subprocess(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result,sample_dir,"samtools_filter")

        # get filtered QC
        self.flagstat(out_file, "filtered")

        # return path to filtered file
        return out_file

    def flagstat(self, bam_file: Path, file_status: str):
        """
        runs flagstat on given BAM file, done before and after sorting, then after filtering
        Params:
            bam_file                        path to the bam file to run flagstat on
            file_status                     step in the samtools process the flagstat is run on, raw, sorted, or filtered
        """
        # get name of sample
        name = bam_file.stem.split("_Aligned")[0]

        # connect to config
        cfg = ConfigLoader(self.config)

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = project / name

        # make sure directories exist
        for dir in [project,sample_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # build command
        cmd = [
            "samtools",
            "flagstat",
            str(bam_file)
        ]

        # run command
        result = subprocess(cmd,capture_output=True,text=True)
        
        # log subprocess
        log_subprocess(result,sample_dir,f"flagstat_{file_status}")

    def cram_file(self, bam_file: Path):
        """
        compressed BAM files to CRAM files, these need the reference.fasta file to decode so don't lose it
        Params:
            bam_file                        path to bam file to compress
        """
        # name of the file
        name = bam_file.stem.split("_Aligned")[0]

        # load config
        cfg = ConfigLoader(self.config)

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = project / name
        # make sure dirs exist
        for dir in [project,sample_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get ref .fasta file
        ref_fasta = cfg.get_path("reference","genome_fasta",base_path=self.root)

        # build out file
        out_file = sample_dir / f"{name}_Aligned.cram"

        # build command
        cmd = [
            "samtools",
            "view",
            "-C",
            "-T", str(ref_fasta),
            "-o", str(out_file),
            str(bam_file)
        ]

        # run command
        result = subprocess(cmd,capture_output=True,text=True)

        # log subprocess
        log_subprocess(result,sample_dir,"compress bam to cram")
