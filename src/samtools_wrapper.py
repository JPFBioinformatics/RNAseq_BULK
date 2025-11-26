# region Imports

from pathlib import Path
import sys, subprocess

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.utils import log_subprocess
from src.config_loader import ConfigLoader

# endregion

class SamtoolsWrapper:
    """
    Wrapper class for common samtools operations on BAM files:
        -sort
        -filter
        -index
        -flatstat (QC for logs)
    pipeline will do these operations in above order by default, flagstat happens multiple times not just at the end
    """

    def __init__(self, cfg: ConfigLoader, root: Path, temp_dir: Path, sample_dir: Path):
        """
        Params:
            cfg                             ConfigLoader object that has loaded config.yaml for this project
            root                            path to root folder of project
            temp_dir                        path to the temp dir for intermediate files
            sample_dir                      directory where sample data will be stored
        """
        self.root = Path(root)
        self.cfg = cfg
        self.temp_dir = Path(temp_dir)
        self.sample_dir = Path(sample_dir)

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
        name = file.stem.split("Aligned")[0]

        # load config
        cfg = self.cfg

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = self.sample_dir
        temp_dir = self.temp_dir / name
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)

        # make sur directories exist
        for dir in [project,sample_dir,ref_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get config values
        threads = cfg.get_threads("samtools")
        sortMemory = cfg.get("tools","samtools","sortMemory")
        
        # build output file
        out_file = temp_dir / f"{name}_Aligned_Sorted.bam"

        # build command
        cmd = [
            "samtools",
            "sort",
            "-@", str(threads),
            "-m", str(sortMemory),
            "-o", str(out_file),
            str(file)
        ]

        # execute command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log results
        log_subprocess(result, sample_dir, "samtools_sort")

        # get sorted QC
        self.flagstat(out_file,"sorted")
        
        # if subprocess was successful then delete input file
        if out_file.exists():
            try:
                file.unlink()
                print(f"Samtools sort complete, deleted input file:\n{file}\n")
            except Exception as e:
                print(f"Warning, could not delete input file:\n{file}\nError:\n{e}\n")
        
        # return path to sorted file
        return out_file

    def index_file(self, file: Path, cram=False):
        """
        Indexes given sorted bam file and produces a bam_file.bam.bai comapinon file in same directory as bam_file, second step of samtools processing
        Params:
            file                              Path to sorted bam/cram file to be indexed (we also filter before indexing by default but filtering is not required)
        """
        # get file name
        name = file.stem.split("_Aligned")[0]

        # load config
        cfg = self.cfg

        # get config dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = self.sample_dir
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)

        # make sure directories exist
        for dir in [project,sample_dir,ref_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get values from config
        ext = "bai"
        threads = cfg.get_threads("samtools")
        if cram:
            genome = cfg.get_path("reference","genome_fasta",base_path=ref_dir)
            ext = "crai"
        
        # out file location
        out_file = self.temp_dir / name / f"{name}.{ext}"

        # build command
        cmd = [
            "samtools",
            "index",
            "-@", str(threads),
            "-o", str(out_file),
        ]

        # if we are indexing a cram file then add genome fasta location
        if cram:
            cmd.extend(["-T",genome])

        # append input file
        cmd.append(str(file))

        # run command
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log results
        log_subprocess(result, sample_dir, "samtools_index")

        return out_file

    def filter_file(self, file: Path):
        """
        Filters input bam file to remove maps of poor quality
        Params:
            file                            Path to bam file you want to filter, must be aligned and ideally sorted as well
        Returns:
            Path to the sorted filtered bam file (in the temp dir, we save later if specified)
        """
        # get sample name
        name = file.stem.split("_Aligned")[0]

        # load config
        cfg = self.cfg

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = self.sample_dir
        temp_dir = self.temp_dir / name
        ref_dir = cfg.get_path("reference","ref_dir",base_path=self.root)

        # make sure directories exist
        for dir in [project,sample_dir,temp_dir,ref_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get values from config
        threads = cfg.get_threads("samtools")
        minMapQuality = cfg.get("tools","samtools","minMapQuality")
        filter_f = cfg.get("tools","samtools","filter_f")
        filter_F = cfg.get("tools", "samtools", "filter_F")

        # build output file
        out_file = temp_dir / f"{name}_Aligned_Sorted_Filtered.bam"

        # build base command
        cmd  = [
            "samtools",
            "view",
            "-b",
            "-@", str(threads),
            "-q", str(minMapQuality),
            "-o", str(out_file)
        ]

        # append additional filters if specified
        if filter_f:
            cmd.extend(["-f",filter_f])
        if filter_F:
            cmd.extend(["-F", filter_F])

        # append input file
        cmd.append(str(file))
        
        # run subprocess
        result = subprocess.run(cmd, capture_output=True, text=True)

        # log subprocess
        log_subprocess(result,sample_dir,"samtools_filter")

        # get filtered QC
        self.flagstat(out_file, "filtered")
        
        # if subprocess was successful then delete input file
        if out_file.exists():
            try:
                file.unlink()
                print(f"Samtools filter complete, deleted input file:\n{file}\n")
            except Exception as e:
                print(f"Warning, could not delete input file:\n{file}\nError:\n{e}\n")
        
        return out_file

    def flagstat(self, file: Path, file_status: str):
        """
        runs flagstat on given BAM file, done before and after sorting, then after filtering
        Params:
            file                            path to the file to run flagstat on
            file_status                     step in the samtools process the flagstat is run on, raw, sorted, or filtered
        """

        # connect to config
        cfg = self.cfg

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = self.sample_dir

        # make sure directories exist
        for dir in [project,sample_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # build command
        cmd = [
            "samtools",
            "flagstat",
            str(file)
        ]

        # run command
        result = subprocess.run(cmd,capture_output=True,text=True)
        
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
        cfg = self.cfg

        # get dirs
        project = cfg.get_path("project","name",base_path=self.root)
        sample_dir = self.sample_dir
        temp_dir = self.temp_dir
        # make sure dirs exist
        for dir in [project,sample_dir,temp_dir]:
            dir.mkdir(parents=True,exist_ok=True)

        # get ref .fasta file
        ref_fasta = cfg.get_path("reference","genome_fasta",base_path=self.root)

        # build out file
        out_file = temp_dir / f"{name}.cram"

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
        result = subprocess.run(cmd,capture_output=True,text=True)

        # log subprocess
        log_subprocess(result,sample_dir,"cram")
        
        # if successful then delete input bam file
        if out_file.exists():
            try:
                bam_file.unlink
                print(f"Cram file generated successfully, deleted input bam:\n{bam_file}\n")
            except Exception as e:
                print(f"Warning, could not delete origonal BAM file:\n{bam_file}\nError:\n{e}\n")

        # return output cram file
        return out_file
