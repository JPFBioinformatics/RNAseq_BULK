# region Imports

from pathlib import Path
import argparse,shutil,sys

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.config_loader import ConfigLoader
from src.utils import generate_paired_samples, get_scratch, find_name
from src.fastp_wrapper import QCTrimmer
from src.star_wrapper import STARIndexBuilder, STARAligner
from src.samtools_wrapper import SamtoolsWrapper
from src.featureCounts_wrapper import FeatureCountsWrapper

# endregion

def parse_args():
    """
    Accepts CLI arguments passed after caling main.py
    ex: python main.py --config /path/to/config.yaml --sample Sample1 --steps trim
        will run the "trim" section of main.py on sample 1
        allows for easier parallelization on HPC cluster
    """
    parser = argparse.ArgumentParser(description="Rna-seq Bulk")

    parser.add_argument(
        "--root",
        required=True,
        help="Path to root directory of pipeline, typically /home/user/RNAseq_bulk"
    )

    parser.add_argument(
        "--indir",
        required=True,
        help="Path to directory where input files are stored"
    )

    parser.add_argument(
        "--sample1",
        required=False,
        help="Forward read sample name exactly as seen in --indir (with extensions)"
    )

    parser.add_argument(
        "--sample2",
        required=False,
        help="Reverse read sample name exactly as seen in --indir (with extensions)"
    )

    parser.add_argument(
        "--steps",
        required=False,
        nargs="+",
        default=["trim","align","count"],
        help="Piepeline steps to run (trim, align, or count)"
    )

    return parser.parse_args()

def main():

    """
    Runs the pipeline on given samples
    """

    # ---------------------------------------------------
    # Parse CLI arguments and load config
    # ---------------------------------------------------

    # parse args
    args = parse_args()
    # load config
    config = Path(args.root) / "config.yaml"
    cfg =  ConfigLoader(Path(config))


    # ---------------------------------------------------
    # Set up directories
    # ---------------------------------------------------

    # get run name
    run_name = cfg.get("project","name")

    # create directory to store outputs from this pipeline run
    root_dir = Path(args.root)
    run_dir = root_dir / run_name
    run_dir.mkdir(parents=True,exist_ok=True)

    # create reference dir
    ref_dir = run_dir / "reference"
    ref_dir.mkdir(parents=True,exist_ok=True)

    # get temp dir
    temp_dir = get_scratch(root_dir)

    # ---------------------------------------------------
    # Get sample(s) to process
    # ---------------------------------------------------

    # path to dir containing raw data
    in_path = Path(args.indir)

    # generate paired samples list
    paired_samples = generate_paired_samples(in_path, args.sample1, args.sample2)
    
    # ---------------------------------------------------
    # run pipeline
    # ---------------------------------------------------

    for r1,r2 in paired_samples:
        
        # get name of the sample and its sample dir
        name = find_name(r1,r2)
        sample_dir = run_dir / name

        # --------------------------
        # fastp
        # --------------------------

        if "trim" in args.steps:
            # generate QC object
            qc = QCTrimmer(cfg,root_dir,temp_dir)
            # run QC and trim
            trimmed_r1,trimmed_r2 = qc.run_fastp(r1,r2)

        # --------------------------
        # STAR
        # --------------------------
        
        if "align" in args.steps:
            
            # check if index is already built
            star_idx_dir = cfg.get_path("reference","star_index",base_path = ref_dir)
            # if not then build index
            if not star_idx_dir.exists() or not any(star_idx_dir.iterdir()):
                index_builder = STARIndexBuilder(cfg,root_dir,temp_dir)
                index_builder.build_index()

            # instantiate aligner
            aligner = STARAligner(cfg, root_dir, temp_dir)
            # input trimmed/untrimmed data based on user specifications and align
            if "trim" in args.steps:
                aligned_file = aligner.align(trimmed_r1, trimmed_r2, cleanup=True)
            else:
                aligned_file = aligner.align(r1,r2)

        # --------------------------
        # samtools
        # --------------------------
        # need an aligned bam file for further processing
        if "align" in args.steps:

            # isntantiate wrapper
            st = SamtoolsWrapper(cfg, root_dir, temp_dir)

            # sort and filter
            sorted_file = st.sort_file(aligned_file)
            clean_file = st.filter_file(sorted_file)

            # see if we are saving as a cram or bam
            save_type = cfg.get("project","save_type")

            # save file/index to sample dir if specified
            if cfg.get("project","save_files"):

                # if save type is cram then cram the file and save cram/crai
                if save_type == "cram":

                    # cram and index
                    cram_file = st.cram_file(clean_file)
                    cram_idx = st.index_file(cram_file,cram=True)

                    # get genome.fasta file path
                    genome = cfg.get_path("reference","genome_fasta",base_path=ref_dir)
                    # build new path to save copy at
                    genome_path = run_dir / genome.name

                    # copy to new location
                    shutil.copy(cram_file, sample_dir / cram_file.name)
                    shutil.copy(cram_idx, sample_dir / cram_idx.name)
                    shutil.copy(genome, genome_path)

                    # delete old files
                    try:
                        cram_file.unlink()
                        cram_idx.unlink()
                    except Exception as e:
                        print(f"Warning: could not delete origonal file and index after samtools cleaning\n{e}")
                
                # if not cram then save bam/bai
                else:

                    # index bam file
                    bam_idx = st.index_file(clean_file)
                    # save loation of new clean file
                    new_clean = sample_dir / clean_file.name

                    # copy bam/bai to new location
                    shutil.copy(clean_file, new_clean)
                    shutil.copy(bam_idx, sample_dir / bam_idx.name)
                    
                    # delete old files
                    try:
                        clean_file.unlink()
                        bam_idx.unlink()
                    except Exception as e:
                        print(f"Warning: could not delete origonal file and index after samtools cleaning\n{e}")

                    # update clean_file location
                    clean_file = new_clean


        # --------------------------
        # featureCounts
        # --------------------------
        if "count" in args.steps:
            # instantiate fc object
            fc = FeatureCountsWrapper(cfg,root_dir)
            # count features
            fc.count_features(clean_file)



if __name__ == "__main__":
    main()