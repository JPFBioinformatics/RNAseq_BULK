# region Imports

from pathlib import Path
import argparse,shutil,sys,json

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
from src.counts import Counts

# endregion

cfg = ConfigLoader(root_dir / "config.yaml")
counter = Counts(root_dir, cfg)

counter.summarize_counts()

matrix = counter.preprocess_full()

counter.pca(matrix, 3)

counter.plot_pca()