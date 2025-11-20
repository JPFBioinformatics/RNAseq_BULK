import subprocess
import sys
from pathlib import Path

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.config_loader import ConfigLoader

"""
This script downlaods reference genome fa and gtf file for mice (mm39) from UCSC genome browser and unzips them for use
places these files where you specify in out, by default in the reference
"""

# load config
cfg = ConfigLoader(root_dir / "config.yaml")
# location of ref dir
ref_dir = cfg.get_path("reference","ref_dir",base_path=root_dir)
ref_dir.mkdir(parents=True,exist_ok=True)

# get base names for files
fa_base = cfg.get("reference","genome_fasta")
gtf_base = cfg.get("reference","gtf_file")

# get urls
fa_url = cfg.get("reference","genome_url")
gtf_url = cfg.get("reference","gtf_url")

# generate names for fa and gtf out files
fa_name = f"{fa_base}.gz"
gtf_name = f"{gtf_base}.gz"

# out files for fa ang gtf
fa_out = ref_dir / fa_name
gtf_out = ref_dir / gtf_name

# downlaod files
for url,out  in [(fa_url,fa_out), (gtf_url,gtf_out)]:

    # downlaod zipped file
    subprocess.run(["wget","-O",str(out),str(url)],check=True)

    # unzip file
    subprocess.run(["gunzip",str(out)],check=True)