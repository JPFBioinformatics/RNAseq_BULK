# loads config.yaml and allows easy retrevial of commonly used values

from pathlib import Path
import yaml
import os

class ConfigLoader:
    def __init__(self, config_file: Path):
        """
        Loads yaml file nd stores it as a dictionary as self.config
        params:
            config_file:            Path to config file, should be a relative path and it depends on where you run the script from in this project folder
        """
        self.config_path = Path(config_file)

        if not self.config_path.exists():
            raise FileNotFoundError(f"Config file at {config_file} not found")
        
        with open(self.config_path, "r") as f:
            self.config = yaml.safe_load(f)

    def get(self, *keys: str, default=None):
        """
        accesses nested values from config dict
        params:
            keys:                   list of keys in order of accessing for config structure
                example:            cfg.get("params","star","threads") returns number of threads tasked to star package
        """

        value = self.config
        
        # iterates over all keys given going into each key subsection at each iteartion
        for key in keys:
            # ensures key is a dict
            if not isinstance(value,dict):
                print("ConfigLoader.get() stopped early")   # for debugging, remove when ready for production
                return default
            # resets value (dict) to the value under key, if key is a subdict name it reaturns that dict
            value = value.get(key, default)
        
        # return value queried for
        return value
    
    def get_threads(self, tool_name: str):
        """
        returns the number of threads used for a sepcifeid tool
        overridden by SLURM if running on HPC cluster
        params:
            tool_name:             string name of the tool, such as "star" or "fastqc"
        """
        # threads listed in yaml
        threads = self.get("params", tool_name, "threads", defalut=1)
        # if SLURM specifies a number of threads override at runtime
        threads = int(os.environ.get("SLURM_CPUS_PER_TASK", threads))

        return threads
    
    def get_path(self, *keys: str, base_path: Path = None, must_exist=False):
        """
        Returns Path object for the path specified in the config within the specified base_dir
        Params:
            keys:                   list of keys to go through to find the desired path
                example:            cfg.get_path("data_dirs", "raw") to get the path to the raw file data directory
            base_path:              Path object to concatenate the path we are retreiving, if we want to put subdir within another dir then we would make the dir the base_path
            must_exist:             default to False for directories that the script will create (such as those in data_dirs), true if it needs to exist before (like raw dir)
        """
        p = Path(self.get(*keys))

        if base_path and isinstance(base_path,Path):
            p = base_path / p

        if must_exist and not p.exists():
            raise FileNotFoundError(f"Path {p} not found for keys {keys}")
        
        return p

