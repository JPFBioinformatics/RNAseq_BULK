# region Imports

from pathlib import Path
import yaml, os

# endregion

class ConfigLoader:
    """
    loads config.yaml and gives easy access to useful information
    """
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

            # see if key is not in config

            # ensures key is a dict
            if not isinstance(value,dict):
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
        threads = self.get("params", tool_name, "threads", default=1)
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

    def check_bools(self):

        # empty list to hold improperly formatted bools
        errors = []

        # list of fields that must be boolean
        bool_fields = {
            "save_files",
            "specify_adapter",
            "outReadsUnmapped",
            "chromosome",
            "isPairedEnd",
            "largestOverlap",
            "ignoreDup",
            "countFraction"
        }

        # recursive function to enter parent dicts
        def recurse(value, path=""):
            # check each key value pair
            for k,v in value.items():
                # get string representation of current value being observed
                current_path = f"{path}.{k}" if path else k
                # if value is a dict then go another layer deeper
                if isinstance(v,dict):
                    recurse(v,current_path)
                # if not a dict, check if it is in bool_fields and if it is properly formatted as bool
                else:
                    # if imrpoperly formatted then append it to errors list
                    if k in bool_fields and not isinstance(v,bool):
                        errors.append(current_path)
        
        # run recursive method on loaded config dict
        recurse(self.config)
        
        # output error/all good message
        if errors:
            raise ValueError(
                f"Invalid boolean fields found in config.yaml, reformat to True/False:\n"+
                "\n".join(f" - {e}" for e in errors)
                )
        else:
            print("All boolean fields valid, continuing pipeline")

