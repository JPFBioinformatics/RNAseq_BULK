# region Imports

import sys,json
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from pathlib import Path

# location of pipeline root dir
root_dir = Path(__file__).resolve().parent.parent
# tell python to look here for modules
sys.path.insert(0, str(root_dir))

from src.config_loader import ConfigLoader

# endregion

class Counts:
    """
    Summariezes counts data and does basic preliminary analysis (PCA)
    """

    def __init__(self, root: Path, cfg: ConfigLoader):

        self.root = Path(root)
        self.cfg = cfg

        name = cfg.get("project","name")
        self.data_dir = Path(root) / name

    @staticmethod
    def parse_count(file: Path):
        """
        Parses a single featureCounts count file and generates a summary dict
        Params:
            file                        Path to the file you want to extract to dict
        Returns:
            a dict that contians [geneID]:count key:value pairs for this file
        """

        # init output dict
        counts = {}
        # header flag
        header = False

        with open(file, "r") as f:
            for line in f:

                # skip over starting comments/header
                if line.startswith("#"):
                    continue
                if not header:
                    header = True
                    continue
                
                # split the line into parts at each tab
                parts = line.strip().split("\t")
                counts[parts[0]] = int(parts[-1])

        return counts

    @staticmethod
    def preprocess_filter(matrix: np.ndarray, threshold: int=10):
        """
        filters out count values lower than an iput threshold, removing them from the matrix
        Params:
            matrix                          2d numpy array to filter
            threshold                       min counts threshold (remove counts lower than this)
        Returns
            filtered_matrix                 matrix with low expressing values removed
        """

        keep_values = (matrix >= threshold).sum(axis=0) >= (0.5 * matrix.shape[0])
        filtered_matrix = matrix[:,keep_values]

        return filtered_matrix
    
    @staticmethod
    def preprocess_cpm(matrix: np.ndarray):
        """
        Normailzes counts matrix so each is represented in cpm (counts per million), so each sample (column) will sum to 1 million this helps to account for diffreences
        in library sizes/number of counts per sample
        Params:
            matrix                      counts matrix you want to normalize, rows = samples cols = genes
        Returns:
            cpm                         normalized matrix with counts now represented in cpm
        """
        # genearte float matrix to store values in
        float_matrix = matrix.astype(float)

        # get the library size of each sample (total counts for that column)
        library_sizes = float_matrix.sum(axis=1,keepdims=True)

        # calculate counts per million by dividing each count by the library size then multiplying by 1 million
        cpm = (float_matrix / library_sizes) * 1e6

        return cpm
    
    @staticmethod
    def preprocess_log2(matrix: np.ndarray):
        """
        Takes log2(cell) for each cell of the matrix to normalize count values, add 1 to each value to avoid doing log2 of 0
        Params:
            matrix                          2d count matrix you want to take the log of
        Returns:
            log2_matrix                     matrix that has been log2 tranformed
        """
        log2_matrix = np.log2(matrix+1)

        return log2_matrix
    
    @staticmethod
    def preprocess_zscore(matrix: np.ndarray):
        """
        Calculates a zscore for each value in the matrix representing its standard deviations from the mean value for each column to correct for column
        speciic variance and prevent columns with generally high values from dominating the analysis
        Params:
            matrix                          2d count matrix you want to generate zscores of
        Returns:
            z_matrix                        matrix of z-score values generated from input matrix
        """
        scaler = StandardScaler()
        z_matrix = scaler.fit_transform(matrix)

        return z_matrix

    def preprocess_full(self):
        """
        Fully preprocesses the counts matrix, filering, converting to CPM, taking log2 correction and finally converting to zscore in that order
        Returns:
            zscore                          fully preprocessed counts matrix
        """
        # get raw matrix
        raw = self.matrix

        # filter, convert to cpm then log2 and finally to zscores
        filtered = Counts.preprocess_filter(raw, 10)
        cpm = Counts.preprocess_cpm(filtered)
        log2 = Counts.preprocess_log2(cpm)
        zscore = Counts.preprocess_zscore(log2)

        # save values
        self.raw_matrix = self.matrix
        self.matrix = zscore

        return zscore
    
    def summarize_counts(self):
        """
        Grabs all count data and summarieze into a metadata dict and counts matrix
        """

        # find data dir
        data_dir = self.data_dir
        
        # generate a list of all counts.txt files from data_dir
        count_files = list(data_dir.rglob("*_counts.txt"))

        # initialize gene/sample maps and lists
        gene_map = {}
        genes = []
        sample_map = {}
        samples = []

        # generate gene map/genes
        first_counts = self.parse_count(count_files[0])
        for idx,(key,_) in enumerate(first_counts.items()):
            gene_map[key] = idx
            genes.append(key)

        # initialize counts matrix (rows = samples columns = genes)
        counts = np.zeros((len(count_files),len(gene_map)))

        # grab counts from each sample run
        for idx,file in enumerate(count_files):

            # get file name and add file to map/list
            file_name = file.name.split("_")[0]
            sample_map[file_name] = idx
            samples.append(file_name)

            # get dict of this file's gene:count
            data = self.parse_count(file)

            # add each count value to the correct position in the counts matrix
            for key,value in data.items():
                col = gene_map[key]
                counts[idx][col] = value

        # generate summary json for the data
        summary = {
            "gene_map": gene_map,
            "genes": genes,
            "sample_map": sample_map,
            "samples": samples
        }

        # save matrix and summary dict to object
        self.matrix = counts
        self.metadata = summary

    def save_counts(self):
        """
        saves the count matrix and count metadata to data directory for this run as npy and json files respectively
        """
        data_dir = self.data_dir

        # save numpy array
        np.save(data_dir/"counts_matrix.npy",self.matrix)

        # save summary json
        output_file = data_dir / "counts_metadata.json"
        with open(output_file, "w") as f:
            json.dump(self.metadata,f,indent=4)

    def pca(self, matrix: np.ndarray, comps: int=2):
        """
        Preforms PCA on the matrix, saves the number of principal compoenents specified by user and saves to object
        Params:
            comps                           Number of principal components you want to keep for this data matrix
            matrix                          The matrix you want to process
        """
        # throw error if too many PCs are asked for
        if comps > min(len(self.metadata["samples"]), len(self.metadata["genes"])):
            raise ValueError(f"Too many PCs requested, comps must be <= min(number_samples, number_genes)")


        # calculate principal components
        pca = PCA(n_components=comps)
        pc_scores = pca.fit_transform(matrix)

        # get variance explaiend by each pc
        variance = pca.explained_variance_ratio_

        # save values to this object
        self.pc_scores = pc_scores
        self.pc_variance = variance

    def plot_pca(self, pcx: int=0, pcy: int=1, save_path: Path=None):
        """
        Geneartes a scatter plot for PCA as a pdf file in the data directory for this run for later use
        Params:
            pcx/pcy                         which principal components you want to graph (0 being most prominant, 1 being next etc...)
            save_path                       optional path to where you want to save the figure if not the default section
        """

        # check to make sure pcs are there
        if self.pc_scores is None:
            raise RuntimeError("No PCA components found, please run Counts.pca() first")
        
        # get your x and y values
        x = self.pc_scores[:,pcx]
        y = self.pc_scores[:,pcy]

        # create initial plot
        plt.figure(figsize=(7, 6))
        plt.scatter(x, y, s=70)

        # get variance for these pcs
        vx = self.pc_variance[pcx] * 100
        vy = self.pc_variance[pcy] * 100

        # additional plot options
        plt.xlabel(f"PC{pcx+1} ({vx:.2f}% variance)")
        plt.ylabel(f"PC{pcy+1} ({vy:.2f}% variance)")
        plt.title(f"PC{pcy+1} vs PC{pcx+1}")
        plt.grid(True)

        # handle saving
        if save_path:
            out_path = save_path
        else:
            out_path = self.data_dir / f"PCA_PC{pcy+1}_vs_PC{pcx+1}"

        plt.savefig(out_path, bbox_inches="tight")

        plt.show
        
        



