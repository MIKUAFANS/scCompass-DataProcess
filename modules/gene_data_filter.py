import os
import numpy as np
import scanpy as sc
import pandas as pd;
from scipy.stats import zscore


class Filter:
    def __init__(self, ref_path=None, min_cells=3, min_features=200, min_protein_mt_number=6, percent_mt_max=15):
        """
        Initialize the Filter class with default parameters.
        """
        self.ref_path = ref_path
        self.min_cells = min_cells
        self.min_features = min_features
        self.min_protein_mt_number = min_protein_mt_number
        self.percent_mt_max = percent_mt_max

        self.latin_maps = {
            "human": "Homo_sapiens", "mouse": "Mus_musculus", "monkey": "Macaca_mulatta", "zhu": "Sus_scrofa",
            "dashu": "Rattus_norvegicus", "mianyang": "Ovis_aries", "ji": "Gallus_gallus", "ma": "Equus_caballus",
            "guoying": "Drosophila_melanogaster", "banmayu": "Danio_rerio", "yang": "Capra_hircus",
            "xianchong": "Caenorhabditis_elegans", "niu": "Bos_taurus", "quan": "Canis_lupus_familiarisgsd"
        }

    @staticmethod
    def load_gz_data(in_path):
        """
        Load data from gzipped files.
        """
        adata = sc.read_mtx(os.path.join(in_path, 'matrix.mtx.gz')).T
        adata.obs['cell_id'] = pd.read_csv(os.path.join(in_path, 'barcodes.tsv.gz'), header=None)[0].tolist()
        adata_features = pd.read_csv(os.path.join(in_path, 'features.tsv.gz'), header=None)[0].str.split('\t').str[1]
        adata.var['gene_name'] = adata_features
        adata.var.index = adata.var['gene_name']

        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
        adata.var['n_cells'] = (adata.X > 0).sum(axis=0).A1
        adata.X = adata.X.A
        adata.X[np.isnan(adata.X)] = 0
        return adata

    def load_gene_name_id_map_dict(self, species_name):
        """
        Load the gene name-ID map for a specific species.
        """
        path = os.path.join(self.ref_path, "Multispecies_all_gene_list", f"{species_name}_all_genelist.xlsx")
        df = pd.read_excel(path).dropna(subset=["Gene name"])
        return dict(zip(df["Gene name"], df["Gene stable ID"]))

    def load_special_gene_list(self, species_name):
        """
        Load special gene lists for mitochondria, proteins, and miRNA genes.
        """
        base = self.ref_path
        protein_path = os.path.join(base, "protein_coding_gene_list", f"{species_name}_protein_coding.txt")
        miRNA_path = os.path.join(base, "miRNA_gene_list", f"{species_name}_miRNA.txt")
        mt_path = os.path.join(base, "MT_gene_list", f"{species_name}_MT.xlsx")

        protein_list = [line.split()[1] for line in open(protein_path)]
        miRNA_list = [line.split()[1] for line in open(miRNA_path)]
        mt_list = pd.read_excel(mt_path).iloc[:, 2].tolist()

        return protein_list, miRNA_list, mt_list

    def filter_cells(self, adata, gene_list):
        """
        Filter cells based on various criteria.
        """
        sc.pp.filter_cells(adata, min_counts=self.min_features)
        sc.pp.filter_genes(adata, min_cells=self.min_cells)
        mask = adata.var.index.isin(gene_list)
        return adata[:, mask]

    def filter_mt_percentage(self, adata, mt_list):
        """
        Filter cells based on mitochondrial gene percentage.
        """
        adata.obs['mito_total'] = adata[:, adata.var_names.str.startswith('MT-')].X.sum(axis=1)
        adata.obs['percent_mito'] = adata.obs['mito_total'] / adata.X.sum(axis=1)
        return adata[adata.obs['percent_mito'] < self.percent_mt_max]

    @staticmethod
    def filter_gene_variance(adata, mt_list):
        """
        Filter cells with gene expression outside 3 standard deviations from the mean.
        """
        total_counts = adata.X.sum(axis=1).squeeze()
        valid = total_counts > 0
        adata = adata[valid, :]

        mito_counts = adata[:, adata.var.index.isin(mt_list)].X.sum(axis=1).squeeze()
        mito_percentage = mito_counts / total_counts

        keep_cells = (
                (zscore(total_counts) > -3) & (zscore(total_counts) < 3) &
                (zscore(mito_percentage) > -3) & (zscore(mito_percentage) < 3)
        )
        return adata[keep_cells, :]

    @staticmethod
    def write_logs(path, step, cells, success):
        """
        Write logs to the specified path.
        """
        with open(os.path.join(path, "logs.txt"), 'a') as f:
            f.write(f"step: {step}: {cells}\n{'success\n' if success else ''}")

    def transform(self, afile, **kwargs):
        """
        Main function to filter data.
        """
        specie = kwargs.get('specie')
        species_name = self.latin_maps.get(specie)
        output_dir = kwargs.get('output_dir', os.path.join(os.getcwd(), "filtered_data"))
        output_dir = os.path.join(output_dir, specie)
        tmpfile = os.path.join(output_dir, f"{os.path.basename(afile).split('.')[0]}.tmp")

        if os.path.exists(output_dir):
            print(f"Ignored, already processed: {afile}")
            if os.path.exists(tmpfile):
                os.unlink(tmpfile)
            return

        os.makedirs(output_dir, exist_ok=True)
        adata = self.load_gz_data(afile)
        gene_map = self.load_gene_name_id_map_dict(species_name)
        protein_list, miRNA_list, mt_list = self.load_special_gene_list(species_name)

        adata = self.filter_cells(adata, gene_map.keys())
        adata = self.filter_mt_percentage(adata, mt_list)
        adata = self.filter_gene_variance(adata, mt_list)

        df = adata.to_df()
        df.to_csv(os.path.join(output_dir, f"{specie}_filtered.csv"))
        self.write_logs(output_dir, "7", str(adata.shape[0]), True)
        print(f"Data exported for {afile}")

    def __call__(self, data, **kwargs):
        """
        Entry function to execute the transform method.
        """
        return self.transform(data, **kwargs)

# if __name__ == "__main__":
#     m_map = Filter()
#
#     animal = ""
#     # input_dir
#     animal_dir = ""
#     # reference path
#     ref_path = ""
#
#     for file in glob.glob(animal_dir):
#         m_map.transform(afile=file,
#                         specie=animal,
#                         ref_path=ref_path)
