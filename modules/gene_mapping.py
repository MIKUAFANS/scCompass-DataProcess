import os
import numpy as np
import pandas as pd
import datetime


class GeneMapping:
    def __init__(self, ref_path, output_dir):
        """
        Initialize the class with reference and output paths.
        """
        self.ref_path = ref_path
        self.output_dir = output_dir

        self.latin_maps = {
            "human": "Homo_sapiens", "mouse": "Mus_musculus", "monkey": "Macaca_mulatta", "zhu": "Sus_scrofa",
            "dashu": "Rattus_norvegicus", "mianyang": "Ovis_aries", "ji": "Gallus_gallus", "ma": "Equus_caballus",
            "guoying": "Drosophila_melanogaster", "banmayu": "Danio_rerio", "xianchong": "Caenorhabditis_elegans",
            "niu": "Bos_taurus", "quan": "Canis_lupus_familiarisgsd"
        }

    def load_sorted_core_gene_id_list(self, species_name):
        """
        Load and sort the core gene list for the given species.
        """
        path = os.path.join(self.ref_path, "filter_gene_list", f"{species_name}_ensemble_filter_genelist.txt")
        core_gene_id_list = []

        with open(path, 'r') as f:
            for line in f:
                gene_id = line.split()[1]
                core_gene_id_list.append(gene_id)

        return sorted(core_gene_id_list)

    @staticmethod
    def write_logs(out_path, step_num, cell_num, success):
        """
        Write logs to the specified path.
        """
        log_path = os.path.join(out_path, "logs.txt")
        log_entry = f"step: {step_num}: {cell_num}\n" if cell_num != "-1" else ""
        log_entry += "success\n" if success else ""

        with open(log_path, 'a') as f:
            f.write(log_entry)

    def transform(self, afile, **kwargs):
        """
        Main function to map gene data.
        """
        specie = kwargs.get('specie')
        output_dir = kwargs.get('output_dir', self.output_dir)
        count_file = os.path.basename(afile)

        input_file = f"{afile}/{count_file}.csv"
        outfile_dir = f"{output_dir}/{specie}/{count_file}/"
        tmpfile = f"{output_dir}/{specie}/{count_file}.tmp"
        mapping_file = f"{output_dir}/{specie}_mapping.txt"

        print(f"Input file: {input_file}")
        print(f"Output directory: {outfile_dir}")
        print(f"Temporary file: {tmpfile}")

        # Check if the input file exists
        if not os.path.exists(input_file):
            print(f"[Error] Input file does not exist: {input_file}")
            return

        # Check if the annotation process is already completed
        if os.path.exists(outfile_dir):
            print(f"Ignored, already processed: {afile}")
            if os.path.exists(tmpfile):
                os.unlink(tmpfile)
            return

        os.makedirs(outfile_dir, exist_ok=True)

        if os.path.exists(tmpfile):
            print(f"Ignored, already in progress: {afile}")
            return

        with open(tmpfile, 'w') as f:
            f.write(afile)

        # Start processing
        START_TIME = datetime.datetime.now()
        if os.path.exists(mapping_file):
            core_gene_id_list = list(np.loadtxt(mapping_file, dtype=str, delimiter=","))
        else:
            print("Core gene list not found, loading...")
            species_name = self.latin_maps.get(specie)
            core_gene_id_list = self.load_sorted_core_gene_id_list(species_name)
            np.savetxt(mapping_file, np.array(core_gene_id_list, dtype=str), fmt='%s', delimiter=",")

        data_pd = pd.read_csv(input_file)
        print(f"Original rows and columns: {data_pd.shape[0]}, {data_pd.shape[1]}")

        # Map core gene indices
        positions = [
            core_gene_id_list.index(col) if col in core_gene_id_list else None
            for col in data_pd.columns
        ]
        target_columns, columns_to_copy = zip(
            *[(pos, idx) for idx, pos in enumerate(positions) if pos is not None]
        )

        # Copy relevant data
        temp_array_to_copy = data_pd.iloc[:, list(columns_to_copy)].values

        # Initialize output matrix
        target_matrix = np.zeros((len(data_pd), len(core_gene_id_list)), int)
        target_matrix[:, target_columns] = temp_array_to_copy

        # Export the processed data
        END_TIME = datetime.datetime.now()
        print(f"Data export started. Time cost: {(END_TIME - START_TIME).seconds} seconds")
        np.savetxt(
            os.path.join(outfile_dir, f"{count_file}.csv"),
            target_matrix, fmt='%d', delimiter=","
        )
        print(f"Data exported: Rows: {target_matrix.shape[0]}, Columns: {target_matrix.shape[1]}")

        if os.path.exists(tmpfile):
            os.unlink(tmpfile)

        self.write_logs(outfile_dir, "7", "-1", True)

    def __call__(self, data, **kwargs):
        """
        Entry function to execute the transform method.
        """
        return self.transform(data, **kwargs)

# if __name__ == "__main__":
#     m_map = GeneMapping()
#
#     animal = ""
#     # input_dir
#     animal_dir = ""
#     output_dir = ""
#
#     for file in glob.glob(animal_dir):
#         m_map.transform(afile=file,
#                         output_dir=output_dir,
#                         specie=animal)
