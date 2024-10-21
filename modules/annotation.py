import glob
import os
import sys
import scanpy as sc
import pandas as pd
import anndata as ad
import datetime
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset


class Annotation:
    def __init__(self, annotation_path, python_module_path):
        """
        Initialize the class and set required paths.
        """
        # Initialize path variables
        self.annotation_path = annotation_path
        self.cell_annotation = CellAnnotation(model_path=self.annotation_path)

        # Set Python module paths
        sys.path.append(python_module_path)
        print("scimilarity imported!")

    def write_logs(self, out_path, step_num, cell_num, success):
        """
        Write logs to the specified path.
        """
        write_logs_path = os.path.join(out_path, "logs.txt")
        outtxt = f"step: {step_num}: {cell_num}\n" if cell_num != "-1" else ""
        outtxt += "success\n" if success else ""

        with open(write_logs_path, 'a') as f:
            f.write(outtxt)

    def transform(self, afile, **kwargs):
        """
        Main function to perform the annotation task.
        """
        specie = kwargs.get('specie')
        output_dir = kwargs.get('output_dir', os.path.join(os.getcwd(), "annotated_data"))
        count_file = os.path.basename(afile)
        input_file = f"{afile}/{count_file}.csv"
        outfile_dir = f"{output_dir}/{specie}/{count_file}"
        tmpfile = f"{output_dir}/{specie}/{count_file}.tmp"

        # Check if the input directory exists
        if not os.path.exists(afile):
            print(f"[Error] Input path does not exist: {afile}")
            return

        print(f"Annotation Start, [Specie]: {specie}\n[Input Dir]: {afile} [Size]: {os.path.getsize(afile)}\n")

        if not os.path.exists(input_file):
            print(f"[Error] Input file does not exist: {input_file}")
            return

        print(f"Output Dir: {outfile_dir}")

        if os.path.exists(outfile_dir):
            print(f"Ignored, already completed: {afile}")
            if os.path.exists(tmpfile):
                os.unlink(tmpfile)
            return

        if os.path.exists(tmpfile):
            print(f"Ignored, in progress: {afile}")
            return

        os.makedirs(outfile_dir, exist_ok=True)

        with open(tmpfile, 'w') as f:
            f.write(afile)

        # Start the annotation process
        START_TIME = datetime.datetime.now()
        try:
            df = pd.read_csv(input_file, header=0, index_col=0)
            adata = ad.AnnData(df)
            print("Original rows, columns: ", adata.shape[0], adata.shape[1])

            if adata.shape[0] == 0:
                self.write_logs(outfile_dir, "7", "-1", True)
                return

            # Align, normalize, and perform dimensionality reduction
            adata = align_dataset(adata, self.cell_annotation.gene_order)
            adata.layers["counts"] = adata.X.copy()
            adata = lognorm_counts(adata)
            adata.obsm['X_scimilarity'] = self.cell_annotation.get_embeddings(adata.X)

            sc.pp.neighbors(adata, use_rep='X_scimilarity')
            sc.tl.umap(adata)

            # Generate annotation results
            predictions, nn_idxs, nn_dists, nn_stats = self.cell_annotation.get_predictions_kNN(
                adata.obsm['X_scimilarity']
            )
            type_res = pd.DataFrame(predictions.values)

            END_TIME = datetime.datetime.now()
            print(f"Time Cost: {(END_TIME - START_TIME).seconds} seconds")

            # Export results
            type_res.to_csv(os.path.join(outfile_dir, f"{count_file}_cell_type.csv"), index=True)
            print("Data exported")

        finally:
            if os.path.exists(tmpfile):
                os.unlink(tmpfile)
            self.write_logs(outfile_dir, "7", "-1", True)

    def __call__(self, data, **kwargs):
        """
        Entry function to execute the transform method.
        """
        return self.transform(data, **kwargs)


# if __name__ == "__main__":
#     m_map = Annotation()
#
#     animal = ""
#     # input_dir
#     animal_dir = ""
#     # output_dir
#     output_dir = ""
#
#     for file in glob.glob(animal_dir):
#         m_map.transform(afile=file,
#                         output_dir=output_dir,
#                         specie=animal)
