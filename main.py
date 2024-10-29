import glob
import os

from modules import Filter
from modules import AnnotationHuman, AnnotationMouse, AnnotationOtherSpecie
from modules import GeneMapping
from modules import SpeciesDataProcessor


def data_filtering_pipeline(animal_dir, species_name, output_dir):
    # Step: Gene Data Filtering
    filter_ref_path = os.path.join(os.getcwd(), "gene_data", "cell_ref")
    filter_instance = Filter(ref_path=filter_ref_path)
    for file in glob.glob(animal_dir):
        filter_instance(data=file, specie=species_name,
                        output_dir=output_dir)  # Assuming filter_data is the method to start filtering


def data_annotation_pipeline(filtered_data_path, species_name, output_dir):
    # Step: Annotation
    annotation_path = "path_to_annotation_model"
    python_module_path = os.getcwd()
    homologous_gene_dir = None
    if species_name == "human":
        annotation_instance = AnnotationHuman(annotation_path=annotation_path, python_module_path=python_module_path)
    elif species_name == "mouse":
        annotation_instance = AnnotationMouse(python_module_path=python_module_path)
    else:
        homologous_gene_dir = os.path.join(os.getcwd(), "gene_data", "homologous_gene")
        annotation_instance = AnnotationOtherSpecie(annotation_path=annotation_path,
                                                    python_module_path=python_module_path)
    for file in glob.glob(filtered_data_path):
        annotation_instance(data=file, specie=species_name,
                            output_dir=output_dir,
                            homologous_gene_dir=homologous_gene_dir)  # Assuming annotate_data is the method to start annotating


def gene_mapping_pipeline(annotated_data_path, species_name, output_dir):
    # Step: Gene Mapping
    mapping_ref_path = os.path.join(os.getcwd(), "gene_data", "cell_ref")
    mapping_instance = GeneMapping(ref_path=mapping_ref_path, output_dir=output_dir)
    for file in glob.glob(annotated_data_path):
        mapping_instance(data=file, specie=species_name)  # Assuming map_data is the method to start mapping


def gene_merging_pipeline(filter_dir, mapping_dir, species_name, output_dir, metadata_path):
    # Step: Gene Merging
    merge_instance = SpeciesDataProcessor(species_name, filter_dir, mapping_dir, output_dir, metadata_path)

    merge_instance()  # Assuming process is the method to start merging


def main(animal_dir, species_name):
    try:
        filtered_data_path = "path_to_filtered_data"
        annotated_data_path = "path_to_annotated_data"
        mapping_dir = "path_to_mapping_dir"
        metadata_path = "path_to_metadata"
        filtered_data_output_dir = "path_to_filtered_data_output"
        annotated_data_output_dir = "path_to_annotated_data_output"
        mapping_output_dir = "path_to_mapping_output"
        merge_output_dir = "path_to_merge_output"

        data_filtering_pipeline(animal_dir, species_name, filtered_data_output_dir)
        data_annotation_pipeline(filtered_data_path, species_name, annotated_data_output_dir)
        gene_mapping_pipeline(annotated_data_path, species_name, mapping_output_dir)
        gene_merging_pipeline(filtered_data_path, mapping_dir, species_name, merge_output_dir, metadata_path)
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    for animal in list():
        main(animal_dir=f"path_to_{animal}_data", species_name=animal)
