import numpy as np
import pandas as pd
import os
import re
import csv
import subprocess

# data filter path
DATA_FILTER_DIR = ""
# data mapping path
DATA_MAPPING_DIR = ""
# save path
DATA_MERGE_DIR = ""

animals = ["human", "mouse", "dashu", "zhu", "monkey", "quan", "mianyang", "niu", "ji", "guoying", "xianchong",
           "banmayu", "ma"]


### log function
def write_logs(out_path, outtxt):
    write_logs_path = os.path.join(out_path, "logs.txt")
    with open(write_logs_path, 'a') as f:
        f.write(outtxt)


def get_sample_organ_array(path, specie_name):
    ### get the dirs
    specie_pd = pd.read_excel(os.path.join(path, specie_name), index_col=0).fillna(0)
    organ_filter_pd = specie_pd[specie_pd['Organ'] != 0]

    return organ_filter_pd


def append_line(text_str, target_file):
    with open(target_file, "a") as file:
        file.write(text_str + "\n")


if __name__ == "__main__":
    for animal in animals:
        specie_name = animal

        sample_organ_array = get_sample_organ_array(
            "metadata_path", animal)

        sample_num = 0
        cell_num = 0
        specie_filter_dir = "{}/{}".format(DATA_FILTER_DIR, specie_name)
        specie_mapping_dir = "{}/{}".format(DATA_MAPPING_DIR, specie_name)
        for index, row in sample_organ_array.iterrows():
            sample_name = index
            organ_name = row["Organ"]
            sample_cell_type_dir = "{}/{}".format(specie_filter_dir, sample_name)
            sample_gene_mapping_dir = "{}/{}".format(specie_mapping_dir, sample_name)
            if os.path.exists(sample_cell_type_dir) and os.path.exists(sample_gene_mapping_dir):
                ### get the file dirs and read the data
                sample_num += 1

                out_str = "sample_num:{} start, sample_cell_type_dir: {}\n".format(sample_num, sample_cell_type_dir)
                out_path = "{}/{}/".format(DATA_MERGE_DIR, specie_name)
                os.makedirs(out_path, exist_ok=True)
                print(out_str)
                write_logs(out_path, out_str)

                cell_type_file = "{}/{}_cell_type.csv".format(sample_cell_type_dir, sample_name)
                cell_mapping_file = "{}/{}.csv".format(sample_gene_mapping_dir, sample_name)
                try:
                    cell_types = pd.read_csv(cell_type_file, header=None)
                    cell_mapping_data = np.loadtxt(cell_mapping_file, delimiter=',', dtype=int)
                except pd.errors.EmptyDataError:
                    with open("empty_data.txt", "a") as f:
                        f.write(f"{cell_type_file}\n")
                    continue

                organ_name = re.sub(r'\s', '', organ_name)

                ###### pandas
                cell_types_data = pd.read_csv(cell_type_file, header=None)
                grouped = cell_types_data.groupby(1)
                for key, group in grouped:
                    print(f"Category: {key}")
                    ### filter the string
                    current_cell_type = re.sub(r'\s', '', key)
                    current_cell_type = current_cell_type.replace(",", "-").replace("#", "-").replace("%",
                                                                                                      "-").replace(
                        "$", "-").replace("^", "-").replace("&", "-").replace("_", "-")
                    target_merge_file = "{}/{}/{}_{}.csv".format(DATA_MERGE_DIR, specie_name, organ_name,
                                                                 current_cell_type)

                    current_cell_type_indexes = np.array(group[0])
                    current_cell_type_data = cell_mapping_data[current_cell_type_indexes, :]

                    # 打开CSV文件以追加模式
                    with open(target_merge_file, 'a', newline='') as file:
                        writer = csv.writer(file, delimiter=',')
                        # 追加多行数据
                        writer.writerows(current_cell_type_data)
        print(sample_num)
