import glob
import sys, os
import argparse
import subprocess
import os
import sys
import json
import time
import subprocess
import pickle
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.stats import zscore
import anndata as ad
import datetime


class GeneMapping(object):
    def transform(self, afile, **kwargs):
        # sys.path.append('/mnt/cstr/celldata/paper/codes/zhangran/')
        # 载入文章所提供的参考数据集
        START_TIME = datetime.datetime.now()
        specie = kwargs.get('specie')
        output_dir = kwargs.get('output_dir')
        count_file = afile.split("/")[-1]

        input_file = "{}/{}.csv".format(afile, count_file)
        outfile_dir = "{}{}/{}".format(output_dir, specie, count_file)
        tmpfile = "{}{}/{}.tmp".format(output_dir, specie, count_file)
        mapping_file = "{}{}_mapping.txt".format(output_dir, specie)

        print("input_file ", input_file)
        print("outfile_dir ", outfile_dir)
        print("tmpfile ", tmpfile)

        ref_path = "/mnt/cstr/celldata/input_data/cell_ref/"
        latin_maps = {"human": "Homo_sapiens", "mouse": "Mus_musculus", "monkey": "Macaca_mulatta", "zhu": "Sus_scrofa",
                      "dashu": "Rattus_norvegicus", "mianyang": "Ovis_aries", "ji": "Gallus_gallus",
                      "ma": "Equus_caballus",
                      "guoying": "Drosophila_melanogaster", "banmayu": "Danio_rerio",
                      "xianchong": "Caenorhabditis_elegans", "niu": "Bos_taurus", "quan": "Canis_lupus_familiarisgsd"}

        def load_sorted_core_gene_id_list(base_dir, species_name):
            '''加载核心基因,并排序
            '''
            path = os.path.join(base_dir, "filter_gene_list", species_name + "_ensemble_filter_genelist.txt")

            core_gene_id_list = []
            with open(path, 'r') as f:
                for line in f:
                    gene_id = line.split()[1]
                    core_gene_id_list.append(gene_id)

            return sorted(core_gene_id_list)

        def write_logs(out_path, step_num, cell_num, success):
            write_logs_path = os.path.join(out_path, "logs.txt")

            outtxt = ""
            if cell_num != "-1":
                outtxt = outtxt + "step: " + step_num + ": " + cell_num + "\n"
            if success:
                outtxt = outtxt + "success\n"

            with open(write_logs_path, 'a') as f:
                f.write(outtxt)

        if os.path.exists(afile):
            print("Annotation Start, [Specie]:", specie, "\n[Input Dir]:", afile, "[Size]:", os.path.getsize(afile),
                  '\n')

        ######  Determine whether the input file is existed
        if os.path.exists(input_file):
            print("[Input File]:" + input_file)
        else:
            return

        #### Determine whether the annotation process is completed
        if os.path.exists(outfile_dir):
            print("ignore, done: " + afile)
            if os.path.exists(tmpfile):
                os.unlink(tmpfile)
            return
        else:
            os.makedirs(outfile_dir)

        if os.path.exists(tmpfile):
            print("ignore, busy: " + afile)
            return

        if not os.path.exists(outfile_dir):
            os.makedirs(outfile_dir)

        with open(tmpfile, 'w') as f:
            f.write(afile)

        ########################################## annotation process ##########################################
        # try:
        # # 读取count.csv和cell_type.csv文件
        ### get core gene list
        if os.path.exists(mapping_file):
            core_gene_id_list = list(np.loadtxt(mapping_file, str, delimiter=","))
            # print("mapping file exists")
        else:
            print("core gene list not exists!")
            core_gene_id_list = load_sorted_core_gene_id_list(ref_path, latin_maps[specie])
            np.savetxt(mapping_file, np.array(core_gene_id_list, str), fmt='%s', delimiter=",")

        data_pd = pd.read_csv(input_file)
        print("原始行、列数: ", data_pd.shape[0], data_pd.shape[1])

        ### given a specific file, get the mapping index for core gene list
        positions = [core_gene_id_list.index(element) if element in core_gene_id_list else None for element in
                     data_pd.columns]
        target_columns, columns_to_copy = zip(
            *[(value, index) for index, value in enumerate(positions) if value is not None])

        ### get the data to copy
        temp_array_to_copy = data_pd.iloc[:, list(columns_to_copy)].values

        ### init output matrix
        target_matrix = np.zeros((len(data_pd), len(core_gene_id_list)), int)
        target_matrix[:, target_columns] = temp_array_to_copy

        # # # 将处理后的数据重新写入csv文件
        # ########
        END_TIME = datetime.datetime.now()
        print("开始导出数据")
        print("Time Cost total: {}".format((END_TIME - START_TIME).seconds))
        np.savetxt(os.path.join(outfile_dir, count_file + ".csv"), target_matrix, fmt='%d', delimiter=",")
        print("数据已导出, 数据行数:{},列数:{}".format(target_matrix.shape[0], target_matrix.shape[1]))
        # matrix_df.to_csv(os.path.join(outfile_dir, count_file + ".csv"), index=False, header=False)
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)
        write_logs(outfile_dir, "7", "-1", True)

    '''
    入口函数
    '''

    def __call__(self, data, **kwargs):
        return self.transform(data, **kwargs)


if __name__ == "__main__":
    m_map = GeneMapping()

    animal = ""
    # input_dir
    animal_dir = ""
    output_dir = ""

    for file in glob.glob(animal_dir):
        m_map.transform(afile=file,
                        output_dir=output_dir,
                        specie=animal)
