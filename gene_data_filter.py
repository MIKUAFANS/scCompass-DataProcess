import glob
import sys, os
import argparse
# import dboxmr
import subprocess
import os
import json
import time
import subprocess
import pickle
import numpy as np
import scanpy as sc
import pandas as pd;
from scipy.stats import zscore


class Filter(object):
    def transform(self, afile, **kwargs):
        print("start")

        specie = kwargs.get('specie')

        latin_maps = {"human": "Homo_sapiens", "mouse": "Mus_musculus", "monkey": "Macaca_mulatta", "zhu": "Sus_scrofa",
                      "dashu": "Rattus_norvegicus", "mianyang": "Ovis_aries", "ji": "Gallus_gallus",
                      "ma": "Equus_caballus",
                      "guoying": "Drosophila_melanogaster", "banmayu": "Danio_rerio", "yang": "Capra_hircus",
                      "xianchong": "Caenorhabditis_elegans", "niu": "Bos_taurus", "quan": "Canis_lupus_familiarisgsd"}

        min_cells = 3
        min_features = 200
        min_protein_mt_number = 6
        percent_mt_max = 15

        ref_path = kwargs.get("ref_path")

        def load_gz_data(in_path):
            """
            加载数据
            """
            adata = sc.read_mtx(os.path.join(in_path, 'matrix.mtx.gz'))
            adata = adata.T

            adata_bc = pd.read_csv(os.path.join(in_path, 'barcodes.tsv.gz'), header=None)
            adata_features = pd.read_csv(os.path.join(in_path, 'features.tsv.gz'), header=None)
            # print(adata_features)
            adata_features = adata_features[0].str.split('\t').str[1].tolist()
            # adata_features = [name.upper() for name in adata_features] # 不能直接全部变成大写

            adata.obs['cell_id'] = list(adata_bc[0])
            adata.var['gene_name'] = adata_features
            adata.var.index = adata.var['gene_name']

            adata.obs['n_counts'] = adata.X.sum(axis=1).A1
            adata.var['n_cells'] = (adata.X > 0).sum(axis=0).A1

            # nan
            adata.X = adata.X.A
            adata.X[np.isnan(adata.X)] = 0
            return adata

        def load_gene_name_id_map_dict(base_dir, species_name):
            """
            从文件列表中读取指定物种 全部基因名列表
            """
            path = os.path.join(base_dir, "Multispecies_all_gene_list", species_name + "_all_genelist.xlsx")

            res_dict = dict()
            df = pd.read_excel(path)
            df = df[~pd.isna(df["Gene name"])]
            res_dict = dict(zip(df["Gene name"], df["Gene stable ID"]))

            return res_dict

        def notin_genelist_filter(adata, gene_list):
            """
            删除在基因名列表中的列
            """
            mask = adata.var.index.isin(gene_list)
            adata = adata[:, mask]
            return adata

        def load_special_gene_list(base_dir, species_name):
            """
            从文件列表中读取指定物种线粒体基因、蛋白质基因、miRNA基因的基因名
            """
            # 读入的内容为 基因原始名称
            protein_file_path = os.path.join(base_dir, "protein_coding_gene_list", species_name + "_protein_coding.txt")
            miRNA_file_path = os.path.join(base_dir, "miRNA_gene_list", species_name + "_miRNA.txt")
            mt_file_path = os.path.join(base_dir, "MT_gene_list", species_name + "_MT.xlsx")

            protein_list, miRNA_list, mt_list = list(), list(), list()

            # protein name
            with open(protein_file_path, "r") as f:
                for line in f:
                    protein_list.append(line.split()[1])

            # miRNA
            with open(miRNA_file_path, "r") as f:
                for line in f:
                    miRNA_list.append(line.split()[1])

            # mitochondria
            df = pd.read_excel(mt_file_path)
            mt_list = df.iloc[:, 2].to_list()

            return protein_list, miRNA_list, mt_list

        def MT_filter_2(adata, percent_mt_max, mt_list):
            """
            过滤 MT占比小于15% 的行
            """
            mito_genes = adata.var_names.str.startswith('MT-')  # 线粒体基因以 MT- 开头
            adata.obs['mito_total'] = adata[:, mito_genes].X.sum(axis=1)
            adata.obs['percent_mito'] = adata.obs['mito_total'] / adata.X.sum(axis=1)
            adata = adata[adata.obs['percent_mito'] < percent_mt_max, :]
            return adata

        def MT_filter_1(adata, gene_list, min_protein_mt_number):
            """
            过滤 蛋白质或miRNA的基因数少于7个的行
            """
            indices = [k in gene_list for k in adata.var.index]
            f_adata = adata[:, indices]
            data = f_adata.X.toarray()
            mask = np.count_nonzero(data, axis=1) > min_protein_mt_number
            return adata[mask]

        def MT_filter_3(adata, mito_list):
            """
            过滤基因总量 在所有细胞平均值三个标准差之外 的行
            """
            total_counts = adata.X.sum(axis=1)
            total_counts = np.array(total_counts).squeeze()
            idx = total_counts > 0
            adata = adata[idx, :]
            total_counts = total_counts[idx]
            gene_name = adata.var['gene_name'].tolist()
            index = [element in mito_list for element in gene_name]

            mito_adata = adata[:, index]
            # 计算每个细胞的线粒体基因表达之和
            mito_counts = mito_adata.X.sum(axis=1)
            mito_counts = np.array(mito_counts).squeeze()
            mito_percentage = mito_counts / total_counts

            # 将每个变量的值进行标准化，计算z-score
            total_counts_zscore = zscore(total_counts)
            mito_percentage_zscore = zscore(mito_percentage)

            # 筛选出符合条件的细胞
            keep_cells = ((total_counts_zscore > -3) & (total_counts_zscore < 3) & (mito_percentage_zscore > -3) & (
                    mito_percentage_zscore < 3))

            # 仅保留符合条件的细胞
            adata = adata[keep_cells, :]
            return adata

        def write_logs(out_path, step_num, cell_num, success):
            write_logs_path = os.path.join(out_path, "logs.txt")

            outtxt = ""
            if cell_num != "-1":
                outtxt = outtxt + "step: " + step_num + ": " + cell_num + "\n"
            if success:
                outtxt = outtxt + "success\n"

            with open(write_logs_path, 'a') as f:
                f.write(outtxt)

        print("start", afile, "size", os.path.getsize(afile))

        # Save the file path output directory
        tmpfile = afile.replace("input_data_R/{}/".format(specie), "control/data_filter_1/{}/".format(specie)) + ".tmp"
        outfile_dir = afile.replace("input_data_R/{}/".format(specie), "control/data_filter_1/{}/".format(specie))
        outfile_dir = outfile_dir + "/"

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

        # with open(tmpfile, 'w') as f:
        #     f.write(afile)

        ########################################## load data ##########################################
        adata = load_gz_data(afile)
        print("原始行、列数: ", adata.shape[0], adata.shape[1])
        gene_name_id_map_dict = load_gene_name_id_map_dict(ref_path, latin_maps[specie])
        protein_list, miRNA_list, MT_list = load_special_gene_list(ref_path, latin_maps[specie])
        write_logs(outfile_dir, "0", str(adata.shape[0]), False)
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

        ########################################## 过滤 表达基因数小于 200 的行 ##########################################
        sc.pp.filter_cells(adata, min_counts=min_features)
        print("过滤表达基因数小于200， 行、列数: ", adata.shape[0], adata.shape[1])
        write_logs(outfile_dir, "1", str(adata.shape[0]), False)
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

        ########################################## 过滤 细胞数小于3 的行 ##########################################
        sc.pp.filter_genes(adata, min_cells=min_cells)
        print("过滤细胞数小于3，行、列数: ", adata.shape[0], adata.shape[1])
        write_logs(outfile_dir, "2", str(adata.shape[0]), False)
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

            # ########################################## 删除 不在基因列表内的 列 ##########################################
        adata = notin_genelist_filter(adata, list(gene_name_id_map_dict.keys()))
        print("删除不在基因列表内，行、列数: ", adata.shape[0], adata.shape[1])
        write_logs(outfile_dir, "3", str(adata.shape[0]), False)
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

            ########################################## 过滤 蛋白质 或 RNA 数目少于7个的行 ##########################################
        adata = MT_filter_1(adata, protein_list + miRNA_list, min_protein_mt_number)
        print("过滤蛋白质或miRNA的基因数少于7个，行、列数: ", adata.shape[0], adata.shape[1])
        write_logs(outfile_dir, "4", str(adata.shape[0]), False)
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

        ########################################## 过滤 MT占比小于15% 的行 ##########################################
        adata = MT_filter_2(adata, percent_mt_max, MT_list)
        print("过滤MT占比小于15% 行、列数: ", adata.shape[0], adata.shape[1])
        write_logs(outfile_dir, "5", str(adata.shape[0]), False)
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

        # ########################################## 过滤基因总量在所有细胞平均值三个标准差之外 的行 ##########################################
        adata = MT_filter_3(adata, MT_list)
        print("过滤基因总量在所有细胞平均值三个标准差之外，行、列数: ", adata.shape[0], adata.shape[1])
        write_logs(outfile_dir, "6", str(adata.shape[0]), False)
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

        print("开始导出数据")
        df = adata.to_df()
        df.to_csv(os.path.join(outfile_dir, outfile_dir[-11: -1] + ".csv"), index=True)
        print("数据已导出")

        if os.path.exists(tmpfile):
            os.unlink(tmpfile)

        write_logs(outfile_dir, "7", "-1", True)

    '''
    入口函数
    '''

    def __call__(self, data, **kwargs):
        return self.transform(data, **kwargs)


if __name__ == "__main__":
    m_map = Filter()

    animal = ""
    # input_dir
    animal_dir = ""
    # reference path
    ref_path = ""

    for file in glob.glob(animal_dir):
        m_map.transform(afile=file,
                        specie=animal,
                        ref_path=ref_path)
