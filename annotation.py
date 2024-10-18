import glob
import sys, os
import argparse
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
from scimilarity.utils import lognorm_counts
from scimilarity import CellAnnotation, align_dataset


class Annotation(object):
    def transform(self, afile, **kwargs):
        # sys.path.append('/mnt/fc34/pyjobs/paper2/annotation/')
        sys.path.append('/mnt/cstr/celldata/control/')
        print("scimilarity imported!")
        # 载入文章所提供的参考数据集
        annotation_path = '/mnt/cstr/celldata/paper/codes/zhangran/annotation_model_v1'
        # annotation_path = '/mnt/fc34/pyjobs/paper2/annotation/annotation_model_v1'
        ca = CellAnnotation(model_path=annotation_path)
        START_TIME = datetime.datetime.now()
        specie = kwargs.get('specie')
        # input_dir = kwargs.get('input_dir')
        output_dir = kwargs.get('output_dir')
        count_file = afile.split("/")[-1]

        input_file = "{}/{}.csv".format(afile, count_file)
        outfile_dir = "{}{}/{}".format(output_dir, specie, count_file)
        tmpfile = "{}{}/{}.tmp".format(output_dir, specie, count_file)

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
        print("Output Dir: ", outfile_dir)
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
        # 读取数据
        df = pd.read_csv(input_file, header=0, index_col=0)

        # 创建一个Anndata对象
        adata = ad.AnnData(df)
        print("原始行、列数: ", adata.shape[0], adata.shape[1])
        if adata.shape[0] == 0:
            write_logs(outfile_dir, "7", "-1", True)
            return

            # 与参考基因集进行匹配排序
        adata = align_dataset(adata, ca.gene_order)

        # 将原始计数矩阵存储到layer中
        adata.layers["counts"] = adata.X.copy()

        # 标准化
        adata = lognorm_counts(adata)

        # 降维
        adata.obsm['X_scimilarity'] = ca.get_embeddings(adata.X)
        sc.pp.neighbors(adata, use_rep='X_scimilarity')
        sc.tl.umap(adata)

        # 注释
        predictions, nn_idxs, nn_dists, nn_stats = ca.get_predictions_kNN(adata.obsm['X_scimilarity'])
        type_res = pd.DataFrame(predictions.values)

        END_TIME = datetime.datetime.now()
        print("Time Cost: {}".format((END_TIME - START_TIME).seconds))
        print("开始导出数据")
        type_res.to_csv(os.path.join(outfile_dir, count_file + "_cell_type.csv"), index=True)
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
    m_map = Annotation()

    animal = ""
    # input_dir
    animal_dir = ""
    # output_dir
    output_dir = ""

    for file in glob.glob(animal_dir):
        m_map.transform(afile=file,
                        output_dir=output_dir,
                        specie=animal)
