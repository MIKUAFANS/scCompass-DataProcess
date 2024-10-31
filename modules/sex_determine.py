import os
from abc import ABC, abstractmethod

import pandas as pd


class SexDetermine(ABC):
    def __init__(self, x_gene_list, y_gene_list):
        self.x_gene_list = x_gene_list
        self.y_gene_list = y_gene_list

    def data_process(self, sex_determine_data):
        x_gene_sum = sex_determine_data.loc[:, sex_determine_data.columns == "Xist"].sum(axis=1)
        total_sum = sex_determine_data.sum(axis=1)
        x_ratio = x_gene_sum / total_sum

        a = sex_determine_data.loc[:, sex_determine_data.columns.isin(self.x_gene_list["Gene name"])]
        b = sex_determine_data.loc[:, sex_determine_data.columns.isin(self.y_gene_list["Gene name"])]
        c = b.loc[:, ~b.columns.isin(a.columns)]
        y_gene_sum = sex_determine_data.loc[:, sex_determine_data.columns.isin(c.columns)].sum(axis=1)
        total_sum1 = sex_determine_data.sum(axis=1)
        y_ratio = y_gene_sum / total_sum1

        sample_y_gene = pd.DataFrame(y_ratio, columns=['ratioY'])
        sample_x_gene = pd.DataFrame(x_ratio, columns=['ratioX'])

        return sample_x_gene, sample_y_gene

    @abstractmethod
    def determine(self, sex_determine_data):
        pass


class HumanSexDetermine(SexDetermine):
    def __init__(self, x_gene_list=None, y_gene_list=None):
        if x_gene_list is None:
            x_gene_list = pd.read_excel(
                os.path.join(os.getcwd(), "gene_data", "sex_determine_data", "Homo_sapiens_chrX.xlsx"))
        if y_gene_list is None:
            y_gene_list = pd.read_csv(
                os.path.join(os.getcwd(), "gene_data", "sex_determine_data", "human_Ygenelist.txt"), sep=",")
        super().__init__(x_gene_list, y_gene_list)

    def determine(self, sex_determine_data):
        sample_x_gene, sample_y_gene = self.data_process(sex_determine_data)
        gender_determine_result = []

        for idx in range(sample_x_gene.shape[0]):
            if (sample_y_gene.loc[idx, 'ratioY'] >= 0.001188).all():
                gender = "Male"
            elif (sample_x_gene.loc[idx, 'ratioX'] <= 0.000070).all() and (sample_y_gene.loc[idx, 'ratioY'] > 0.000106).all():
                gender = "Male"
            elif (sample_x_gene.loc[idx, 'ratioX'] >= 0.000001).all() and (sample_y_gene.loc[idx, 'ratioY'] <= 0.000106).all():
                gender = "Female"
            elif (sample_x_gene.loc[idx, 'ratioX'] > 0.000070).all() and (0.000106 < sample_y_gene.loc[idx, 'ratioY']).all() and (
                    sample_y_gene.loc[idx, 'ratioY'] <= 0.001188).all():
                gender = "Mix"
            else:
                gender = "unknown"

            gender_determine_result.append(gender)

        return gender_determine_result


class MouseSexDetermine(SexDetermine):
    def __init__(self, x_gene_list=None, y_gene_list=None):
        if x_gene_list is None:
            x_gene_list = pd.read_excel(
                os.path.join(os.getcwd(), "gene_data", "sex_determine_data", "Mus_musculus_chrX.xlsx"))
        if y_gene_list is None:
            y_gene_list = pd.read_csv(
                os.path.join(os.getcwd(), "gene_data", "sex_determine_data", "mouse_Ygenelist.txt"), sep=",")
        super().__init__(x_gene_list, y_gene_list)

    def determine(self, sex_determine_data):
        sample_x_gene, sample_y_gene = self.data_process(sex_determine_data)
        gender_determine_result = []

        for idx in range(sample_x_gene.shape[0]):
            if (sample_y_gene.loc[idx, 'ratioY'] >= 0.000065).all():
                gender = "Male"
            elif (sample_x_gene.loc[idx, 'ratioX'] <= 0.000008).all() and (sample_y_gene.loc[idx, 'ratioY'] > 0.000004).all() and (
                    sample_y_gene.loc[idx, 'ratioY'] < 0.000065).all():
                gender = "Male"
            elif (sample_x_gene.loc[idx, 'ratioX'] >= 0.000555).all() and (sample_y_gene.loc[idx, 'ratioY'] > 0.000004).all() and (
                    sample_y_gene.loc[idx, 'ratioY'] < 0.000065).all():
                gender = "Female"
            elif (sample_x_gene.loc[idx, 'ratioX'] > 0).all() and (sample_y_gene.loc[idx, 'ratioY'] <= 0.000004).all():
                gender = "Female"
            elif (0.000008 < sample_x_gene.loc[idx, 'ratioX']).all() and (sample_x_gene.loc[idx, 'ratioX'] < 0.000555).all() and (
                    0.000004 < sample_y_gene.loc[idx, 'ratioY']).all() and (sample_y_gene.loc[idx, 'ratioY'] < 0.000065).all():
                gender = "Mix"
            else:
                gender = "unknown"

            gender_determine_result.append(gender)

        return gender_determine_result
