import os
import collections


class Prepare:

    @staticmethod
    def snps_sequence1(start_snp, end_snp):
        snps_sequence = []
        for i in range(start_snp, end_snp + 1):
            snps_sequence.append(i)
        return snps_sequence

    @staticmethod
    def scaling_rank(n, number):
        summary = 0
        l = []
        for i in range(number):
            temp = (1 - summary) * n
            l.append(temp)
            summary += temp
        return l

    @staticmethod
    def dict_pair(start_snp, end_snp):
        """
        dict_pair = {
            0: [(0, 1), (0, 2), (0, 3)..., (0, end_snp)],
            1: [(1, 0), (1, 2), (1, 3)..., (1, end_snp)],
            2: [(2, 0), (2, 1), (2, 3)..., (2, end_snp)],
            ...
            end_snp: [(end_snp, 0),(end_snp, 1),..., (end_snp, end_snp-1)]
        }
        """
        dict_pair = collections.defaultdict(list)
        for i in range(start_snp, end_snp + 1):
            for j in range(start_snp, i):
                dict_pair[i].append((i, j))
            for j in range(i + 1, end_snp + 1):
                dict_pair[i].append((i, j))
        return dict_pair

    @staticmethod
    def dict_prob1(start_snp, end_snp):
        """
        dict_prob = {
            0: [1, 1, 1, ..., 1, 1, 1],
            1: [1, 1, 1, ..., 1, 1, 1],
            2: [1, 1, 1, ..., 1, 1, 1]
            ...
            end_snp: [1, 1, 1, ..., 1, 1, 1]
        }
        """
        dict_prob = collections.defaultdict(list)
        if len(dict_prob) != 0:
            print(dict_prob)
        for i in range(start_snp, end_snp):
            for x in range(end_snp - start_snp - 1):
                dict_prob[i].append(1)
        return dict_prob

    @staticmethod
    def data_matrix(file_order, individuals, snps_number):
        # 原始数据集转存data_set矩阵，及最后一列class存入label
        if file_order < 10:
            file = open(os.getcwd() + '\data_set\\2SNP_MAF0.4_MAF0.6_H0.1_models_D{}_S1000_1000_EDM-5_00'.format(
                snps_number) + str(file_order) + '.txt', 'r')
        elif file_order < 100:
            file = open(os.getcwd() + '\data_set\\2SNP_MAF0.4_MAF0.6_H0.1_models_D{}_S1000_1000_EDM-5_0'.format(
                snps_number) + str(file_order) + '.txt', 'r')
        else:
            file = open(os.getcwd() + '\data_set\\2SNP_MAF0.4_MAF0.6_H0.1_models_D{}_S1000_1000_EDM-5_'.format(
                snps_number) + str(file_order) + '.txt', 'r')
        data_list = file.readlines()[1:individuals + 1]  # [1:2001],   [1,n]意思是从第二行开始，一共读出（N-1）行
        label = []
        data_set = []
        for line in data_list:
            line = line.split('\t')
            line[snps_number].replace('\n', '')
            label.append(line[snps_number])
            del line[snps_number]
            data_set.append(line)
        file.close()
        return label, data_set

    @staticmethod
    def load_weka_ranking(weka_ranking_file_order):
        """
        Read weka snps ranking, example:
        dict_weka = {
            100: 1.0143
            99: 1.0109
            10: 1.01065
            37: 1.01045
            ...
        }
        Normalize the ReliefF weights so that they lie between 0 and 2
        """
        file = open(os.getcwd() + '\weka_ranking\\' + str(weka_ranking_file_order) + '.txt', 'r')
        list2 = file.readlines()
        dict_weka = {}
        for line in list2:
            x = line.split()
            dict_weka[int(x[1]) - 1] = float(
                x[0]) + 1
        file.close()
        return dict_weka





