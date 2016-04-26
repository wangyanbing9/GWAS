import os
from prepare import Prepare
from ant import Ant
from chi_square_test import ChiSquareTest


def gwas(start_snp, end_snp, individuals, ants_number, update_times, factor, threshold):
    """
     individuals: sum of cases and control.
     factor:  0<factor<1 user-adjustable distribution.if factor is near 1,snps with high scores are only slightly
              more likely to be selected.
     threshold <-> P-Value: 26.12<->0.001；20.09<->0.01；15.51<->.0.05
    """
    snps_number = end_snp - start_snp + 1
    dict_pair = Prepare.dict_pair(start_snp, end_snp)
    print('dict_pair done')

    n = 0
    for order in range(1, 101):
        dict_prob = Prepare.dict_prob1(start_snp, end_snp + 1)
        print('dict_prob done')
        label, data_set = Prepare.data_matrix(order, individuals, snps_number)
        dict_weka = Prepare.load_weka_ranking(order)
        snps_sequence = Prepare.snps_sequence1(start_snp, end_snp)
        snps_scores = [0] * snps_number
        # Rank dict on weight,return a list of tuples
        list_weka = sorted(dict_weka.items(), key=lambda weka: weka[1], reverse=True)
        rank = Prepare.scaling_rank(factor, snps_number)
        for snp_weight_tuple, new_weight in zip(list_weka, rank):
            # allocate scaling rank weight to corresponding snp
            snps_scores[snp_weight_tuple[0] - start_snp] = new_weight
        selected_pairs = {}
        for x in range(update_times):
            list_k1_threshold = []  # 存储chi2_statistics > 阈值的pair在dict_prob[first_snp]中的索引位置
            list_k2_threshold = []  # 存储chi2_statistics > 阈值的pair在dict_prob[second_snp]中的索引位置
            list_chi2_statistics = []  # 存储chi2_statistics > 阈值的pair的卡方值，在update的时候使用
            list_k1 = []  # 存储chi2_statistics < 阈值的pair在dict_prob[first_snp]中的索引位置
            list_k2 = []  # 存储chi2_statistics < 阈值的pair在dict_prob[second_snp]中的索引位置
            list_first_snp_threshold = []  # 存储chi2_statistics > 阈值的first_snp
            list_second_snp_threshold = []  # 存储chi2_statistics > 阈值的second_snp
            list_first_snp = []  # 存储chi2_statistics < 阈值的first_snp
            list_second_snp = []  # 存储chi2_statistics < 阈值的second_snp
            for y in range(ants_number):
                first_snp = Ant.ant_first_snp(snps_sequence, snps_scores)
                second_snp = Ant.ant_second_snp(first_snp, dict_pair, dict_prob)
                # print(first_snp, second_snp)
                if second_snp > first_snp:
                    k1 = second_snp - start_snp - 1
                    k2 = first_snp - start_snp
                else:
                    k1 = second_snp - start_snp
                    k2 = first_snp - start_snp - 1
                case_row, control_row = ChiSquareTest.chi2_table(data_set, label, first_snp, second_snp, individuals)
                sum_colum, A2_NcNr, chi2_statistics = ChiSquareTest.chi2_test(case_row, control_row, individuals)
                if chi2_statistics > threshold:
                    list_k1_threshold.append(k1)
                    list_k2_threshold.append(k2)
                    list_chi2_statistics.append(chi2_statistics)
                    list_first_snp_threshold.append(first_snp)
                    list_second_snp_threshold.append(second_snp)
                else:
                    list_k1.append(k1)
                    list_k2.append(k2)
                    list_first_snp.append(first_snp)
                    list_second_snp.append(second_snp)
            # 更新条件概率和snps_scores.
            # 除以 10 的目的是将卡方值转换到和weight同样的量级.snps_score_original[x]引入expert information
            for m, i, j in zip(list_first_snp_threshold, list_k1_threshold, list_chi2_statistics):
                dict_prob[m][i] += (int(j) / 10)
            for m, i, j in zip(list_second_snp_threshold, list_k2_threshold, list_chi2_statistics):
                dict_prob[m][i] += (int(j) / 10)
            for m, i in zip(list_first_snp, list_k1):
                dict_prob[m][i] = 0
            for m, i in zip(list_second_snp, list_k2):
                dict_prob[m][i] = 0
            for first, second, i in zip(list_first_snp_threshold, list_second_snp_threshold, list_k1_threshold):
                selected_pairs[(first, second)] = dict_prob[first][i]
            print(selected_pairs)
        result = sorted(selected_pairs.items(), key=lambda x: x[1], reverse=True)
        for x in result[0:20]:
            if x[0] == (999, 988) or x[0] == (998, 999):
                n += 1
                break
        f = open('result_conProb-H0.1-{ants}_{updates}-fac={factor}.txt'.format(
            ants=ants_number, updates=update_times, factor=factor), 'a')
        f.writelines('\n')
        f.writelines(str(x[0]) + '\t' + str(x[1]) + '\n' for x in result[0:20])
        f.writelines('----------------------------------')
        f.close()
    f = open('result_conProb-H0.1-{ants}_{updates}-fac={factor}.txt'.format(
        ants=ants_number, updates=update_times, factor=factor), 'a')
    f.writelines('\nThe pow of this ACO algorithm is ' + str(n) + '%')
    f.close()
if __name__ == '__main__':
    gwas(0, 99, 2000, 250, 20, 0.1, 15.51)
