class ChiSquareTest:

    @staticmethod
    def chi2_table(data_set, label, a, b, individuals, i=0):
        """
        创建计算卡方值的表格，2行9列，即分别计算在case 和control 中 SNPi,SNPj各个组合的频数
        """
        case_row = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        control_row = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        list1 = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
        while True:
            j = 0
            for x in list1:
                if data_set[i][a] + data_set[i][b] == str(x[0]) + str(x[1]) and int(label[i]) == 1:
                    case_row[j] += 1
                    break
                elif data_set[i][a] + data_set[i][b] == str(x[0]) + str(x[1]) and int(label[i]) == 0:
                    control_row[j] += 1
                    break
                j += 1
            i += 1
            if i > individuals - 1:
                return case_row, control_row

    @staticmethod
    def chi2_test(case_row, control_row, individuals, i=0):
        """
        计算出列和行的total，计算卡方统计量
        """
        sum_colum = [0] * 9
        A2_NcNr = [0] * 18
        while True:
            sum_colum[i] = float(case_row[i]) + float(control_row[i])  # 将chi2_table列求和，存入sum_list
            i += 1
            if i > 8:
                break
        # 以下两个while计算了chi2_table中每个元素对应的A2_NcNr，并计算了chi2统计量的值
        # 计算case那一行每个元素的A2_NcNr
        i = 0
        while True:
            if case_row[i] == 0:
                A2_NcNr[i] = 0
            else:
                A2_NcNr[i] = case_row[i] * case_row[i] / (sum(case_row) * sum_colum[i])
            i += 1
            if i > 8:
                break
        # 计算control那一行每个元素的A2_NcNr
        i = 0
        while True:
            if case_row[i] == 0:
                A2_NcNr[i + 9] = 0
            else:
                A2_NcNr[i + 9] = control_row[i] * control_row[i] / (sum(control_row) * sum_colum[i])
            i += 1
            if i > 8:
                break
        chi2_statistics = individuals * (sum(A2_NcNr) - 1)
        # print chi_statistics
        return sum_colum, A2_NcNr, chi2_statistics
